use std::sync::Mutex;

use base64::Engine;
use tauri::ipc::Response;

use pcat_pipeline::cpr::{self, CprFrame};
use pcat_pipeline::stretched_cpr;
use crate::state::AppState;

// ---------------------------------------------------------------------------
// Legacy result types — kept for backward-compatible commands
// ---------------------------------------------------------------------------

#[derive(serde::Serialize)]
pub struct CprCommandResult {
    /// Base64-encoded f32 little-endian bytes of the CPR image
    pub image_base64: String,
    /// [height, width] of the image
    pub shape: [usize; 2],
    /// Arc-length positions in mm for each column
    pub arclengths: Vec<f64>,
}

#[derive(serde::Serialize)]
pub struct CrossSectionCommandResult {
    /// Base64-encoded f32 little-endian bytes of the cross-section image
    pub image_base64: String,
    /// Size of the square image
    pub pixels: usize,
    /// Arc-length position in mm
    pub arc_mm: f64,
    /// Equivalent-circle lumen diameter in mm (FWHM per-ray scan)
    pub vessel_diameter_mm: f64,
    /// Lumen boundary polygon [x, y] in pixel coords
    pub vessel_wall: Vec<[f64; 2]>,
}

// ---------------------------------------------------------------------------
// Helper: clone CprFrame out of AppState for use on a blocking thread
// ---------------------------------------------------------------------------

fn clone_frame(frame_ref: &CprFrame) -> CprFrame {
    CprFrame {
        positions: frame_ref.positions.clone(),
        tangents: frame_ref.tangents.clone(),
        normals: frame_ref.normals.clone(),
        binormals: frame_ref.binormals.clone(),
        arclengths: frame_ref.arclengths.clone(),
    }
}

// ---------------------------------------------------------------------------
// Phase 1: Build frame
// ---------------------------------------------------------------------------

/// Build and cache the CPR frame from a centerline.
/// Called once when the centerline changes.
///
/// - `centerline_mm`: Dense centerline points in [z, y, x] mm.
/// - `pixels_wide`: Number of arc-length samples (output columns).
#[tauri::command]
pub async fn build_cpr_frame(
    centerline_mm: Vec<[f64; 3]>,
    pixels_wide: usize,
    state: tauri::State<'_, Mutex<AppState>>,
) -> Result<(), String> {
    if centerline_mm.len() < 2 {
        return Err("centerline must have at least 2 points".into());
    }
    if pixels_wide < 2 {
        return Err("pixels_wide must be at least 2".into());
    }

    let frame = tokio::task::spawn_blocking(move || {
        CprFrame::from_centerline(&centerline_mm, pixels_wide)
    })
    .await
    .map_err(|e| format!("build_cpr_frame task failed: {e}"))?;

    let mut guard = state.lock().map_err(|e| format!("lock poisoned: {e}"))?;
    guard.cpr_frame = Some(std::sync::Arc::new(frame));

    Ok(())
}

// ---------------------------------------------------------------------------
// Phase 2: Raw binary IPC commands (new, fast)
// ---------------------------------------------------------------------------

/// Render a straightened CPR image. Returns raw binary:
///   [width: u32 LE][height: u32 LE][n_arclengths: u32 LE]
///   [arclengths: n * f64 LE]
///   [image: width*height * f32 LE]
#[tauri::command]
pub async fn render_cpr_image(
    rotation_deg: f64,
    width_mm: f64,
    pixels_high: usize,
    slab_mm: f64,
    state: tauri::State<'_, Mutex<AppState>>,
) -> Result<Response, String> {
    if pixels_high < 2 {
        return Err("pixels_high must be at least 2".into());
    }

    let (volume_data, spacing, origin, direction, frame) = {
        let guard = state.lock().map_err(|e| format!("lock poisoned: {e}"))?;
        let vol = guard.volume.as_ref()
            .ok_or_else(|| "no volume loaded".to_string())?;
        let frame_ref = guard.cpr_frame.as_ref()
            .ok_or_else(|| "no CPR frame built -- call build_cpr_frame first".to_string())?;
        (vol.data.clone(), vol.spacing, vol.origin, vol.direction, clone_frame(frame_ref))
    };

    let result = tokio::task::spawn_blocking(move || {
        frame.render_cpr(&volume_data, spacing, origin, &direction, rotation_deg, width_mm, pixels_high, slab_mm)
    })
    .await
    .map_err(|e| format!("render_cpr_image task failed: {e}"))?;

    // Pack binary: header + arclengths + image
    let n_arc = result.arclengths.len();
    let n_pixels = result.image.len();
    let header_size = 12; // 3 x u32
    let arc_size = n_arc * 8; // f64
    let img_size = n_pixels * 4; // f32
    let mut bytes = Vec::with_capacity(header_size + arc_size + img_size);

    bytes.extend_from_slice(&(result.pixels_wide as u32).to_le_bytes());
    bytes.extend_from_slice(&(result.pixels_high as u32).to_le_bytes());
    bytes.extend_from_slice(&(n_arc as u32).to_le_bytes());
    bytes.extend_from_slice(bytemuck::cast_slice::<f64, u8>(&result.arclengths));
    bytes.extend_from_slice(bytemuck::cast_slice::<f32, u8>(&result.image));

    Ok(Response::new(bytes))
}

/// Render a stretched CPR image. Returns raw binary (same format as straightened):
///   [width: u32 LE][height: u32 LE][n_arclengths: u32 LE]
///   [arclengths: n * f64 LE]
///   [image: width*height * f32 LE]
#[tauri::command]
pub async fn render_stretched_cpr_image(
    rotation_deg: f64,
    width_mm: f64,
    pixels_wide: usize,
    pixels_high: usize,
    slab_mm: f64,
    state: tauri::State<'_, Mutex<AppState>>,
) -> Result<Response, String> {
    if pixels_wide < 2 || pixels_high < 2 {
        return Err("output dimensions must be at least 2".into());
    }

    let (volume_data, spacing, origin, direction, frame) = {
        let guard = state.lock().map_err(|e| format!("lock poisoned: {e}"))?;
        let vol = guard.volume.as_ref()
            .ok_or_else(|| "no volume loaded".to_string())?;
        let frame_ref = guard.cpr_frame.as_ref()
            .ok_or_else(|| "no CPR frame built -- call build_cpr_frame first".to_string())?;
        (vol.data.clone(), vol.spacing, vol.origin, vol.direction, clone_frame(frame_ref))
    };

    // Grow the vertical viewport so the vessel's out-of-plane depth
    // excursion fits. `pixels_high` is the caller's minimum; see
    // `effective_stretched_pixels_high` for the sizing rule.
    let effective_high = {
        let geom = stretched_cpr::compute_stretched_geometry(
            &frame.positions,
            pixels_wide,
            rotation_deg,
        );
        effective_stretched_pixels_high(&geom, pixels_high)
    };

    let result = tokio::task::spawn_blocking(move || {
        frame.render_stretched(
            &volume_data, spacing, origin, &direction,
            rotation_deg, width_mm,
            pixels_wide, effective_high, slab_mm,
        )
    })
    .await
    .map_err(|e| format!("render_stretched_cpr_image task failed: {e}"))?;

    // Same binary format as straightened CPR
    let n_arc = result.arclengths.len();
    let n_pixels = result.image.len();
    let header_size = 12;
    let arc_size = n_arc * 8;
    let img_size = n_pixels * 4;
    let mut bytes = Vec::with_capacity(header_size + arc_size + img_size);

    bytes.extend_from_slice(&(result.pixels_wide as u32).to_le_bytes());
    bytes.extend_from_slice(&(result.pixels_high as u32).to_le_bytes());
    bytes.extend_from_slice(&(n_arc as u32).to_le_bytes());
    bytes.extend_from_slice(bytemuck::cast_slice::<f64, u8>(&result.arclengths));
    bytes.extend_from_slice(bytemuck::cast_slice::<f32, u8>(&result.image));

    Ok(Response::new(bytes))
}

/// Render batch cross-sections. Returns raw binary:
///   [n_sections: u32 LE]
///   For each section:
///     [pixels: u32 LE][arc_mm: f64 LE][diameter_mm: f64 LE][n_wall: u32 LE]
///     [image: pixels*pixels * f32 LE]
///     [wall: n_wall * 2 * f32 LE]  (x, y pairs in pixel coords)
#[tauri::command]
pub async fn render_cross_sections(
    position_fractions: Vec<f64>,
    rotation_deg: f64,
    width_mm: f64,
    pixels: usize,
    state: tauri::State<'_, Mutex<AppState>>,
) -> Result<Response, String> {
    if pixels < 2 {
        return Err("output size must be at least 2".into());
    }
    for &frac in &position_fractions {
        if !(0.0..=1.0).contains(&frac) {
            return Err(format!("position_fraction must be in [0, 1], got {frac}"));
        }
    }

    let (volume_data, spacing, origin, direction, frame) = {
        let guard = state.lock().map_err(|e| format!("lock poisoned: {e}"))?;
        let vol = guard.volume.as_ref()
            .ok_or_else(|| "no volume loaded".to_string())?;
        let frame_ref = guard.cpr_frame.as_ref()
            .ok_or_else(|| "no CPR frame built -- call build_cpr_frame first".to_string())?;
        (vol.data.clone(), vol.spacing, vol.origin, vol.direction, clone_frame(frame_ref))
    };

    let results = tokio::task::spawn_blocking(move || {
        frame.render_cross_sections(
            &volume_data, spacing, origin, &direction,
            &position_fractions, rotation_deg, width_mm, pixels,
        )
    })
    .await
    .map_err(|e| format!("render_cross_sections task failed: {e}"))?;

    // Pack: header + per-section data (see doc comment above for layout).
    let n_sections = results.len();
    let mut bytes = Vec::with_capacity(
        4 + results
            .iter()
            .map(|r| 4 + 8 + 8 + 4 + r.image.len() * 4 + r.vessel_wall.len() * 2 * 4)
            .sum::<usize>(),
    );

    bytes.extend_from_slice(&(n_sections as u32).to_le_bytes());
    for r in &results {
        let n_wall = r.vessel_wall.len() as u32;
        bytes.extend_from_slice(&(r.pixels as u32).to_le_bytes());
        bytes.extend_from_slice(&r.arc_mm.to_le_bytes());
        bytes.extend_from_slice(&r.vessel_diameter_mm.to_le_bytes());
        bytes.extend_from_slice(&n_wall.to_le_bytes());
        bytes.extend_from_slice(bytemuck::cast_slice::<f32, u8>(&r.image));

        // Flatten wall polygon into an [x0, y0, x1, y1, ...] f32 stream.
        let mut wall_flat: Vec<f32> = Vec::with_capacity(r.vessel_wall.len() * 2);
        for [x, y] in &r.vessel_wall {
            wall_flat.push(*x as f32);
            wall_flat.push(*y as f32);
        }
        bytes.extend_from_slice(bytemuck::cast_slice::<f32, u8>(&wall_flat));
    }

    Ok(Response::new(bytes))
}

// ---------------------------------------------------------------------------
// Phase 3: Projection info (lightweight JSON)
// ---------------------------------------------------------------------------

#[derive(serde::Serialize)]
pub struct CprProjectionInfo {
    pub total_arc_mm: f64,
    pub total_proj_arc_mm: f64,
    pub half_width_mm: f64,
    pub projection_normal: [f64; 3],
    pub mid_height_point: [f64; 3],
    pub dy_mm: f64,
    pub pixels_wide: usize,
    pub pixels_high: usize,
    /// Lookup table for `worldToStretchedCpr`. Uniformly sampled in projected
    /// arc-length; sized independently of the render resolution so that the
    /// frontend's per-seed segment search stays cheap.
    pub proj_col_pts: Vec<[f64; 3]>,
    pub arclengths: Vec<f64>,
    pub positions: Vec<[f64; 3]>,
    pub normals: Vec<[f64; 3]>,
}

/// Upper bound on `proj_col_pts` length used by `get_cpr_projection_info`.
/// Keep this small — the frontend projection math is insensitive to the exact
/// length, and the old PCA projection info was effectively O(1). Rendering
/// still uses the full requested `pixels_wide`; only the frontend lookup table
/// is capped.
const PROJECTION_INFO_MAX_COLS: usize = 128;

/// Vertical margin (mm) around the vessel's depth excursion when auto-sizing
/// `pixels_high`. Keep this in one place so renderer + projection-info agree.
const STRETCHED_VERTICAL_MARGIN_MM: f64 = 6.0;
/// Hard cap on auto-grown `pixels_high` to bound IPC / frontend paint cost.
const STRETCHED_MAX_PIXELS_HIGH: usize = 1024;

/// Compute the effective `pixels_high` that the stretched renderer will use,
/// given a caller-supplied minimum. Grows to fit the vessel's out-of-plane
/// depth excursion so oblique rotations no longer clip off the panel.
fn effective_stretched_pixels_high(
    geom: &stretched_cpr::StretchedGeometry,
    min_pixels_high: usize,
) -> usize {
    let depth_span_mm = geom.proj_max - geom.proj_min;
    let needed = ((depth_span_mm + 2.0 * STRETCHED_VERTICAL_MARGIN_MM) / geom.dy_mm).ceil()
        as usize;
    min_pixels_high.max(needed).min(STRETCHED_MAX_PIXELS_HIGH)
}

/// Return the projection parameters needed to map 3D seed positions
/// to/from 2D CPR canvas coordinates.
#[tauri::command]
pub async fn get_cpr_projection_info(
    rotation_deg: f64,
    width_mm: f64,
    pixels_wide: usize,
    pixels_high: usize,
    state: tauri::State<'_, Mutex<AppState>>,
) -> Result<CprProjectionInfo, String> {
    let frame = {
        let guard = state.lock().map_err(|e| format!("lock poisoned: {e}"))?;
        let frame_ref = guard.cpr_frame.as_ref()
            .ok_or_else(|| "no CPR frame built".to_string())?;
        clone_frame(frame_ref)
    };

    if pixels_wide < 2 {
        return Err("pixels_wide must be at least 2".into());
    }

    // Cap the lookup-table resolution for the frontend. The renderer still
    // uses the caller's `pixels_wide` when it renders; this is only the size
    // of the `proj_col_pts` array returned to the frontend for seed/overlay
    // projection, which does not need the full render resolution.
    let lookup_cols = pixels_wide.min(PROJECTION_INFO_MAX_COLS);

    let geom = stretched_cpr::compute_stretched_geometry(
        &frame.positions,
        lookup_cols,
        rotation_deg,
    );

    // `dy_mm` returned to the frontend must match the *render* resolution,
    // not `lookup_cols`. The frontend inverts
    //     row = pixels_high/2 - depth_mm / dy_mm
    // to place seeds vertically; feeding it `geom.dy_mm` (which is
    // total_proj_arc / (lookup_cols - 1)) would over-report pixel spacing by
    // `(pixels_wide - 1) / (lookup_cols - 1)` and visibly pull seed markers
    // away from the rendered vessel.
    let dy_mm_render = geom.total_proj_arc / (pixels_wide - 1) as f64;

    // The renderer may grow `pixels_high` beyond the caller's request to fit
    // the vessel's depth excursion. Report the same effective value here so
    // the frontend's seed/row math lines up with what was actually drawn.
    //
    // `effective_stretched_pixels_high` reads `dy_mm` from the passed geom,
    // but we want the grow-rule keyed to the *render* dy_mm above, not the
    // lookup-table dy_mm. Build a view of the geom with the render dy_mm
    // swapped in so the ceiling division uses the same resolution the
    // renderer will use.
    let mut geom_for_height = geom;
    geom_for_height.dy_mm = dy_mm_render;
    let effective_high = effective_stretched_pixels_high(&geom_for_height, pixels_high);
    let geom = geom_for_height;

    let total_arc = *frame.arclengths.last().unwrap_or(&0.0);

    // Rotated Bishop normals -- still needed for straightened CPR overlays.
    let (rot_normals, _rot_binormals) = frame.rotated_frame(rotation_deg);
    let normals_arr: Vec<[f64; 3]> = rot_normals.iter()
        .map(|n| [n[0], n[1], n[2]])
        .collect();

    let proj_col_pts_arr: Vec<[f64; 3]> = geom.proj_col_pts.iter()
        .map(|v| [v[0], v[1], v[2]])
        .collect();

    Ok(CprProjectionInfo {
        total_arc_mm: total_arc,
        total_proj_arc_mm: geom.total_proj_arc,
        half_width_mm: width_mm,
        projection_normal: [
            geom.projection_normal[0],
            geom.projection_normal[1],
            geom.projection_normal[2],
        ],
        mid_height_point: [
            geom.mid_height_point[0],
            geom.mid_height_point[1],
            geom.mid_height_point[2],
        ],
        dy_mm: dy_mm_render,
        pixels_wide,
        pixels_high: effective_high,
        proj_col_pts: proj_col_pts_arr,
        arclengths: frame.arclengths.clone(),
        positions: frame.positions.clone(),
        normals: normals_arr,
    })
}

// ---------------------------------------------------------------------------
// Legacy commands -- kept for backward compatibility but delegate to new API
// ---------------------------------------------------------------------------

/// Legacy: Compute a CPR image in one call (builds frame + renders).
#[tauri::command]
pub async fn compute_cpr_image(
    centerline_mm: Vec<[f64; 3]>,
    rotation_deg: f64,
    width_mm: f64,
    slab_mm: f64,
    pixels_wide: usize,
    pixels_high: usize,
    state: tauri::State<'_, Mutex<AppState>>,
) -> Result<CprCommandResult, String> {
    if centerline_mm.len() < 2 {
        return Err("centerline must have at least 2 points".into());
    }
    if pixels_wide < 2 || pixels_high < 2 {
        return Err("output dimensions must be at least 2".into());
    }

    let (volume_data, spacing, origin, direction) = {
        let guard = state.lock().map_err(|e| format!("lock poisoned: {e}"))?;
        let vol = guard.volume.as_ref()
            .ok_or_else(|| "no volume loaded".to_string())?;
        (vol.data.clone(), vol.spacing, vol.origin, vol.direction)
    };

    let result = tokio::task::spawn_blocking(move || {
        cpr::compute_cpr(
            &volume_data, &centerline_mm, spacing, origin, &direction,
            width_mm, slab_mm, pixels_wide, pixels_high, rotation_deg,
        )
    })
    .await
    .map_err(|e| format!("CPR task failed: {e}"))?;

    let bytes: &[u8] = bytemuck::cast_slice(&result.image);
    let image_base64 = base64::engine::general_purpose::STANDARD.encode(bytes);

    Ok(CprCommandResult {
        image_base64,
        shape: [result.pixels_high, result.pixels_wide],
        arclengths: result.arclengths,
    })
}

/// Legacy: Compute a single cross-section image.
#[tauri::command]
pub async fn compute_cross_section_image(
    centerline_mm: Vec<[f64; 3]>,
    position_fraction: f64,
    rotation_deg: f64,
    width_mm: f64,
    pixels: usize,
    state: tauri::State<'_, Mutex<AppState>>,
) -> Result<CrossSectionCommandResult, String> {
    if centerline_mm.len() < 2 {
        return Err("centerline must have at least 2 points".into());
    }
    if pixels < 2 {
        return Err("output size must be at least 2".into());
    }
    if !(0.0..=1.0).contains(&position_fraction) {
        return Err(format!(
            "position_fraction must be in [0, 1], got {position_fraction}"
        ));
    }

    let (volume_data, spacing, origin, direction) = {
        let guard = state.lock().map_err(|e| format!("lock poisoned: {e}"))?;
        let vol = guard.volume.as_ref()
            .ok_or_else(|| "no volume loaded".to_string())?;
        (vol.data.clone(), vol.spacing, vol.origin, vol.direction)
    };

    let result = tokio::task::spawn_blocking(move || {
        cpr::compute_cross_section(
            &volume_data, &centerline_mm, spacing, origin, &direction,
            position_fraction, rotation_deg, width_mm, pixels,
        )
    })
    .await
    .map_err(|e| format!("cross-section task failed: {e}"))?;

    let bytes: &[u8] = bytemuck::cast_slice(&result.image);
    let image_base64 = base64::engine::general_purpose::STANDARD.encode(bytes);

    Ok(CrossSectionCommandResult {
        image_base64,
        pixels: result.pixels,
        arc_mm: result.arc_mm,
        vessel_diameter_mm: result.vessel_diameter_mm,
        vessel_wall: result.vessel_wall,
    })
}

/// Legacy: Batch-compute multiple cross-sections.
#[tauri::command]
pub async fn compute_cross_sections_batch(
    centerline_mm: Vec<[f64; 3]>,
    position_fractions: Vec<f64>,
    rotation_deg: f64,
    width_mm: f64,
    pixels: usize,
    state: tauri::State<'_, Mutex<AppState>>,
) -> Result<Vec<CrossSectionCommandResult>, String> {
    if centerline_mm.len() < 2 {
        return Err("centerline must have at least 2 points".into());
    }
    if pixels < 2 {
        return Err("output size must be at least 2".into());
    }
    for &frac in &position_fractions {
        if !(0.0..=1.0).contains(&frac) {
            return Err(format!(
                "position_fraction must be in [0, 1], got {frac}"
            ));
        }
    }

    let (volume_data, spacing, origin, direction) = {
        let guard = state.lock().map_err(|e| format!("lock poisoned: {e}"))?;
        let vol = guard.volume.as_ref()
            .ok_or_else(|| "no volume loaded".to_string())?;
        (vol.data.clone(), vol.spacing, vol.origin, vol.direction)
    };

    let results = tokio::task::spawn_blocking(move || {
        cpr::compute_cross_sections_batch(
            &volume_data, &centerline_mm, spacing, origin, &direction,
            &position_fractions, rotation_deg, width_mm, pixels,
        )
    })
    .await
    .map_err(|e| format!("batch cross-section task failed: {e}"))?;

    Ok(results
        .into_iter()
        .map(|r| {
            let bytes: &[u8] = bytemuck::cast_slice(&r.image);
            let image_base64 = base64::engine::general_purpose::STANDARD.encode(bytes);
            CrossSectionCommandResult {
                image_base64,
                pixels: r.pixels,
                arc_mm: r.arc_mm,
                vessel_diameter_mm: r.vessel_diameter_mm,
                vessel_wall: r.vessel_wall,
            }
        })
        .collect())
}
