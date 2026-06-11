use std::collections::HashMap;
use std::path::Path;
use std::sync::{Arc, Mutex};

use ndarray::Array2;
use serde::{Deserialize, Serialize};
use tauri::{Emitter, Manager};

use pcat_pipeline::active_contour::{
    compute_gradient_field, evolve_snake as pipeline_evolve_snake, init_circular_contour,
    insert_control_point, SnakeParams,
};
use pcat_pipeline::annotation::{self, AnnotationBatchParams, AnnotationTarget};
use pcat_pipeline::cpr::CprFrame;
use pcat_pipeline::mmd::{self, MmdResult, WlAnchor, WlCalibration};
use pcat_pipeline::radial_angular::{self, CrossSectionSurface, RadialAngularParams};
use pcat_pipeline::roi;

use crate::state::AppState;

// ---------------------------------------------------------------------------
// Result types returned to the frontend
// ---------------------------------------------------------------------------

#[derive(Serialize)]
pub struct SnakeResult {
    pub points: Vec<[f64; 2]>,
    pub iterations: usize,
    pub max_displacement: f64,
    pub converged: bool,
}

#[derive(Serialize)]
pub struct MmdSummary {
    pub method: String,
    pub iterations: usize,
    pub converged: bool,
    pub n_voxels: usize,
    pub mean_water_frac: f64,
    pub mean_lipid_frac: f64,
    pub mean_iodine_frac: f64,
    pub mean_calcium_frac: f64,
}

// ---------------------------------------------------------------------------
// Progress event payload
// ---------------------------------------------------------------------------

#[derive(Serialize, Clone)]
struct AnnotationProgress {
    stage: String,
    progress: f64,
}

// ---------------------------------------------------------------------------
// Helper: clone CprFrame for use on blocking thread
// ---------------------------------------------------------------------------

fn clone_frame(frame: &CprFrame) -> CprFrame {
    CprFrame {
        positions: frame.positions.clone(),
        tangents: frame.tangents.clone(),
        normals: frame.normals.clone(),
        binormals: frame.binormals.clone(),
        arclengths: frame.arclengths.clone(),
    }
}

// ---------------------------------------------------------------------------
// Commands
// ---------------------------------------------------------------------------

/// Generate annotation targets for the proximal RCA (or any vessel segment).
///
/// Returns cross-section images, auto-detected vessel walls, and initial snake
/// boundaries for each of N evenly-spaced cross-sections along the centerline.
#[tauri::command]
pub async fn generate_annotation_targets(
    centerline_mm: Vec<[f64; 3]>,
    ostium_mm: Option<[f64; 3]>,
    state: tauri::State<'_, Mutex<AppState>>,
) -> Result<Vec<AnnotationTarget>, String> {
    if centerline_mm.len() < 2 {
        return Err("centerline must have at least 2 points".into());
    }

    // Extract volume data and CprFrame under lock, then release.
    let (volume_data, spacing, origin, direction, frame) = {
        let guard = state.lock().map_err(|e| format!("lock poisoned: {e}"))?;
        let vol = guard
            .volume
            .as_ref()
            .ok_or_else(|| "no volume loaded".to_string())?;
        let frame = guard
            .cpr_frame
            .as_ref()
            .ok_or_else(|| "no CPR frame built — call build_cpr_frame first".to_string())?;
        (vol.data.clone(), vol.spacing, vol.origin, vol.direction, Arc::clone(frame))
    };

    let frame_clone = clone_frame(&frame);

    // Run heavy computation on a blocking thread.
    let targets = tokio::task::spawn_blocking(move || {
        let params = AnnotationBatchParams {
            ostium_zyx: ostium_mm,
            ..AnnotationBatchParams::default()
        };
        annotation::generate_annotation_batch(
            &frame_clone,
            &volume_data,
            spacing,
            origin,
            &direction,
            &params,
        )
    })
    .await
    .map_err(|e| format!("generate_annotation_targets task failed: {e}"))?;

    // Store targets in state.
    {
        let mut guard = state.lock().map_err(|e| format!("lock poisoned: {e}"))?;
        guard.annotation_targets = Some(targets.clone());
        // Reset contours and finalization state for fresh annotation session.
        guard.snake_contours.clear();
        guard.finalized.clear();
    }

    Ok(targets)
}

/// Initialize the snake contour on a specific cross-section.
///
/// If `init_radius_mm` is provided, creates a circular contour at that radius;
/// otherwise uses the target's pre-computed `init_boundary`.
///
/// Returns the initial contour points as `[x, y]` in pixel coordinates.
#[tauri::command]
pub async fn init_snake(
    target_index: usize,
    init_radius_mm: Option<f64>,
    state: tauri::State<'_, Mutex<AppState>>,
) -> Result<Vec<[f64; 2]>, String> {
    let mut guard = state.lock().map_err(|e| format!("lock poisoned: {e}"))?;

    let targets = guard
        .annotation_targets
        .as_ref()
        .ok_or_else(|| "no annotation targets generated".to_string())?;

    let target = targets
        .get(target_index)
        .ok_or_else(|| format!("target_index {target_index} out of range (0..{})", targets.len()))?;

    let contour = if let Some(radius_mm) = init_radius_mm {
        // Create a circular contour at the requested radius.
        let mm_per_pixel = 2.0 * target.width_mm / target.pixels as f64;
        let center_px = target.pixels as f64 / 2.0;
        let radius_px = radius_mm / mm_per_pixel;
        let n_points = target.init_boundary.len(); // match point count
        init_circular_contour(center_px, center_px, radius_px, n_points)
    } else {
        // Use the target's pre-computed init_boundary.
        target.init_boundary.clone()
    };

    guard.snake_contours.insert(target_index, contour.clone());
    guard.finalized.insert(target_index, false);

    Ok(contour)
}

/// Evolve the active contour on a cross-section for `n_iterations` steps.
///
/// Returns the updated contour points and convergence information.
#[tauri::command]
pub async fn evolve_snake(
    target_index: usize,
    n_iterations: usize,
    state: tauri::State<'_, Mutex<AppState>>,
) -> Result<SnakeResult, String> {
    // Extract what we need under the lock, then release for the heavy computation.
    let (mut contour, image_vec, pixels) = {
        let guard = state.lock().map_err(|e| format!("lock poisoned: {e}"))?;

        let targets = guard
            .annotation_targets
            .as_ref()
            .ok_or_else(|| "no annotation targets generated".to_string())?;

        let target = targets
            .get(target_index)
            .ok_or_else(|| format!("target_index {target_index} out of range"))?;

        let contour = guard
            .snake_contours
            .get(&target_index)
            .ok_or_else(|| format!("no snake initialized for target {target_index} — call init_snake first"))?
            .clone();

        (contour, target.image.clone(), target.pixels)
    };

    // Run evolution on a blocking thread.
    let result = tokio::task::spawn_blocking(move || {
        // Convert flat Vec<f32> to Array2<f32>.
        let image = Array2::from_shape_vec((pixels, pixels), image_vec)
            .map_err(|e| format!("image reshape failed: {e}"))?;

        let (grad_x, grad_y) = compute_gradient_field(&image);

        let params = SnakeParams::default();
        let info = pipeline_evolve_snake(
            &mut contour,
            &grad_x,
            &grad_y,
            &image,
            &params,
            n_iterations,
        );

        Ok::<_, String>(SnakeResult {
            points: contour,
            iterations: info.iterations,
            max_displacement: info.max_displacement,
            converged: info.converged,
        })
    })
    .await
    .map_err(|e| format!("evolve_snake task failed: {e}"))??;

    // Store the updated contour back in state.
    {
        let mut guard = state.lock().map_err(|e| format!("lock poisoned: {e}"))?;
        guard
            .snake_contours
            .insert(target_index, result.points.clone());
    }

    Ok(result)
}

/// Replace the snake contour points for a cross-section (e.g. after user drag).
///
/// Also marks the contour as not finalized.
#[tauri::command]
pub async fn update_snake_points(
    target_index: usize,
    points: Vec<[f64; 2]>,
    state: tauri::State<'_, Mutex<AppState>>,
) -> Result<(), String> {
    let mut guard = state.lock().map_err(|e| format!("lock poisoned: {e}"))?;

    // Validate that annotation targets exist and index is in range.
    let targets = guard
        .annotation_targets
        .as_ref()
        .ok_or_else(|| "no annotation targets generated".to_string())?;

    if target_index >= targets.len() {
        return Err(format!(
            "target_index {target_index} out of range (0..{})",
            targets.len()
        ));
    }

    guard.snake_contours.insert(target_index, points);
    guard.finalized.insert(target_index, false);

    Ok(())
}

/// Add a control point to the snake contour at the given position.
///
/// The point is inserted at the closest edge of the existing contour.
/// Returns the index of the inserted point.
#[tauri::command]
pub async fn add_snake_point(
    target_index: usize,
    position: [f64; 2],
    state: tauri::State<'_, Mutex<AppState>>,
) -> Result<usize, String> {
    let mut guard = state.lock().map_err(|e| format!("lock poisoned: {e}"))?;

    let contour = guard
        .snake_contours
        .get_mut(&target_index)
        .ok_or_else(|| {
            format!("no snake initialized for target {target_index} — call init_snake first")
        })?;

    let idx = insert_control_point(contour, position);

    // Mark as not finalized since the contour changed.
    guard.finalized.insert(target_index, false);

    Ok(idx)
}

/// Number of control points kept when adopting the auto-detected vessel wall
/// as an editable contour. 16 points = one every 22.5°, the classical density
/// for smooth closed medical contours — loose enough to drag a single point
/// cleanly on a ~2mm-radius lumen (~13 px apart on the editor canvas) while
/// still capturing mild eccentricity. The underlying vessel_wall polygon
/// keeps its full 360 samples for display.
const VESSEL_WALL_CONTOUR_POINTS: usize = 16;

/// Subsample a dense closed polygon to `target_count` evenly-spaced points
/// (by input index). Input points are assumed roughly uniform — the vessel
/// wall scanner emits 360 samples evenly in angle around the lumen center.
fn resample_polygon(points: &[[f64; 2]], target_count: usize) -> Vec<[f64; 2]> {
    if target_count == 0 || points.len() <= target_count {
        return points.to_vec();
    }
    let n = points.len();
    let mut out = Vec::with_capacity(target_count);
    for i in 0..target_count {
        let src_idx = (i * n) / target_count;
        out.push(points[src_idx]);
    }
    out
}

/// Per-target resampled vessel-wall contour, returned to the frontend so the
/// canvas and backend stay in sync without an extra round-trip.
#[derive(Serialize)]
pub struct AdoptedContour {
    pub target_index: usize,
    pub points: Vec<[f64; 2]>,
}

/// Adopt the auto-detected vessel wall polygon as the finalized contour for a
/// cross-section, bypassing snake evolution entirely.
///
/// Resamples `target.vessel_wall` to `VESSEL_WALL_CONTOUR_POINTS` control
/// points (comfortable dragging density), stores the result in
/// `snake_contours[target_index]`, and marks the entry as finalized. If
/// `all == true`, does the same for every target that has a non-empty
/// vessel wall.
///
/// Returns the resampled contours so the frontend can display the same
/// points the backend now holds.
#[tauri::command]
pub async fn use_vessel_wall_as_contour(
    target_index: Option<usize>,
    all: Option<bool>,
    state: tauri::State<'_, Mutex<AppState>>,
) -> Result<Vec<AdoptedContour>, String> {
    let mut guard = state.lock().map_err(|e| format!("lock poisoned: {e}"))?;

    let targets = guard
        .annotation_targets
        .as_ref()
        .ok_or_else(|| "no annotation targets generated".to_string())?
        .clone();

    let do_all = all.unwrap_or(false);
    let mut out = Vec::new();

    if do_all {
        for (idx, target) in targets.iter().enumerate() {
            if target.vessel_wall.is_empty() {
                continue;
            }
            let contour = resample_polygon(&target.vessel_wall, VESSEL_WALL_CONTOUR_POINTS);
            guard.snake_contours.insert(idx, contour.clone());
            guard.finalized.insert(idx, true);
            out.push(AdoptedContour { target_index: idx, points: contour });
        }
    } else {
        let idx = target_index.ok_or_else(|| {
            "target_index must be provided when all=false".to_string()
        })?;
        let target = targets.get(idx).ok_or_else(|| {
            format!("target_index {idx} out of range (0..{})", targets.len())
        })?;
        if target.vessel_wall.is_empty() {
            return Err(format!(
                "target {idx} has no auto-detected vessel wall"
            ));
        }
        let contour = resample_polygon(&target.vessel_wall, VESSEL_WALL_CONTOUR_POINTS);
        guard.snake_contours.insert(idx, contour.clone());
        guard.finalized.insert(idx, true);
        out.push(AdoptedContour { target_index: idx, points: contour });
    }

    Ok(out)
}

/// Mark a cross-section's contour as finalized.
#[tauri::command]
pub async fn finalize_contour(
    target_index: usize,
    state: tauri::State<'_, Mutex<AppState>>,
) -> Result<(), String> {
    let mut guard = state.lock().map_err(|e| format!("lock poisoned: {e}"))?;

    // Validate the target exists and has a contour.
    let targets = guard
        .annotation_targets
        .as_ref()
        .ok_or_else(|| "no annotation targets generated".to_string())?;

    if target_index >= targets.len() {
        return Err(format!(
            "target_index {target_index} out of range (0..{})",
            targets.len()
        ));
    }

    if !guard.snake_contours.contains_key(&target_index) {
        return Err(format!(
            "no snake contour for target {target_index} — call init_snake first"
        ));
    }

    guard.finalized.insert(target_index, true);

    Ok(())
}

/// Build a 3D ROI mask from every cross-section that currently has a contour
/// and run multi-material decomposition on the masked region.
///
/// There is no "accept" step: any contour present in `snake_contours` (from
/// auto-adoption or manual drag) is included. Returns summary statistics
/// (mean fractions over the mask).
#[tauri::command]
pub async fn run_mmd_on_roi(
    method: String,
    app: tauri::AppHandle,
    state: tauri::State<'_, Mutex<AppState>>,
) -> Result<MmdSummary, String> {
    // Validate method. Noise-aware water/lipid GLS is the only solver; the
    // 3-material direct/PWSQS solvers were dropped (ill-conditioned at
    // soft-tissue contrast).
    if method != "gls" {
        return Err(format!(
            "unknown method '{method}': only 'gls' (noise-aware water/lipid) is supported"
        ));
    }

    // Extract everything we need from state under the lock.
    let (
        targets,
        finalized_contours,
        finalized_indices,
        frame,
        dual_energy,
        volume_shape,
        spacing,
        origin,
        direction,
        existing_calib,
    ) = {
        let guard = state.lock().map_err(|e| format!("lock poisoned: {e}"))?;

        let targets = guard
            .annotation_targets
            .as_ref()
            .ok_or_else(|| "no annotation targets generated".to_string())?;

        // Collect every cross-section that has a contour, sorted by target
        // index (required by build_3d_roi_mask). No finalized filter.
        let mut contour_pairs: Vec<(usize, Vec<[f64; 2]>)> = guard
            .snake_contours
            .iter()
            .map(|(&idx, c)| (idx, c.clone()))
            .collect();

        if contour_pairs.is_empty() {
            return Err("no contours — open the MMD tab to auto-adopt the lumen wall".into());
        }

        contour_pairs.sort_by_key(|(idx, _)| *idx);

        let finalized_contours: Vec<Vec<[f64; 2]>> =
            contour_pairs.iter().map(|(_, c)| c.clone()).collect();

        let finalized_indices: Vec<usize> = contour_pairs
            .iter()
            .map(|(idx, _)| targets[*idx].frame_index)
            .collect();

        let frame = guard
            .cpr_frame
            .as_ref()
            .ok_or_else(|| "no CPR frame built".to_string())?;

        let de = guard
            .dual_energy
            .as_ref()
            .ok_or_else(|| "no dual-energy volume loaded".to_string())?;

        // Build the mask in the dual-energy grid — MMD runs against `de.low`
        // / `de.high`, so the mask shape must match those volumes exactly.
        // If the currently-active `state.volume` is a different series
        // (e.g. CCTA) with different dims, using its shape here would panic
        // in pwsqs with a mask-shape assertion. The CPR frame + contours
        // live in patient mm so they're reusable across grids.
        let de_low_shape = de.low.shape();
        let volume_shape = [de_low_shape[0], de_low_shape[1], de_low_shape[2]];

        // Get cross-section params from the first target.
        let cs_width_mm = targets[0].width_mm;
        let cs_pixels = targets[0].pixels;

        // GLS reuses a prior whole-volume calibration if present; otherwise the
        // blocking task self-calibrates and we store the result below.
        let existing_calib = guard.wl_calibration.clone();

        (
            (cs_width_mm, cs_pixels),
            finalized_contours,
            finalized_indices,
            clone_frame(frame),
            (Arc::clone(&de.low), Arc::clone(&de.high), de.low_energy_kev, de.high_energy_kev),
            volume_shape,
            de.spacing,
            de.origin,
            de.direction,
            existing_calib,
        )
    };

    let (cs_width_mm, cs_pixels) = targets;
    let (low_energy, high_energy, low_kev, high_kev) = dual_energy;
    let method_clone = method.clone();
    let app_clone = app.clone();

    // Run on a blocking thread.
    let (mmd_result, summary, fresh_calib) = tokio::task::spawn_blocking(move || {
        let _ = app_clone.emit(
            "annotation-progress",
            AnnotationProgress {
                stage: "building_roi_mask".into(),
                progress: 0.1,
            },
        );

        // Build a ring ROI in the pericoronary region: starts at the user's
        // finalized lumen contour and extends outward by PVAT_RING_OUTER_MM.
        // Matches the radial-angular surface-plot extent so every analysis
        // view samples the same tissue.
        const PVAT_RING_OUTER_MM: f64 = 10.0;
        let mask = roi::build_3d_roi_mask(
            &finalized_contours,
            &frame,
            &finalized_indices,
            volume_shape,
            spacing,
            origin,
            &direction,
            cs_width_mm,
            cs_pixels,
            PVAT_RING_OUTER_MM,
        );

        let n_voxels = mask.iter().filter(|&&v| v).count();
        if n_voxels == 0 {
            return Err("ROI mask is empty — no voxels inside the finalized contours".to_string());
        }

        let _ = app_clone.emit(
            "annotation-progress",
            AnnotationProgress {
                stage: "decomposing".into(),
                progress: 0.3,
            },
        );

        // Noise-aware GLS water/lipid decomposition on the ROI mask, theoretical
        // anchor (f_l=1 ⇒ pure lipid). Reuse the whole-volume calibration if the
        // user already ran it; otherwise self-calibrate from this volume and
        // persist it (`fresh_calib` is Some only when freshly measured).
        let _ = app_clone.emit(
            "annotation-progress",
            AnnotationProgress { stage: "wl_calibrating".into(), progress: 0.35 },
        );
        let (calib, fresh_calib): (WlCalibration, Option<WlCalibration>) = match existing_calib {
            Some(c) => (c, None),
            None => {
                let c = mmd::self_calibrate(&low_energy, &high_energy, low_kev, high_kev)?;
                (c.clone(), Some(c))
            }
        };
        let result: MmdResult =
            mmd::decompose_volume_gls(&low_energy, &high_energy, &mask, &calib, WlAnchor::Theoretical);

        let _ = app_clone.emit(
            "annotation-progress",
            AnnotationProgress {
                stage: "computing_stats".into(),
                progress: 0.95,
            },
        );

        // Compute mean fractions over the mask.
        let mask_slice = result.mask.as_slice().unwrap();
        let wf_slice = result.water_frac.as_slice().unwrap();
        let lf_slice = result.lipid_frac.as_slice().unwrap();
        let if_slice = result.iodine_frac.as_slice().unwrap();
        let cf_slice = result.calcium_frac.as_slice().unwrap();

        let mut sum_w = 0.0_f64;
        let mut sum_l = 0.0_f64;
        let mut sum_i = 0.0_f64;
        let mut sum_c = 0.0_f64;

        for idx in 0..mask_slice.len() {
            if mask_slice[idx] {
                sum_w += wf_slice[idx] as f64;
                sum_l += lf_slice[idx] as f64;
                sum_i += if_slice[idx] as f64;
                sum_c += cf_slice[idx] as f64;
            }
        }

        let n = n_voxels as f64;
        let summary = MmdSummary {
            method: method_clone,
            iterations: result.iterations,
            converged: result.converged,
            n_voxels,
            mean_water_frac: sum_w / n,
            mean_lipid_frac: sum_l / n,
            mean_iodine_frac: sum_i / n,
            mean_calcium_frac: sum_c / n,
        };

        Ok((result, summary, fresh_calib))
    })
    .await
    .map_err(|e| format!("run_mmd_on_roi task failed: {e}"))??;

    // Store the MmdResult in state; persist a freshly-measured GLS calibration
    // so a later whole-volume Water/Lipid view reuses it.
    {
        let mut guard = state.lock().map_err(|e| format!("lock poisoned: {e}"))?;
        guard.mmd_result = Some(mmd_result);
        if let Some(c) = fresh_calib {
            guard.wl_calibration = Some(c);
        }
    }

    let _ = app.emit(
        "annotation-progress",
        AnnotationProgress {
            stage: "done".into(),
            progress: 1.0,
        },
    );

    Ok(summary)
}

// ---------------------------------------------------------------------------
// Helper: select material array from MmdResult
// ---------------------------------------------------------------------------

fn select_material_array<'a>(
    mmd_result: &'a MmdResult,
    material: &str,
    unit: &str,
) -> Result<&'a ndarray::Array3<f32>, String> {
    match (material, unit) {
        ("water", "fraction") => Ok(&mmd_result.water_frac),
        ("water", "mass") => Ok(&mmd_result.water_mass),
        ("lipid", "fraction") => Ok(&mmd_result.lipid_frac),
        ("lipid", "mass") => Ok(&mmd_result.lipid_mass),
        ("iodine", "fraction") => Ok(&mmd_result.iodine_frac),
        ("iodine", "mass") => Ok(&mmd_result.iodine_mass),
        ("calcium", "fraction") => Ok(&mmd_result.calcium_frac),
        ("calcium", "mass") => Ok(&mmd_result.calcium_mass),
        ("density", _) => Ok(&mmd_result.total_density),
        _ => Err(format!("unknown material/unit: {material}/{unit}")),
    }
}

// ---------------------------------------------------------------------------
// Surface sampling command
// ---------------------------------------------------------------------------

/// Sample radial-angular surface data from the MMD result.
/// Returns surface data for all finalized cross-sections.
#[tauri::command]
pub async fn sample_surfaces(
    material: String,
    unit: String,
    state: tauri::State<'_, Mutex<AppState>>,
) -> Result<Vec<CrossSectionSurface>, String> {
    // Extract everything under the lock, then release.
    let (material_map, frame, targets, finalized_contours, spacing, origin, direction) = {
        let guard = state.lock().map_err(|e| format!("lock poisoned: {e}"))?;

        let mmd_result = guard
            .mmd_result
            .as_ref()
            .ok_or_else(|| "no MMD result — run decomposition first".to_string())?;

        let material_map = select_material_array(mmd_result, &material, &unit)?.clone();

        let targets = guard
            .annotation_targets
            .as_ref()
            .ok_or_else(|| "no annotation targets generated".to_string())?
            .clone();

        // Use every cross-section that has a contour (no finalized filter).
        let finalized_contours: HashMap<usize, Vec<[f64; 2]>> = guard
            .snake_contours
            .iter()
            .map(|(&idx, c)| (idx, c.clone()))
            .collect();

        if finalized_contours.is_empty() {
            return Err("no contours — open the MMD tab to auto-adopt the lumen wall".into());
        }

        let frame = guard
            .cpr_frame
            .as_ref()
            .ok_or_else(|| "no CPR frame built".to_string())?;
        let frame = clone_frame(frame);

        // Sample in the dual-energy grid (the material maps live there),
        // not the currently-active `state.volume` which may be a different
        // series (CCTA, CaScore) after patient-wide loading.
        let de = guard
            .dual_energy
            .as_ref()
            .ok_or_else(|| "no dual-energy volume loaded".to_string())?;

        (
            material_map,
            frame,
            targets,
            finalized_contours,
            de.spacing,
            de.origin,
            de.direction,
        )
    };

    // Run sampling on a blocking thread.
    let surfaces = tokio::task::spawn_blocking(move || {
        let params = RadialAngularParams::default();
        radial_angular::sample_radial_angular(
            &material_map,
            &frame,
            &targets,
            &finalized_contours,
            spacing,
            origin,
            &direction,
            &params,
        )
    })
    .await
    .map_err(|e| format!("sample_surfaces task failed: {e}"))?;

    Ok(surfaces)
}

// ---------------------------------------------------------------------------
// MMD overlay for a single cross-section
// ---------------------------------------------------------------------------

/// Get the current MMD result for a specific material as a flat array for overlay rendering.
///
/// Extracts the material values for the cross-section slice at the given target,
/// for rendering as a colormap overlay in the editor.
#[tauri::command]
pub async fn get_mmd_overlay(
    target_index: usize,
    material: String,
    unit: String,
    state: tauri::State<'_, Mutex<AppState>>,
) -> Result<Vec<f32>, String> {
    let guard = state.lock().map_err(|e| format!("lock poisoned: {e}"))?;

    let mmd_result = guard
        .mmd_result
        .as_ref()
        .ok_or_else(|| "no MMD result — run decomposition first".to_string())?;

    let material_map = select_material_array(mmd_result, &material, &unit)?;

    let targets = guard
        .annotation_targets
        .as_ref()
        .ok_or_else(|| "no annotation targets generated".to_string())?;

    let target = targets
        .get(target_index)
        .ok_or_else(|| format!("target_index {target_index} out of range (0..{})", targets.len()))?;

    let frame = guard
        .cpr_frame
        .as_ref()
        .ok_or_else(|| "no CPR frame built".to_string())?;

    // Overlay samples the MMD material map, which lives in the dual-energy
    // grid — pull geometry from state.dual_energy, not state.volume (the
    // currently-displayed series may be a different volume entirely).
    let de = guard
        .dual_energy
        .as_ref()
        .ok_or_else(|| "no dual-energy volume loaded".to_string())?;

    let spacing = de.spacing;
    let origin = de.origin;
    let direction = de.direction;

    let frame_idx = target.frame_index;
    if frame_idx >= frame.n_cols() {
        return Err(format!("frame_index {frame_idx} out of range"));
    }

    let pos_mm = frame.positions[frame_idx];
    let normal = frame.normals[frame_idx];
    let binormal = frame.binormals[frame_idx];

    let pixels = target.pixels;
    let width_mm = target.width_mm;
    let inv_spacing = [1.0 / spacing[0], 1.0 / spacing[1], 1.0 / spacing[2]];

    let mut overlay = vec![f32::NAN; pixels * pixels];

    for row in 0..pixels {
        for col in 0..pixels {
            // Convert pixel (row, col) to offset in mm from center.
            // Row direction: row=0 -> +normal, row=pixels-1 -> -normal
            let offset_n = width_mm * (1.0 - 2.0 * row as f64 / (pixels as f64 - 1.0));
            // Col direction: col=0 -> +binormal, col=pixels-1 -> -binormal
            let offset_b = width_mm * (1.0 - 2.0 * col as f64 / (pixels as f64 - 1.0));

            let wz = pos_mm[0] + offset_n * normal[0] + offset_b * binormal[0];
            let wy = pos_mm[1] + offset_n * normal[1] + offset_b * binormal[1];
            let wx = pos_mm[2] + offset_n * normal[2] + offset_b * binormal[2];

            // Convert world coords to voxel indices (IOP-aware).
            let [vz, vy, vx] = pcat_pipeline::types::patient_to_voxel(
                [wz, wy, wx],
                origin,
                inv_spacing,
                &direction,
            );

            // Only voxels inside the ROI mask have meaningful material
            // fractions (solver only runs there). Elsewhere the backing
            // array is zero — return NaN so the UI can fall back to CT.
            let mask_shape = mmd_result.mask.shape();
            let zi = vz.round() as i64;
            let yi = vy.round() as i64;
            let xi = vx.round() as i64;
            let inside_volume = zi >= 0
                && yi >= 0
                && xi >= 0
                && (zi as usize) < mask_shape[0]
                && (yi as usize) < mask_shape[1]
                && (xi as usize) < mask_shape[2];
            let inside_mask = inside_volume
                && mmd_result.mask[[zi as usize, yi as usize, xi as usize]];
            if !inside_mask {
                overlay[row * pixels + col] = f32::NAN;
                continue;
            }

            let val = pcat_pipeline::interp::trilinear(material_map, vz, vy, vx);
            overlay[row * pixels + col] = val;
        }
    }

    Ok(overlay)
}

// ---------------------------------------------------------------------------
// Save / Load annotation state
// ---------------------------------------------------------------------------

/// Saved annotation state for one patient.
#[derive(Serialize, Deserialize)]
pub struct AnnotationStateJson {
    pub centerline_mm: Vec<[f64; 3]>,
    pub snake_contours: HashMap<usize, Vec<[f64; 2]>>,
    pub finalized: HashMap<usize, bool>,
    pub mmd_method: Option<String>,
    pub mmd_iterations: Option<usize>,
    pub mmd_converged: Option<bool>,
}

/// Sanitize a path string into a safe filename component.
fn sanitize_for_filename(s: &str) -> String {
    s.chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() || c == '-' || c == '_' || c == '.' {
                c
            } else {
                '_'
            }
        })
        .collect::<String>()
        .replace("..", "_")
}

/// Save the current annotation state for the given patient (DICOM folder path).
/// Stores under `app_data_dir/annotations/<sanitized-name>.json`.
#[tauri::command]
pub async fn save_annotations(
    app: tauri::AppHandle,
    dicom_path: String,
    centerline_mm: Vec<[f64; 3]>,
    state: tauri::State<'_, Mutex<AppState>>,
) -> Result<String, String> {
    // Read snake_contours and finalized from state.
    let (snake_contours, finalized, mmd_meta) = {
        let guard = state.lock().map_err(|e| format!("lock poisoned: {e}"))?;
        let sc = guard.snake_contours.clone();
        let fin = guard.finalized.clone();
        let mmd_meta = guard.mmd_result.as_ref().map(|r| {
            (
                "pwsqs".to_string(), // method is not stored in MmdResult, default
                r.iterations,
                r.converged,
            )
        });
        (sc, fin, mmd_meta)
    };

    let annotation_state = AnnotationStateJson {
        centerline_mm,
        snake_contours,
        finalized,
        mmd_method: mmd_meta.as_ref().map(|(m, _, _)| m.clone()),
        mmd_iterations: mmd_meta.as_ref().map(|(_, i, _)| *i),
        mmd_converged: mmd_meta.as_ref().map(|(_, _, c)| *c),
    };

    let dir = app
        .path()
        .app_data_dir()
        .expect("app data dir")
        .join("annotations");
    let _ = std::fs::create_dir_all(&dir);

    let key = Path::new(&dicom_path)
        .file_name()
        .map(|f| f.to_string_lossy().to_string())
        .unwrap_or_else(|| sanitize_for_filename(&dicom_path));
    let path = dir.join(format!("{}.json", sanitize_for_filename(&key)));

    let json = serde_json::to_string_pretty(&annotation_state)
        .map_err(|e| format!("serialize failed: {e}"))?;
    std::fs::write(&path, &json).map_err(|e| format!("write failed: {e}"))?;

    Ok(path.to_string_lossy().to_string())
}

/// Load saved annotation state for the given patient. Restores into AppState.
/// Returns the loaded state for the frontend to display, or None if no save exists.
#[tauri::command]
pub async fn load_annotations(
    app: tauri::AppHandle,
    dicom_path: String,
    state: tauri::State<'_, Mutex<AppState>>,
) -> Result<Option<AnnotationStateJson>, String> {
    let dir = app
        .path()
        .app_data_dir()
        .expect("app data dir")
        .join("annotations");
    let key = Path::new(&dicom_path)
        .file_name()
        .map(|f| f.to_string_lossy().to_string())
        .unwrap_or_else(|| sanitize_for_filename(&dicom_path));
    let path = dir.join(format!("{}.json", sanitize_for_filename(&key)));

    if !path.exists() {
        return Ok(None);
    }

    let data = std::fs::read_to_string(&path).map_err(|e| format!("read failed: {e}"))?;
    let annotation_state: AnnotationStateJson =
        serde_json::from_str(&data).map_err(|e| format!("parse failed: {e}"))?;

    // Restore snake_contours and finalized into AppState.
    {
        let mut guard = state.lock().map_err(|e| format!("lock poisoned: {e}"))?;
        guard.snake_contours = annotation_state.snake_contours.clone();
        guard.finalized = annotation_state.finalized.clone();
    }

    Ok(Some(annotation_state))
}

// ---------------------------------------------------------------------------
// CSV export of MMD surface data
// ---------------------------------------------------------------------------

/// Material key → (material name, unit name) for select_material_array.
const MATERIAL_KEYS: &[(&str, &str, &str)] = &[
    ("lipid_frac", "lipid", "fraction"),
    ("lipid_mass", "lipid", "mass"),
    ("water_frac", "water", "fraction"),
    ("water_mass", "water", "mass"),
    ("iodine_frac", "iodine", "fraction"),
    ("iodine_mass", "iodine", "mass"),
    ("calcium_frac", "calcium", "fraction"),
    ("calcium_mass", "calcium", "mass"),
    ("total_density", "density", "fraction"), // unit doesn't matter for density
];

/// Export current MMD surface data as CSV string.
/// Includes one row per (target_index, theta, r) sample point with all materials.
#[tauri::command]
pub async fn export_mmd_csv(
    patient_id: String,
    state: tauri::State<'_, Mutex<AppState>>,
) -> Result<String, String> {
    // Extract everything we need from state under one lock.
    let (mmd_result_arrays, frame, targets, finalized_contours, spacing, origin, direction) = {
        let guard = state.lock().map_err(|e| format!("lock poisoned: {e}"))?;

        let mmd_result = guard
            .mmd_result
            .as_ref()
            .ok_or_else(|| "no MMD result — run decomposition first".to_string())?;

        // Clone all 9 material arrays.
        let arrays: Vec<(&str, ndarray::Array3<f32>)> = MATERIAL_KEYS
            .iter()
            .map(|(key, mat, unit)| {
                let arr = select_material_array(mmd_result, mat, unit)
                    .expect("known material key")
                    .clone();
                (*key, arr)
            })
            .collect();

        let targets = guard
            .annotation_targets
            .as_ref()
            .ok_or_else(|| "no annotation targets generated".to_string())?
            .clone();

        let finalized_contours: HashMap<usize, Vec<[f64; 2]>> = guard
            .snake_contours
            .iter()
            .map(|(&idx, c)| (idx, c.clone()))
            .collect();

        if finalized_contours.is_empty() {
            return Err("no contours — open the MMD tab to auto-adopt the lumen wall".into());
        }

        let frame = guard
            .cpr_frame
            .as_ref()
            .ok_or_else(|| "no CPR frame built".to_string())?;
        let frame = clone_frame(frame);

        // CSV surfaces use the same radial-angular sampler as sample_surfaces;
        // geometry must come from the dual-energy grid.
        let de = guard
            .dual_energy
            .as_ref()
            .ok_or_else(|| "no dual-energy volume loaded".to_string())?;

        (
            arrays,
            frame,
            targets,
            finalized_contours,
            de.spacing,
            de.origin,
            de.direction,
        )
    };

    // Sample all 9 material surfaces on a blocking thread.
    let csv_string = tokio::task::spawn_blocking(move || {
        let params = RadialAngularParams::default();

        // Sample each material → Vec<CrossSectionSurface>.
        let mut all_surfaces: Vec<(&str, Vec<CrossSectionSurface>)> = Vec::new();
        for (key, array) in &mmd_result_arrays {
            let surfaces = radial_angular::sample_radial_angular(
                array,
                &frame,
                &targets,
                &finalized_contours,
                spacing,
                origin,
                &direction,
                &params,
            );
            all_surfaces.push((key, surfaces));
        }

        // All material surfaces have the same grid layout (same targets × same theta × same r),
        // so we can iterate them together.

        if all_surfaces.is_empty() || all_surfaces[0].1.is_empty() {
            return Err("no surface data to export".to_string());
        }

        let n_surfaces = all_surfaces[0].1.len(); // number of cross-sections

        let mut csv = String::with_capacity(1024 * 1024);

        // Header
        csv.push_str("patient_id,target_index,arc_mm,theta_deg,r_mm,lipid_frac,lipid_mass,water_frac,water_mass,iodine_frac,iodine_mass,calcium_frac,calcium_mass,total_density\n");

        // Build a lookup from key to index in all_surfaces for ordered access.
        let key_order = [
            "lipid_frac",
            "lipid_mass",
            "water_frac",
            "water_mass",
            "iodine_frac",
            "iodine_mass",
            "calcium_frac",
            "calcium_mass",
            "total_density",
        ];
        let key_to_idx: HashMap<&str, usize> = all_surfaces
            .iter()
            .enumerate()
            .map(|(i, (k, _))| (*k, i))
            .collect();

        for cs_idx in 0..n_surfaces {
            let ref_surface = &all_surfaces[0].1[cs_idx];
            let arc_mm = ref_surface.arc_mm;
            let n_theta = ref_surface.n_theta;
            let n_radial = ref_surface.n_radial;

            // Find the target index from arc_mm matching.
            // The surfaces are produced for finalized targets in order, so we need to
            // look up which target this corresponds to.
            let target_idx = targets
                .iter()
                .position(|t| (t.arc_mm - arc_mm).abs() < 1e-6)
                .unwrap_or(cs_idx);

            for i_theta in 0..n_theta {
                let theta = ref_surface.theta_deg[i_theta];
                for i_r in 0..n_radial {
                    let r = ref_surface.r_mm[i_r];
                    let flat_idx = i_theta * n_radial + i_r;

                    // Check if any material has a valid (non-NaN) value at this point.
                    let ref_val = ref_surface.surface[flat_idx];
                    if ref_val.is_nan() {
                        continue; // Skip NaN entries (beyond contour boundary).
                    }

                    // Collect values in key_order.
                    let mut vals = [f64::NAN; 9];
                    let mut any_nan = false;
                    for (order_idx, &key) in key_order.iter().enumerate() {
                        if let Some(&surf_idx) = key_to_idx.get(key) {
                            let v = all_surfaces[surf_idx].1[cs_idx].surface[flat_idx];
                            if v.is_nan() {
                                any_nan = true;
                                break;
                            }
                            vals[order_idx] = v as f64;
                        }
                    }

                    if any_nan {
                        continue;
                    }

                    use std::fmt::Write;
                    let _ = writeln!(
                        csv,
                        "{},{},{:.4},{:.2},{:.4},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6}",
                        patient_id,
                        target_idx,
                        arc_mm,
                        theta,
                        r,
                        vals[0], vals[1], vals[2], vals[3],
                        vals[4], vals[5], vals[6], vals[7], vals[8],
                    );
                }
            }
        }

        Ok(csv)
    })
    .await
    .map_err(|e| format!("export_mmd_csv task failed: {e}"))??;

    Ok(csv_string)
}
