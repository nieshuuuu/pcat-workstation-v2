use ndarray::Array3;
use serde::Serialize;

use crate::active_contour::init_circular_contour;
use crate::cpr::CprFrame;

/// A single cross-section prepared for annotation.
#[derive(Debug, Clone, Serialize)]
pub struct AnnotationTarget {
    /// Cross-section image (HU), row-major, pixels x pixels.
    pub image: Vec<f32>,
    /// Image dimension in pixels.
    pub pixels: usize,
    /// Physical width of the cross-section in mm (from center to edge on each side).
    pub width_mm: f64,
    /// Arc-length position along the centerline in mm.
    pub arc_mm: f64,
    /// Frame column index this cross-section corresponds to.
    pub frame_index: usize,
    /// Auto-detected vessel wall contour as [x,y] pixel coordinates.
    pub vessel_wall: Vec<[f64; 2]>,
    /// Vessel equivalent radius in mm (sqrt(area/pi)).
    pub vessel_radius_mm: f64,
    /// Initial snake boundary: circle at ~2x vessel equivalent radius from center, in [x,y] pixel coords.
    pub init_boundary: Vec<[f64; 2]>,
}

/// Parameters for generating annotation targets.
#[derive(Debug, Clone)]
pub struct AnnotationBatchParams {
    /// Number of cross-sections (default: 20).
    pub n_sections: usize,
    /// Spacing between cross-sections in mm (default: 2.0).
    pub section_spacing_mm: f64,
    /// Cross-section image half-width in mm (default: 15.0).
    pub width_mm: f64,
    /// Cross-section image pixel size (default: 128).
    pub pixels: usize,
    /// Number of points on snake contour (default: 72).
    pub n_snake_points: usize,
    /// Optional ostium position in patient (mm) coords, order [z, y, x].
    /// When set, sections start at the arc-length nearest this point.
    pub ostium_zyx: Option<[f64; 3]>,
}

impl Default for AnnotationBatchParams {
    fn default() -> Self {
        Self {
            n_sections: 20,
            section_spacing_mm: 2.0,
            width_mm: 15.0,
            pixels: 128,
            n_snake_points: 72,
            ostium_zyx: None,
        }
    }
}

/// Resolve the starting arc-length (mm) for a batch of annotation cross-sections.
///
/// If `ostium_zyx` is `Some`, finds the frame column nearest that patient-space
/// point and returns its arc-length, clamped so that at least one section fits
/// within the remaining centerline. If `None`, returns 0.0 (start at first
/// centerline waypoint — legacy behaviour).
///
/// Coordinate note: `CprFrame::positions` is stored in `[z, y, x]` patient mm
/// (see `cpr.rs::CprFrame`). `ostium_zyx` is in the same order, so we compare
/// component-wise without any swap.
fn resolve_start_arc(
    frame: &CprFrame,
    ostium_zyx: Option<[f64; 3]>,
    section_spacing_mm: f64,
) -> f64 {
    let Some(ostium) = ostium_zyx else {
        return 0.0;
    };

    let n_cols = frame.n_cols();
    if n_cols == 0 {
        return 0.0;
    }

    // Argmin squared distance over frame columns.
    let mut best_idx = 0usize;
    let mut best_d2 = f64::INFINITY;
    for (i, p) in frame.positions.iter().enumerate() {
        // Both p and ostium are [z, y, x].
        let dz = p[0] - ostium[0];
        let dy = p[1] - ostium[1];
        let dx = p[2] - ostium[2];
        let d2 = dx * dx + dy * dy + dz * dz;
        if d2 < best_d2 {
            best_d2 = d2;
            best_idx = i;
        }
    }

    // Clamp so we always produce at least one section when possible.
    // If the ostium is downstream of the last frame column, start_arc falls
    // back to (total_arc − section_spacing_mm) so we still emit one target
    // rather than silently producing zero.
    let total_arc = frame.arclengths[n_cols - 1];
    let raw = frame.arclengths[best_idx];
    raw.min((total_arc - section_spacing_mm).max(0.0))
}

/// Generate annotation targets for a segment of the centerline.
///
/// Takes a precomputed CprFrame and produces `n_sections` cross-section images
/// at equal spacing, with auto-detected vessel walls and initial snake boundaries.
pub fn generate_annotation_batch(
    frame: &CprFrame,
    volume: &Array3<f32>,
    spacing: [f64; 3],
    origin: [f64; 3],
    direction: &[f64; 9],
    params: &AnnotationBatchParams,
) -> Vec<AnnotationTarget> {
    let n_cols = frame.n_cols();
    if n_cols < 2 {
        return Vec::new();
    }

    let total_arc = frame.arclengths[n_cols - 1];
    let n_angles = 360usize;
    let mm_per_pixel = 2.0 * params.width_mm / params.pixels as f64;
    let center_px = params.pixels as f64 / 2.0;

    // If an ostium is provided, shift the start of the batch to its arc-length;
    // otherwise start at 0 (first centerline waypoint) as before.
    let start_arc_mm = resolve_start_arc(frame, params.ostium_zyx, params.section_spacing_mm);

    let mut targets = Vec::with_capacity(params.n_sections);

    for section in 0..params.n_sections {
        let arc_mm = start_arc_mm + section as f64 * params.section_spacing_mm;

        // Stop if we exceed the centerline length.
        if arc_mm > total_arc {
            break;
        }

        // Find nearest frame column index for this arc-length.
        let frame_index = find_nearest_arc_index(&frame.arclengths, arc_mm);

        // Render cross-section at this position.
        let position_frac = frame_index as f64 / (n_cols - 1) as f64;
        let cs = frame.render_cross_section(
            volume,
            spacing,
            origin,
            direction,
            position_frac,
            0.0, // rotation_deg = 0
            params.width_mm,
            params.pixels,
        );

        // Auto-detect vessel wall + equivalent radius from the 2-D cross-section.
        let geom = crate::vessel_wall::compute_vessel_geometry(
            &cs.image,
            params.pixels,
            params.width_mm,
            n_angles,
        );
        let vessel_wall = geom.wall;
        let vessel_radius_mm = geom.diameter_mm / 2.0;

        // Initial snake boundary: circle at 2 * r_eq from center (1 diameter out from wall).
        let init_radius_mm = (2.0 * vessel_radius_mm).min(params.width_mm * 0.9);
        let init_radius_px = init_radius_mm / mm_per_pixel;
        let init_boundary =
            init_circular_contour(center_px, center_px, init_radius_px, params.n_snake_points);

        targets.push(AnnotationTarget {
            image: cs.image,
            pixels: params.pixels,
            width_mm: params.width_mm,
            arc_mm: cs.arc_mm,
            frame_index,
            vessel_wall,
            vessel_radius_mm,
            init_boundary,
        });
    }

    targets
}

/// Find the index in `arclengths` closest to `target_arc`.
fn find_nearest_arc_index(arclengths: &[f64], target_arc: f64) -> usize {
    let mut best_idx = 0;
    let mut best_diff = f64::INFINITY;
    for (i, &s) in arclengths.iter().enumerate() {
        let diff = (s - target_arc).abs();
        if diff < best_diff {
            best_diff = diff;
            best_idx = i;
        }
    }
    best_idx
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cpr::CprFrame;
    use ndarray::Array3;

    /// Create a synthetic volume with a contrast-enhanced tube (cylinder) along z-axis
    /// centered at (y=32, x=32). The tube has HU=300 inside, -100 outside.
    fn make_vessel_phantom() -> (Array3<f32>, [f64; 3], [f64; 3]) {
        let size = 64;
        let mut vol = Array3::<f32>::from_elem((size, size, size), -100.0);
        let center_y = 32.0_f64;
        let center_x = 32.0_f64;
        let vessel_radius = 3.0_f64; // 3mm radius vessel

        for z in 0..size {
            for y in 0..size {
                for x in 0..size {
                    let dy = y as f64 - center_y;
                    let dx = x as f64 - center_x;
                    let dist = (dy * dy + dx * dx).sqrt();
                    if dist <= vessel_radius {
                        vol[[z, y, x]] = 300.0;
                    }
                }
            }
        }

        let spacing = [1.0, 1.0, 1.0];
        let origin = [0.0, 0.0, 0.0];
        (vol, spacing, origin)
    }

    /// Build a CprFrame along the center of the phantom vessel.
    fn make_frame(n_cols: usize) -> CprFrame {
        let centerline_mm: Vec<[f64; 3]> = (2..62)
            .map(|z| [z as f64, 32.0, 32.0])
            .collect();
        CprFrame::from_centerline(&centerline_mm, n_cols)
    }

    #[test]
    fn test_correct_count() {
        let (vol, spacing, origin) = make_vessel_phantom();
        let frame = make_frame(200);
        let params = AnnotationBatchParams::default();
        let targets = generate_annotation_batch(&frame, &vol, spacing, origin, &crate::types::IDENTITY_DIRECTION, &params);
        assert_eq!(
            targets.len(),
            params.n_sections,
            "Should produce exactly {} targets, got {}",
            params.n_sections,
            targets.len()
        );
    }

    #[test]
    fn test_arc_length_spacing() {
        let (vol, spacing, origin) = make_vessel_phantom();
        let frame = make_frame(200);
        let params = AnnotationBatchParams::default();
        let targets = generate_annotation_batch(&frame, &vol, spacing, origin, &crate::types::IDENTITY_DIRECTION, &params);

        // Verify consecutive targets are spaced at approximately section_spacing_mm.
        for i in 1..targets.len() {
            let delta = targets[i].arc_mm - targets[i - 1].arc_mm;
            assert!(
                (delta - params.section_spacing_mm).abs() < 1.0,
                "Spacing between target {} and {}: {:.2} mm, expected ~{:.1} mm",
                i - 1,
                i,
                delta,
                params.section_spacing_mm
            );
        }
    }

    #[test]
    fn test_vessel_wall_exists() {
        let (vol, spacing, origin) = make_vessel_phantom();
        let frame = make_frame(200);
        let params = AnnotationBatchParams::default();
        let targets = generate_annotation_batch(&frame, &vol, spacing, origin, &crate::types::IDENTITY_DIRECTION, &params);

        for (i, t) in targets.iter().enumerate() {
            assert!(
                !t.vessel_wall.is_empty(),
                "Target {} should have a non-empty vessel wall contour",
                i
            );
            assert!(
                t.vessel_radius_mm > 0.0,
                "Target {} should have positive vessel radius, got {}",
                i,
                t.vessel_radius_mm
            );
        }
    }

    #[test]
    fn test_init_boundary_larger_than_vessel() {
        let (vol, spacing, origin) = make_vessel_phantom();
        let frame = make_frame(200);
        let params = AnnotationBatchParams::default();
        let targets = generate_annotation_batch(&frame, &vol, spacing, origin, &crate::types::IDENTITY_DIRECTION, &params);

        let mm_per_pixel = 2.0 * params.width_mm / params.pixels as f64;
        let center_px = params.pixels as f64 / 2.0;

        for (i, t) in targets.iter().enumerate() {
            // Average radius of init_boundary in pixels from center.
            let avg_init_r_px: f64 = t
                .init_boundary
                .iter()
                .map(|p| {
                    let dx = p[0] - center_px;
                    let dy = p[1] - center_px;
                    (dx * dx + dy * dy).sqrt()
                })
                .sum::<f64>()
                / t.init_boundary.len() as f64;

            let vessel_r_px = t.vessel_radius_mm / mm_per_pixel;

            assert!(
                avg_init_r_px > vessel_r_px,
                "Target {}: init_boundary avg radius ({:.1} px) should be > vessel radius ({:.1} px)",
                i,
                avg_init_r_px,
                vessel_r_px
            );
        }
    }

    #[test]
    fn test_image_dimensions_correct() {
        let (vol, spacing, origin) = make_vessel_phantom();
        let frame = make_frame(200);
        let params = AnnotationBatchParams::default();
        let targets = generate_annotation_batch(&frame, &vol, spacing, origin, &crate::types::IDENTITY_DIRECTION, &params);

        for (i, t) in targets.iter().enumerate() {
            assert_eq!(
                t.image.len(),
                t.pixels * t.pixels,
                "Target {}: image length {} != pixels^2 = {}",
                i,
                t.image.len(),
                t.pixels * t.pixels
            );
            assert_eq!(t.pixels, params.pixels);
        }
    }

    #[test]
    fn test_resolve_start_arc_none_returns_zero() {
        // No ostium → legacy behaviour (start at 0).
        let frame = make_frame(64);
        assert_eq!(resolve_start_arc(&frame, None, 2.0), 0.0);
    }

    #[test]
    fn test_resolve_start_arc_snaps_to_nearest_column() {
        // Straight centerline along +z axis from z=0 to z=40 (at y=0, x=0).
        // CprFrame stores positions as [z, y, x] in mm.
        let points: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [40.0, 0.0, 0.0]];
        let frame = CprFrame::from_centerline(&points, 64);

        // Ostium at z=5 mm (in [z, y, x] order). The nearest frame column's
        // arc-length should be ≈ 5 mm.
        let start = resolve_start_arc(&frame, Some([5.0, 0.0, 0.0]), 2.0);

        // Frame has ~64 columns spanning ~40 mm, so arc quantum ≈ 0.625 mm.
        // Expect snap to ≈ 5 mm within 1.5 × quantum tolerance.
        let arc_per_col = 40.0 / 64.0;
        assert!(
            (start - 5.0).abs() < arc_per_col * 1.5,
            "resolve_start_arc returned {}, expected ≈5 (tolerance {:.3})",
            start,
            arc_per_col * 1.5,
        );
    }

    #[test]
    fn test_resolve_start_arc_clamps_past_end() {
        // Ostium placed far downstream of the centerline end.
        // Should clamp to (total_arc − section_spacing_mm).
        let points: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [40.0, 0.0, 0.0]];
        let frame = CprFrame::from_centerline(&points, 64);

        // An ostium at z=1000 mm maps to the last frame column (farthest
        // reachable), whose arc-length is < total_arc. The clamp should never
        // produce a value > total_arc − section_spacing_mm.
        let start = resolve_start_arc(&frame, Some([1000.0, 0.0, 0.0]), 2.0);
        let total_arc = frame.arclengths[frame.n_cols() - 1];
        assert!(
            start <= total_arc - 2.0 + 1e-9,
            "start {} should be clamped ≤ total_arc − spacing = {}",
            start,
            total_arc - 2.0,
        );
        assert!(start >= 0.0, "start should be non-negative, got {}", start);
    }

    #[test]
    fn test_batch_starts_at_ostium() {
        // Full batch: check the first target's arc_mm reflects the ostium.
        let (vol, spacing, origin) = make_vessel_phantom();
        let frame = make_frame(200);

        // make_frame builds a centerline at y=32, x=32, z in [2, 62).
        // Put the ostium at z=10 mm (patient mm), in [z, y, x] order.
        let params = AnnotationBatchParams {
            n_sections: 5,
            section_spacing_mm: 2.0,
            ostium_zyx: Some([10.0, 32.0, 32.0]),
            ..AnnotationBatchParams::default()
        };
        let targets = generate_annotation_batch(
            &frame,
            &vol,
            spacing,
            origin,
            &crate::types::IDENTITY_DIRECTION,
            &params,
        );
        assert!(!targets.is_empty(), "expected ≥1 target when ostium set");

        // The first target's arc_mm should land near the ostium's arc-length,
        // which corresponds to ~(z=10 − centerline_start_z=2) = 8 mm along the
        // straight centerline, within one frame column quantum.
        let total_arc = frame.arclengths[frame.n_cols() - 1];
        let arc_per_col = total_arc / frame.n_cols() as f64;
        let first_arc = targets[0].arc_mm;
        assert!(
            (first_arc - 8.0).abs() < arc_per_col * 2.0,
            "first arc_mm = {}, expected ≈8 (tolerance {:.3})",
            first_arc,
            arc_per_col * 2.0,
        );
    }
}
