use std::collections::HashMap;

use ndarray::Array3;
use serde::Serialize;

use crate::annotation::AnnotationTarget;
use crate::cpr::CprFrame;
use crate::interp::trilinear;

/// Surface data for one cross-section: material values sampled on a polar (theta, r) grid.
#[derive(Clone, Serialize)]
pub struct CrossSectionSurface {
    /// Arc-length position along the centerline in mm.
    pub arc_mm: f64,
    /// Angular bin centers in degrees (e.g. [0, 22.5, 45, ..., 337.5]).
    pub theta_deg: Vec<f64>,
    /// Radial positions in mm from vessel wall.
    pub r_mm: Vec<f64>,
    /// Material values on the (theta, r) grid: [n_theta x n_radial].
    /// Row-major: surface[i_theta * n_radial + i_r].
    pub surface: Vec<f32>,
    /// Number of angular bins.
    pub n_theta: usize,
    /// Number of radial steps.
    pub n_radial: usize,
    /// Maximum radius per angular bin (outer contour boundary), in mm.
    pub max_r_per_theta: Vec<f64>,
}

/// Parameters for radial-angular sampling.
#[derive(Debug, Clone)]
pub struct RadialAngularParams {
    /// Number of angular bins (default: 16).
    pub n_theta: usize,
    /// Radial step size in mm (default: 0.5).
    pub radial_step_mm: f64,
    /// Maximum radial distance from wall in mm (default: 20.0).
    pub max_radius_mm: f64,
}

impl Default for RadialAngularParams {
    fn default() -> Self {
        Self {
            // Moderately dense sampling (7.5° steps, 0.4 mm rings → ~48×25) so
            // the surface reads as a smooth term-structure rather than the coarse
            // 16×20 lattice, without paying for the 72×40 grid most of whose
            // detail the viewer's smoother erased anyway. The viewer denoises in
            // PHYSICAL units (mostly angular, to preserve the thin radial
            // fat↔muscle boundary), so it works from whatever grid this provides
            // and does not depend on a specific cell count.
            n_theta: 48,
            radial_step_mm: 0.4,
            // 10 mm outward from the lumen wall covers the pericoronary
            // fat band (typically < 6 mm) plus a margin for context.
            max_radius_mm: 10.0,
        }
    }
}

/// Sample material values on polar grids around cross-section centerline positions.
///
/// For each annotation target, samples the specified material map at (theta, r)
/// grid points in the cross-section plane. r=0 sits on the user's finalized
/// contour (auto-adopted lumen wall by default) and r extends outward to
/// `params.max_radius_mm`. Samples that fall outside the CT volume are
/// returned as NaN (from `trilinear`). There is no independent "outer
/// boundary" any more — the r axis is a fixed span so sections can be
/// compared at the same radial depth.
///
/// # Arguments
///
/// * `material_map` - 3D volume of material values (e.g. lipid_frac or lipid_mass from MmdResult)
/// * `frame` - Precomputed CprFrame with Bishop frame geometry
/// * `annotation_targets` - The annotation targets (for vessel wall + contour info)
/// * `finalized_contours` - Map of target_index -> finalized contour points [x,y] pixel coords
/// * `spacing` - [sz, sy, sx] voxel spacing in mm
/// * `origin` - [oz, oy, ox] volume origin in mm
/// * `params` - Sampling parameters
pub fn sample_radial_angular(
    material_map: &Array3<f32>,
    frame: &CprFrame,
    annotation_targets: &[AnnotationTarget],
    finalized_contours: &HashMap<usize, Vec<[f64; 2]>>,
    spacing: [f64; 3],
    origin: [f64; 3],
    direction: &[f64; 9],
    params: &RadialAngularParams,
) -> Vec<CrossSectionSurface> {
    let n_radial = (params.max_radius_mm / params.radial_step_mm).ceil() as usize;
    let r_mm: Vec<f64> = (0..n_radial)
        .map(|i| i as f64 * params.radial_step_mm)
        .collect();
    let theta_deg: Vec<f64> = (0..params.n_theta)
        .map(|i| 360.0 * i as f64 / params.n_theta as f64)
        .collect();

    let inv_spacing = [1.0 / spacing[0], 1.0 / spacing[1], 1.0 / spacing[2]];

    let mut surfaces = Vec::new();

    for (target_idx, target) in annotation_targets.iter().enumerate() {
        // Only process targets that have a user-finalized contour — that's
        // the lumen-wall reference we sample outward from.
        let wall_contour = match finalized_contours.get(&target_idx) {
            Some(c) if !c.is_empty() => c,
            _ => continue,
        };

        let frame_idx = target.frame_index;
        if frame_idx >= frame.n_cols() {
            continue;
        }

        let center_px = target.pixels as f64 / 2.0;
        let mm_per_pixel = 2.0 * target.width_mm / target.pixels as f64;

        // Get frame vectors at this position.
        let pos_mm = frame.positions[frame_idx];
        let normal = frame.normals[frame_idx];
        let binormal = frame.binormals[frame_idx];

        // max_r_per_theta now records max_radius_mm uniformly since r extends
        // to a fixed outward distance. Kept for backward compatibility with
        // consumers that used it as a per-angle mask cue.
        let mut max_r_per_theta_vec = Vec::with_capacity(params.n_theta);
        let mut surface_data = vec![f32::NAN; params.n_theta * n_radial];

        for (i_theta, &theta_d) in theta_deg.iter().enumerate() {
            let theta_rad = theta_d.to_radians();
            let cos_t = theta_rad.cos();
            let sin_t = theta_rad.sin();

            // Ray direction in pixel coords from center: (cos, -sin) because y-axis is inverted.
            let ray_dx = cos_t;
            let ray_dy = -sin_t;

            // Distance from image center to the wall contour along this ray.
            // Convert to mm so r is measured in physical units outward from
            // the wall.
            let r_wall_px = ray_contour_intersection(
                center_px, center_px, ray_dx, ray_dy, wall_contour,
            );
            let r_wall_mm = r_wall_px * mm_per_pixel;

            max_r_per_theta_vec.push(params.max_radius_mm);

            // Direction in 3D world coords for this angle in the cross-section plane.
            // The cross-section image maps: row -> normal direction, col -> binormal direction.
            // Moving from center at angle theta in pixel space:
            //   pixel offset = (cos_t, -sin_t) -> (dx_col, dy_row)
            //   world offset = dy_row * (-normal) + dx_col * (-binormal)
            // Since row increases downward (decreasing normal offset)
            // and col increases rightward (decreasing binormal offset).
            //
            // In the cross-section image (see render_cross_section):
            //   row direction: top = +normal, bottom = -normal
            //   col direction: left = +binormal, right = -binormal
            // So offset_n = width_mm * (1 - 2*row/(pixels-1))  -> row=0 => +width_mm
            //    offset_b = width_mm * (1 - 2*col/(pixels-1))  -> col=0 => +width_mm
            //
            // A pixel displacement (dx_px, dy_px) from center means:
            //   d_row = dy_px (positive dy_px = lower row = more negative normal)
            //   d_col = dx_px (positive dx_px = higher col = more negative binormal)
            // So the world direction per pixel is:
            //   world_dir = (-dy_px * normal_hat - dx_px * binormal_hat) * mm_per_pixel
            // But we want unit direction in mm, so just normalize:
            let dir_3d = nalgebra::Vector3::new(
                (-ray_dy) * normal[0] + (-ray_dx) * binormal[0],
                (-ray_dy) * normal[1] + (-ray_dx) * binormal[1],
                (-ray_dy) * normal[2] + (-ray_dx) * binormal[2],
            );
            // dir_3d has magnitude 1 since normal and binormal are orthonormal and
            // (-ray_dy)^2 + (-ray_dx)^2 = sin^2 + cos^2 = 1.

            // Sample along the radial direction, from the wall outward.
            for (i_r, &r) in r_mm.iter().enumerate() {
                // World position: center + (r_wall_mm + r) along the radial direction.
                let dist_from_center_mm = r_wall_mm + r;
                let world_z = pos_mm[0] + dist_from_center_mm * dir_3d[0];
                let world_y = pos_mm[1] + dist_from_center_mm * dir_3d[1];
                let world_x = pos_mm[2] + dist_from_center_mm * dir_3d[2];

                // Convert to voxel coords (honors DICOM IOP).
                let [vz, vy, vx] = crate::types::patient_to_voxel(
                    [world_z, world_y, world_x],
                    origin,
                    inv_spacing,
                    direction,
                );

                let val = trilinear(material_map, vz, vy, vx);
                surface_data[i_theta * n_radial + i_r] = val;
            }
        }

        surfaces.push(CrossSectionSurface {
            arc_mm: target.arc_mm,
            theta_deg: theta_deg.clone(),
            r_mm: r_mm.clone(),
            surface: surface_data,
            n_theta: params.n_theta,
            n_radial,
            max_r_per_theta: max_r_per_theta_vec,
        });
    }

    surfaces
}

/// Find the distance (in pixels) from a ray origin to the first intersection with a polygon contour.
///
/// The ray starts at (ox, oy) and goes in direction (dx, dy).
/// The contour is a closed polygon given as [x, y] pixel coordinates.
///
/// Returns the parametric distance `t` (in pixels) of the closest intersection,
/// or `f64::MAX` if no intersection is found.
fn ray_contour_intersection(
    ox: f64,
    oy: f64,
    dx: f64,
    dy: f64,
    contour: &[[f64; 2]],
) -> f64 {
    let n = contour.len();
    if n < 3 {
        return f64::MAX;
    }

    let mut best_t = f64::MAX;

    for i in 0..n {
        let j = (i + 1) % n;

        let p0 = contour[i];
        let p1 = contour[j];

        // Edge direction.
        let ex = p1[0] - p0[0];
        let ey = p1[1] - p0[1];

        // Solve: origin + t * dir = p0 + s * edge
        // [dx, -ex] [t]   [p0x - ox]
        // [dy, -ey] [s] = [p0y - oy]
        let denom = dx * (-ey) - dy * (-ex);
        if denom.abs() < 1e-12 {
            continue; // Parallel.
        }

        let bx = p0[0] - ox;
        let by = p0[1] - oy;

        let t = (bx * (-ey) - by * (-ex)) / denom;
        let s = (dx * by - dy * bx) / denom;

        if t > 1e-6 && s >= 0.0 && s <= 1.0 && t < best_t {
            best_t = t;
        }
    }

    best_t
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::annotation::AnnotationTarget;
    use crate::cpr::CprFrame;
    use ndarray::Array3;
    use std::collections::HashMap;

    /// Create a simple z-axis centerline frame.
    fn make_simple_frame(n_cols: usize) -> CprFrame {
        let points: Vec<[f64; 3]> = (0..60)
            .map(|z| [z as f64 + 2.0, 32.0, 32.0])
            .collect();
        CprFrame::from_centerline(&points, n_cols)
    }

    /// Create a circular contour centered at (cx, cy) with radius r in pixel coords.
    fn circular_contour(cx: f64, cy: f64, r: f64, n_points: usize) -> Vec<[f64; 2]> {
        (0..n_points)
            .map(|i| {
                let theta = 2.0 * std::f64::consts::PI * i as f64 / n_points as f64;
                [cx + r * theta.cos(), cy - r * theta.sin()]
            })
            .collect()
    }

    /// Build a synthetic annotation target with circular vessel wall and given parameters.
    fn make_target(frame_index: usize, arc_mm: f64, pixels: usize, width_mm: f64) -> AnnotationTarget {
        let center_px = pixels as f64 / 2.0;
        let mm_per_pixel = 2.0 * width_mm / pixels as f64;
        let vessel_radius_mm = 3.0;
        let vessel_radius_px = vessel_radius_mm / mm_per_pixel;

        let vessel_wall = circular_contour(center_px, center_px, vessel_radius_px, 72);
        let init_radius_px = 2.0 * vessel_radius_px;
        let init_boundary = circular_contour(center_px, center_px, init_radius_px, 72);

        AnnotationTarget {
            image: vec![0.0f32; pixels * pixels],
            pixels,
            width_mm,
            arc_mm,
            frame_index,
            vessel_wall,
            vessel_radius_mm,
            init_boundary,
        }
    }

    #[test]
    fn test_correct_surface_dimensions() {
        // Uniform volume.
        let vol = Array3::<f32>::from_elem((64, 64, 64), 0.5);
        let spacing = [1.0, 1.0, 1.0];
        let origin = [0.0, 0.0, 0.0];

        let frame = make_simple_frame(100);
        let params = RadialAngularParams {
            n_theta: 16,
            radial_step_mm: 0.5,
            max_radius_mm: 10.0,
        };
        let n_radial_expected = (params.max_radius_mm / params.radial_step_mm).ceil() as usize;

        let target = make_target(50, 25.0, 128, 15.0);
        let targets = vec![target];

        // Finalized contour: circle larger than vessel wall.
        let center_px = 64.0;
        let mm_per_pixel = 2.0 * 15.0 / 128.0;
        let outer_radius_px = 12.0 / mm_per_pixel; // 12mm outer radius
        let outer_contour = circular_contour(center_px, center_px, outer_radius_px, 72);

        let mut finalized = HashMap::new();
        finalized.insert(0, outer_contour);

        let surfaces = sample_radial_angular(
            &vol, &frame, &targets, &finalized, spacing, origin,
            &crate::types::IDENTITY_DIRECTION, &params,
        );

        assert_eq!(surfaces.len(), 1, "Should produce 1 surface");
        let s = &surfaces[0];
        assert_eq!(s.n_theta, 16);
        assert_eq!(s.n_radial, n_radial_expected);
        assert_eq!(s.theta_deg.len(), 16);
        assert_eq!(s.r_mm.len(), n_radial_expected);
        assert_eq!(s.surface.len(), 16 * n_radial_expected);
        assert_eq!(s.max_r_per_theta.len(), 16);
    }

    /// Consistency check (RIGHT path): the real `sample_radial_angular` must
    /// sample the geometrically-correct world point for each (theta, r).
    /// We encode world-x into the volume (value = voxel x; identity direction,
    /// unit spacing, zero origin -> trilinear of a linear ramp is exact), so the
    /// sampled value should equal the world-x of the point the sampler claims to
    /// hit: pos.x + (r_wall_mm + r) * dir_x, with dir = sinθ·normal - cosθ·binormal.
    #[test]
    fn radial_sampling_lands_at_expected_world_position() {
        let (nz, ny, nx) = (80usize, 80, 80);
        let mut vol = Array3::<f32>::zeros((nz, ny, nx));
        for ((_z, _y, x), v) in vol.indexed_iter_mut() {
            *v = x as f32; // value == world-x (identity dir, unit spacing, zero origin)
        }
        let spacing = [1.0, 1.0, 1.0];
        let origin = [0.0, 0.0, 0.0];

        let frame = make_simple_frame(100);
        let fi = 50usize;
        let pos = frame.positions[fi];
        let normal = frame.normals[fi];
        let binormal = frame.binormals[fi];

        let pixels = 128usize;
        let width_mm = 15.0_f64;
        let mm_per_pixel = 2.0 * width_mm / pixels as f64;
        let r_wall_mm = 4.0_f64;
        let center_px = pixels as f64 / 2.0;
        let contour = circular_contour(center_px, center_px, r_wall_mm / mm_per_pixel, 144);
        let mut finalized = HashMap::new();
        finalized.insert(0, contour);

        let target = make_target(fi, 25.0, pixels, width_mm);
        let targets = vec![target];
        let params = RadialAngularParams { n_theta: 8, radial_step_mm: 1.0, max_radius_mm: 8.0 };

        let surfaces = sample_radial_angular(
            &vol, &frame, &targets, &finalized, spacing, origin,
            &crate::types::IDENTITY_DIRECTION, &params,
        );
        let s = &surfaces[0];

        for (it, &theta_d) in s.theta_deg.iter().enumerate() {
            let th = theta_d.to_radians();
            let dir_x = th.sin() * normal[2] - th.cos() * binormal[2];
            for (ir, &r) in s.r_mm.iter().enumerate() {
                let got = s.surface[it * s.n_radial + ir];
                if got.is_nan() { continue; }
                let expected_world_x = pos[2] + (r_wall_mm + r) * dir_x;
                assert!(
                    (got as f64 - expected_world_x).abs() < 1e-2,
                    "θ={theta_d} r={r}: sampler hit world-x {got} but geometry says {expected_world_x}"
                );
            }
        }
    }

    /// Consistency check (LEFT vs RIGHT plane): the 2D overlay
    /// (`get_mmd_overlay`, which maps pixel (row,col) -> world) and the 3D
    /// radial surface (`sample_radial_angular`, which maps (theta, dist) ->
    /// world) MUST place the same physical point at the same world position,
    /// or the two MMD views show misregistered data. Both reduce to
    /// world = pos + off_n·normal + off_b·binormal; this pins that equivalence.
    /// Also documents the one known sub-pixel convention difference: the overlay
    /// uses pixel→mm = 2·w/(pixels-1) while the radial path uses 2·w/pixels
    /// (~0.4–0.8% radial scale; NOT the cause of the LEFT/RIGHT colour mismatch,
    /// which is region + resolution + smoothing + per-pixel noise).
    #[test]
    fn overlay_and_radial_agree_on_world_position() {
        let frame = make_simple_frame(100);
        let fi = 50usize;
        let pos = frame.positions[fi];
        let normal = frame.normals[fi];
        let binormal = frame.binormals[fi];
        let pixels = 128usize;
        let width_mm = 15.0_f64;

        let world_from_offset = |off_n: f64, off_b: f64| {
            [
                pos[0] + off_n * normal[0] + off_b * binormal[0],
                pos[1] + off_n * normal[1] + off_b * binormal[1],
                pos[2] + off_n * normal[2] + off_b * binormal[2],
            ]
        };

        for &(off_n, off_b) in &[(4.3, -2.1), (-6.0, 0.0), (0.0, 7.5), (2.0, 2.0)] {
            // Overlay convention: off_n = w(1 - 2·row/(p-1)) -> recover (row,col),
            // then world = pos + off_n·normal + off_b·binormal.
            let row = (1.0 - off_n / width_mm) * (pixels as f64 - 1.0) / 2.0;
            let col = (1.0 - off_b / width_mm) * (pixels as f64 - 1.0) / 2.0;
            let ov_off_n = width_mm * (1.0 - 2.0 * row / (pixels as f64 - 1.0));
            let ov_off_b = width_mm * (1.0 - 2.0 * col / (pixels as f64 - 1.0));
            let w_overlay = world_from_offset(ov_off_n, ov_off_b);

            // Radial convention: off_n = dist·sinθ, off_b = -dist·cosθ;
            // world = pos + dist·(sinθ·normal - cosθ·binormal).
            let dist = (off_n * off_n + off_b * off_b).sqrt();
            let theta = off_n.atan2(-off_b);
            let (sn, cs) = (theta.sin(), theta.cos());
            let w_radial = [
                pos[0] + dist * (sn * normal[0] - cs * binormal[0]),
                pos[1] + dist * (sn * normal[1] - cs * binormal[1]),
                pos[2] + dist * (sn * normal[2] - cs * binormal[2]),
            ];

            for k in 0..3 {
                assert!(
                    (w_overlay[k] - w_radial[k]).abs() < 1e-9,
                    "off=({off_n},{off_b}) axis {k}: overlay {} vs radial {}",
                    w_overlay[k], w_radial[k]
                );
            }
        }

        // Documented sub-pixel convention difference (the only geometric
        // inconsistency between the two views; well under 1%).
        let mmpp_overlay = 2.0 * width_mm / (pixels as f64 - 1.0);
        let mmpp_radial = 2.0 * width_mm / pixels as f64;
        let rel = (mmpp_overlay - mmpp_radial).abs() / mmpp_radial;
        assert!(rel < 0.01, "pixel→mm convention diff {rel} unexpectedly large");
    }

    #[test]
    fn test_uniform_volume_produces_uniform_surface() {
        let val = 0.5f32;
        let vol = Array3::<f32>::from_elem((64, 64, 64), val);
        let spacing = [1.0, 1.0, 1.0];
        let origin = [0.0, 0.0, 0.0];

        let frame = make_simple_frame(100);
        let params = RadialAngularParams {
            n_theta: 8,
            radial_step_mm: 0.5,
            max_radius_mm: 5.0,
        };

        let target = make_target(50, 25.0, 128, 15.0);
        let targets = vec![target];

        let center_px = 64.0;
        let mm_per_pixel = 2.0 * 15.0 / 128.0;
        let outer_radius_px = 10.0 / mm_per_pixel;
        let outer_contour = circular_contour(center_px, center_px, outer_radius_px, 72);

        let mut finalized = HashMap::new();
        finalized.insert(0, outer_contour);

        let surfaces = sample_radial_angular(
            &vol, &frame, &targets, &finalized, spacing, origin,
            &crate::types::IDENTITY_DIRECTION, &params,
        );

        assert_eq!(surfaces.len(), 1);
        let s = &surfaces[0];

        // Check non-NaN values are approximately 0.5.
        let mut non_nan_count = 0;
        for &v in &s.surface {
            if !v.is_nan() {
                assert!(
                    (v - val).abs() < 0.01,
                    "Expected ~{}, got {}",
                    val,
                    v
                );
                non_nan_count += 1;
            }
        }
        assert!(
            non_nan_count > 0,
            "Should have at least some non-NaN samples"
        );
    }

    #[test]
    fn test_surface_values_bounded() {
        // Volume with values in [0, 1].
        let mut vol = Array3::<f32>::zeros((64, 64, 64));
        for ((z, y, x), v) in vol.indexed_iter_mut() {
            *v = ((z + y + x) as f32 / 192.0).clamp(0.0, 1.0);
        }
        let spacing = [1.0, 1.0, 1.0];
        let origin = [0.0, 0.0, 0.0];

        let frame = make_simple_frame(100);
        let params = RadialAngularParams {
            n_theta: 8,
            radial_step_mm: 0.5,
            max_radius_mm: 5.0,
        };

        let target = make_target(50, 25.0, 128, 15.0);
        let targets = vec![target];

        let center_px = 64.0;
        let mm_per_pixel = 2.0 * 15.0 / 128.0;
        let outer_radius_px = 10.0 / mm_per_pixel;
        let outer_contour = circular_contour(center_px, center_px, outer_radius_px, 72);

        let mut finalized = HashMap::new();
        finalized.insert(0, outer_contour);

        let surfaces = sample_radial_angular(
            &vol, &frame, &targets, &finalized, spacing, origin,
            &crate::types::IDENTITY_DIRECTION, &params,
        );

        assert_eq!(surfaces.len(), 1);
        let s = &surfaces[0];

        for &v in &s.surface {
            if !v.is_nan() {
                assert!(
                    v >= -0.01 && v <= 1.01,
                    "Value {} out of expected [0, 1] range",
                    v
                );
            }
        }
    }

    #[test]
    fn test_r_axis_uniform_across_contours() {
        // The r axis is a fixed grid [0, step, 2*step, ..., max_radius_mm];
        // it does not depend on the contour shape. Two different finalized
        // contours must produce the same r_mm and max_r_per_theta.
        let vol = Array3::<f32>::from_elem((64, 64, 64), 0.5);
        let spacing = [1.0, 1.0, 1.0];
        let origin = [0.0, 0.0, 0.0];

        let frame = make_simple_frame(100);
        let params = RadialAngularParams {
            n_theta: 8,
            radial_step_mm: 0.5,
            max_radius_mm: 4.0,
        };

        let target = make_target(50, 25.0, 128, 15.0);
        let targets = vec![target];

        let center_px = 64.0;
        let mm_per_pixel = 2.0 * 15.0 / 128.0;
        let small = circular_contour(center_px, center_px, 3.0 / mm_per_pixel, 72);
        let large = circular_contour(center_px, center_px, 5.0 / mm_per_pixel, 72);

        let mut finalized_small = HashMap::new();
        finalized_small.insert(0, small);
        let mut finalized_large = HashMap::new();
        finalized_large.insert(0, large);

        let s_small = &sample_radial_angular(
            &vol, &frame, &targets, &finalized_small, spacing, origin,
            &crate::types::IDENTITY_DIRECTION, &params,
        )[0];
        let s_large = &sample_radial_angular(
            &vol, &frame, &targets, &finalized_large, spacing, origin,
            &crate::types::IDENTITY_DIRECTION, &params,
        )[0];

        assert_eq!(s_small.r_mm, s_large.r_mm);
        assert_eq!(s_small.max_r_per_theta, s_large.max_r_per_theta);
        for &m in &s_small.max_r_per_theta {
            assert!((m - params.max_radius_mm).abs() < 1e-9);
        }
    }

    #[test]
    fn test_no_finalized_contour_skips_target() {
        let vol = Array3::<f32>::from_elem((64, 64, 64), 0.5);
        let spacing = [1.0, 1.0, 1.0];
        let origin = [0.0, 0.0, 0.0];

        let frame = make_simple_frame(100);
        let params = RadialAngularParams::default();

        let target = make_target(50, 25.0, 128, 15.0);
        let targets = vec![target];

        // Empty finalized map: no contours finalized.
        let finalized = HashMap::new();

        let surfaces = sample_radial_angular(
            &vol, &frame, &targets, &finalized, spacing, origin,
            &crate::types::IDENTITY_DIRECTION, &params,
        );

        assert_eq!(surfaces.len(), 0, "No surfaces if no finalized contours");
    }

    #[test]
    fn test_ray_contour_intersection_circle() {
        // Verify ray-contour intersection with a known circular contour.
        let cx = 64.0;
        let cy = 64.0;
        let r = 20.0;
        let contour = circular_contour(cx, cy, r, 360);

        // Ray from center going right (positive x, angle=0).
        let t = ray_contour_intersection(cx, cy, 1.0, 0.0, &contour);
        assert!(
            (t - r).abs() < 1.0,
            "Expected intersection at ~{} px, got {} px",
            r,
            t
        );

        // Ray from center going up (angle=90, dy=-1 because y is inverted in pixel space,
        // but our ray direction is (cos, -sin) = (0, -1)).
        let t_up = ray_contour_intersection(cx, cy, 0.0, -1.0, &contour);
        assert!(
            (t_up - r).abs() < 1.0,
            "Expected upward intersection at ~{} px, got {} px",
            r,
            t_up
        );
    }
}
