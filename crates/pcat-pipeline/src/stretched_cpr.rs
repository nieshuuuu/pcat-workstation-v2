//! Horos-style Stretched CPR renderer.
//!
//! Each output column corresponds to one projected arc-length position along the
//! centerline, with the column world position anchored to the mid-height plane.
//! Vertical pixels step along `projection_normal`. Horizontal pixel spacing is
//! `total_proj_arc / (pixels_wide - 1)` -- isotropic by construction.

use nalgebra::Vector3;
use ndarray::Array3;
use rayon::prelude::*;

use crate::interp::trilinear;

// ---------------------------------------------------------------------------
// Public result type
// ---------------------------------------------------------------------------

pub struct StretchedCprResult {
    /// Flattened row-major image, shape (pixels_high, pixels_wide).
    /// Pixels outside the volume bounds are NAN.
    pub image: Vec<f32>,
    pub pixels_wide: usize,
    pub pixels_high: usize,
    /// Pass-through of the input arc-lengths (one per CprFrame column).
    pub arclengths: Vec<f64>,
}

// ---------------------------------------------------------------------------
// Geometry computation (shared by renderer + get_cpr_projection_info)
// ---------------------------------------------------------------------------

/// All pre-computed geometry needed to render one stretched CPR frame.
pub struct StretchedGeometry {
    pub projection_normal: Vector3<f64>,
    pub mid_height_point: Vector3<f64>,
    /// One point per output column on the mid-height plane (3D).
    pub proj_col_pts: Vec<Vector3<f64>>,
    /// One point per output column on the original 3D centerline.
    pub orig_col_pts: Vec<Vector3<f64>>,
    /// In-plane slab direction per column.
    pub slab_dirs: Vec<Vector3<f64>>,
    /// Total projected arc-length (horizontal extent).
    pub total_proj_arc: f64,
    /// Isotropic pixel spacing (mm/pixel).
    pub dy_mm: f64,
    /// Signed projection of all centerline positions onto `projection_normal`;
    /// `max - min` is the vessel's out-of-plane depth span and drives
    /// adaptive `pixels_high` sizing so the vessel does not clip vertically
    /// at oblique rotations.
    pub proj_min: f64,
    pub proj_max: f64,
}

/// Compute the stretched-CPR geometry from a set of pre-sampled centerline positions.
///
/// `rotation_deg` rotates `base_normal` around `curve_direction` via Rodrigues formula.
pub fn compute_stretched_geometry(
    positions: &[[f64; 3]],
    pixels_wide: usize,
    rotation_deg: f64,
) -> StretchedGeometry {
    let n = positions.len();
    assert!(n >= 2, "need at least 2 centerline positions");
    assert!(pixels_wide >= 2, "pixels_wide must be at least 2");

    let pos_vecs: Vec<Vector3<f64>> = positions.iter()
        .map(|p| Vector3::new(p[0], p[1], p[2]))
        .collect();

    // 1. curve_direction
    let curve_direction = {
        let v = pos_vecs[n - 1] - pos_vecs[0];
        if v.norm() > 1e-9 {
            v.normalize()
        } else {
            let v2 = pos_vecs[n - 1] - pos_vecs[n / 2];
            if v2.norm() > 1e-9 { v2.normalize() }
            else { Vector3::new(1.0, 0.0, 0.0) }
        }
    };

    // 2. base_direction (world Z = patient SI); fall back to world Y if nearly parallel
    let base_direction = {
        let candidate = Vector3::new(1.0, 0.0, 0.0);
        if candidate.dot(&curve_direction).abs() > 0.99 {
            Vector3::new(0.0, 1.0, 0.0)
        } else {
            candidate
        }
    };
    let base_normal = base_direction.cross(&curve_direction).normalize();

    // 3. projection_normal via Rodrigues rotation around curve_direction
    let theta = rotation_deg.to_radians();
    let projection_normal = {
        let k = curve_direction;
        let cos_t = theta.cos();
        let sin_t = theta.sin();
        let v = base_normal;
        // Rodrigues: v cos(t) + (k x v) sin(t) + k (k.v)(1 - cos(t))
        (v * cos_t + k.cross(&v) * sin_t + k * k.dot(&v) * (1.0 - cos_t)).normalize()
    };
    // Invariant: projection_normal perp to curve_direction (within fp precision)
    debug_assert!(
        projection_normal.dot(&curve_direction).abs() < 1e-9,
        "projection_normal . curve_direction = {}",
        projection_normal.dot(&curve_direction)
    );

    // 4. mid_height_point
    let centroid = pos_vecs.iter().fold(Vector3::zeros(), |acc, p| acc + p) / n as f64;
    let centroid_dot = centroid.dot(&projection_normal);
    let projections: Vec<f64> = pos_vecs.iter().map(|p| p.dot(&projection_normal)).collect();
    let min_proj = projections.iter().cloned().fold(f64::INFINITY, f64::min);
    let max_proj = projections.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let mid_val = (min_proj + max_proj) / 2.0;
    let mid_height_point = centroid + (mid_val - centroid_dot) * projection_normal;

    // 5. Project each centerline point onto the mid-height plane and build cumulative arc
    let project_to_plane = |p: &Vector3<f64>| {
        p - (p - mid_height_point).dot(&projection_normal) * projection_normal
    };
    let proj_pts: Vec<Vector3<f64>> = pos_vecs.iter().map(|p| project_to_plane(p)).collect();

    let mut proj_arc = vec![0.0f64; n];
    for i in 1..n {
        proj_arc[i] = proj_arc[i - 1] + (proj_pts[i] - proj_pts[i - 1]).norm();
    }
    let total_proj_arc = proj_arc[n - 1];
    assert!(total_proj_arc > 1e-12, "projected arc-length is zero -- degenerate centerline");

    // 6. Sample evenly-spaced column points along the projected arc
    let mut orig_col_pts = Vec::with_capacity(pixels_wide);
    let mut proj_col_pts_3d = Vec::with_capacity(pixels_wide);

    for j in 0..pixels_wide {
        let s = j as f64 / (pixels_wide - 1) as f64 * total_proj_arc;

        let seg_idx = match proj_arc.binary_search_by(|v| v.partial_cmp(&s).unwrap()) {
            Ok(idx) => idx.min(n - 2),
            Err(idx) => if idx == 0 { 0 } else { (idx - 1).min(n - 2) },
        };
        let i0 = seg_idx;
        let i1 = (seg_idx + 1).min(n - 1);
        let seg_len = proj_arc[i1] - proj_arc[i0];
        let frac = if seg_len > 1e-20 {
            ((s - proj_arc[i0]) / seg_len).clamp(0.0, 1.0)
        } else {
            0.0
        };

        orig_col_pts.push(pos_vecs[i0] + frac * (pos_vecs[i1] - pos_vecs[i0]));
        proj_col_pts_3d.push(proj_pts[i0] + frac * (proj_pts[i1] - proj_pts[i0]));
    }

    // 7. Slab direction per column: normalize(projection_normal x tangent_in_plane)
    let mut slab_dirs = Vec::with_capacity(pixels_wide);
    for j in 0..pixels_wide {
        let tangent_in_plane = if j + 1 < pixels_wide {
            let diff = proj_col_pts_3d[j + 1] - proj_col_pts_3d[j];
            if diff.norm() > 1e-12 { diff.normalize() } else { Vector3::zeros() }
        } else {
            let diff = proj_col_pts_3d[j] - proj_col_pts_3d[j - 1];
            if diff.norm() > 1e-12 { diff.normalize() } else { Vector3::zeros() }
        };

        let slab_cross = projection_normal.cross(&tangent_in_plane);
        slab_dirs.push(
            if slab_cross.norm() > 1e-12 { slab_cross.normalize() }
            else { projection_normal }
        );
    }

    let dy_mm = total_proj_arc / (pixels_wide - 1) as f64;

    StretchedGeometry {
        projection_normal,
        mid_height_point,
        proj_col_pts: proj_col_pts_3d,
        orig_col_pts,
        slab_dirs,
        total_proj_arc,
        dy_mm,
        proj_min: min_proj,
        proj_max: max_proj,
    }
}

// ---------------------------------------------------------------------------
// Renderer
// ---------------------------------------------------------------------------

#[allow(clippy::too_many_arguments)]
pub(crate) fn render_stretched(
    positions: &[[f64; 3]],
    _normals: &[Vector3<f64>],
    _binormals: &[Vector3<f64>],
    arclengths: &[f64],
    volume: &Array3<f32>,
    spacing: [f64; 3],
    origin: [f64; 3],
    direction: &[f64; 9],
    width_mm: f64,
    pixels_wide: usize,
    pixels_high: usize,
    slab_mm: f64,
    rotation_deg: f64,
) -> StretchedCprResult {
    let n = positions.len();
    assert!(n >= 2, "need at least 2 centerline positions");
    assert_eq!(arclengths.len(), n);
    assert!(pixels_wide >= 2 && pixels_high >= 2);

    let geom = compute_stretched_geometry(positions, pixels_wide, rotation_deg);
    let inv_spacing = [1.0 / spacing[0], 1.0 / spacing[1], 1.0 / spacing[2]];

    let n_slab_steps = if slab_mm > 0.01 { 5usize } else { 1 };
    let slab_offsets: Vec<f64> = if n_slab_steps > 1 {
        (0..n_slab_steps)
            .map(|k| -slab_mm / 2.0 + slab_mm * k as f64 / (n_slab_steps - 1) as f64)
            .collect()
    } else {
        vec![0.0]
    };

    let mut image = vec![f32::NAN; pixels_high * pixels_wide];

    image
        .par_chunks_mut(pixels_wide)
        .enumerate()
        .for_each(|(row, row_slice)| {
            // row 0 = top = +Y; row (pixels_high-1) = bottom = -Y
            let y_offset_mm = (pixels_high as f64 / 2.0 - row as f64) * geom.dy_mm;

            // Skip rows outside the lateral field of view
            if y_offset_mm.abs() > width_mm + geom.dy_mm {
                return;
            }

            for col in 0..pixels_wide {
                let base_pt = geom.proj_col_pts[col] + y_offset_mm * geom.projection_normal;
                let slab_dir = &geom.slab_dirs[col];

                let mut sum = 0.0f64;
                let mut count = 0u32;

                for &s_off in &slab_offsets {
                    let sample_pt = base_pt + s_off * slab_dir;
                    let [vz, vy, vx] = crate::types::patient_to_voxel(
                        [sample_pt[0], sample_pt[1], sample_pt[2]],
                        origin,
                        inv_spacing,
                        direction,
                    );
                    let val = trilinear(volume, vz, vy, vx);
                    if !val.is_nan() {
                        sum += val as f64;
                        count += 1;
                    }
                }

                row_slice[col] = if count > 0 {
                    (sum / count as f64) as f32
                } else {
                    f32::NAN
                };
            }
        });

    StretchedCprResult {
        image,
        pixels_wide,
        pixels_high,
        arclengths: arclengths.to_vec(),
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// projection_normal must be perpendicular to curve_direction at every rotation angle.
    #[test]
    fn test_projection_normal_perpendicular_to_curve() {
        let positions: Vec<[f64; 3]> = (0..50).map(|z| [z as f64, 32.0, 32.0]).collect();

        for rot in [0.0, 30.0, 45.0, 90.0, 135.0, 180.0, 270.0, 359.9f64] {
            let geom = compute_stretched_geometry(&positions, 64, rot);
            let curve_dir = (Vector3::new(positions[49][0], positions[49][1], positions[49][2])
                - Vector3::new(positions[0][0], positions[0][1], positions[0][2]))
                .normalize();
            let dot = geom.projection_normal.dot(&curve_dir);
            assert!(
                dot.abs() < 1e-9,
                "rot={rot}deg: projection_normal . curve_direction = {dot:.2e}"
            );
        }
    }

    /// Degenerate case: centerline along Y (base_direction fallback branch).
    #[test]
    fn test_projection_normal_perp_axis_aligned_input() {
        let positions: Vec<[f64; 3]> = (0..30).map(|y| [32.0, y as f64, 32.0]).collect();

        for rot in [0.0, 90.0, 180.0f64] {
            let geom = compute_stretched_geometry(&positions, 32, rot);
            let curve_dir = (Vector3::new(positions[29][0], positions[29][1], positions[29][2])
                - Vector3::new(positions[0][0], positions[0][1], positions[0][2]))
                .normalize();
            let dot = geom.projection_normal.dot(&curve_dir);
            assert!(dot.abs() < 1e-9, "Y-axis curve, rot={rot}deg: dot = {dot:.2e}");
        }
    }

    /// Straight Z-axis centerline in a z-gradient volume: mid-row values must be monotonic
    /// along the column axis, spanning most of the centerline's Z range.
    ///
    /// Volume: `vol[z, y, x] = z as f32` so that a straight centerline from z=10 to z=53
    /// (at y=32, x=32) samples increasing Z values as the column advances.  The mid-row
    /// must therefore be monotonic and its range must cover ≥80 % of z_end − z_start.
    #[test]
    fn test_stretched_z_gradient_monotonic() {
        let nz = 64usize;
        let ny = 64usize;
        let nx = 64usize;

        // vol[z, y, x] = z
        let vol = {
            let mut v = Array3::<f32>::zeros((nz, ny, nx));
            for iz in 0..nz {
                for iy in 0..ny {
                    for ix in 0..nx {
                        v[[iz, iy, ix]] = iz as f32;
                    }
                }
            }
            v
        };
        let spacing = [1.0f64, 1.0, 1.0];
        let origin = [0.0f64, 0.0, 0.0];

        let z_start = 10usize;
        let z_end = 53usize;
        let n_pts = z_end - z_start + 1;
        // Centerline along Z at (y=32, x=32)
        let positions: Vec<[f64; 3]> =
            (z_start..=z_end).map(|z| [z as f64, 32.0, 32.0]).collect();
        let arclengths: Vec<f64> = (0..n_pts).map(|i| i as f64).collect();

        let normals: Vec<Vector3<f64>> =
            (0..n_pts).map(|_| Vector3::new(0.0, 1.0, 0.0)).collect();
        let binormals: Vec<Vector3<f64>> =
            (0..n_pts).map(|_| Vector3::new(0.0, 0.0, 1.0)).collect();

        let pixels_wide = 64usize;
        let pixels_high = 64usize;

        let result = render_stretched(
            &positions, &normals, &binormals, &arclengths,
            &vol, spacing, origin, &crate::types::IDENTITY_DIRECTION,
            20.0, pixels_wide, pixels_high, 0.5, 0.0,
        );

        // Extract mid-row (non-NaN values only)
        let mid_row = result.pixels_high / 2;
        let mid_row_vals: Vec<f32> = (0..result.pixels_wide)
            .map(|c| result.image[mid_row * result.pixels_wide + c])
            .collect();

        let valid: Vec<f32> = mid_row_vals.iter().copied().filter(|v| !v.is_nan()).collect();
        assert!(
            valid.len() > 10,
            "expected many non-NaN pixels in mid-row, got {}",
            valid.len()
        );

        // Assert monotonicity: either all diffs ≥ 0 or all diffs ≤ 0
        let diffs: Vec<f32> = valid.windows(2).map(|w| w[1] - w[0]).collect();
        let increasing = diffs.iter().all(|&d| d >= -0.5);
        let decreasing = diffs.iter().all(|&d| d <= 0.5);
        assert!(
            increasing || decreasing,
            "mid-row values are not monotonic; first 10 values: {:?}",
            &valid[..valid.len().min(10)]
        );

        // Assert that the value range spans ≥80 % of the centerline's Z range
        let z_range = (z_end - z_start) as f32;
        let val_min = valid.iter().cloned().fold(f32::INFINITY, f32::min);
        let val_max = valid.iter().cloned().fold(f32::NEG_INFINITY, f32::max);
        let val_span = val_max - val_min;
        assert!(
            val_span >= 0.8 * z_range,
            "mid-row value span {val_span:.1} should be ≥80% of z_range {z_range:.1}"
        );
    }

    /// L-shape centerline: vessel (HU=300) tube should appear near mid-Y in most columns.
    #[test]
    fn test_l_shape_vessel_at_mid_row() {
        let leg1: Vec<[f64; 3]> = (0..20).map(|i| [20.0 + i as f64, 20.0, 32.0]).collect();
        let leg2: Vec<[f64; 3]> = (1..21).map(|i| [39.0, 20.0 + i as f64, 32.0]).collect();
        let positions: Vec<[f64; 3]> = leg1.into_iter().chain(leg2).collect();
        let n_pts = positions.len();

        let arclengths: Vec<f64> = {
            let mut s = vec![0.0f64; n_pts];
            for i in 1..n_pts {
                let dz = positions[i][0] - positions[i - 1][0];
                let dy = positions[i][1] - positions[i - 1][1];
                let dx = positions[i][2] - positions[i - 1][2];
                s[i] = s[i - 1] + (dz * dz + dy * dy + dx * dx).sqrt();
            }
            s
        };

        let vol_size = 64usize;
        let tube_radius = 2.0f64;
        let mut vol = Array3::<f32>::zeros((vol_size, vol_size, vol_size));
        for iz in 0..vol_size {
            for iy in 0..vol_size {
                for ix in 0..vol_size {
                    let vz = iz as f64;
                    let vy = iy as f64;
                    let vx = ix as f64;
                    let min_d = positions.iter().map(|p| {
                        ((vz - p[0]).powi(2) + (vy - p[1]).powi(2) + (vx - p[2]).powi(2)).sqrt()
                    }).fold(f64::INFINITY, f64::min);
                    if min_d < tube_radius {
                        vol[[iz, iy, ix]] = 300.0;
                    }
                }
            }
        }

        let normals: Vec<Vector3<f64>> = (0..n_pts).map(|_| Vector3::new(0.0, 1.0, 0.0)).collect();
        let binormals: Vec<Vector3<f64>> = (0..n_pts).map(|_| Vector3::new(0.0, 0.0, 1.0)).collect();

        let pixels_wide = 40usize;
        let pixels_high = 60usize;
        let result = render_stretched(
            &positions, &normals, &binormals, &arclengths,
            &vol, [1.0, 1.0, 1.0], [0.0, 0.0, 0.0], &crate::types::IDENTITY_DIRECTION,
            30.0, pixels_wide, pixels_high, 0.0, 0.0,
        );

        let geom = compute_stretched_geometry(&positions, pixels_wide, 0.0);
        assert!(geom.total_proj_arc > 0.0, "projected arc-length must be positive");

        let mid_row = pixels_high / 2;
        let mut n_col_checked = 0;
        let mut n_col_passed = 0;

        for col in 0..pixels_wide {
            let col_vals: Vec<(usize, f32)> = (0..pixels_high)
                .map(|r| (r, result.image[r * pixels_wide + col]))
                .filter(|(_, v)| !v.is_nan())
                .collect();
            if col_vals.is_empty() { continue; }

            let bright_near_mid = col_vals.iter().any(|(r, v)| {
                *v > 200.0 && (*r as isize - mid_row as isize).abs() <= 8
            });
            let dark_far = col_vals.iter().any(|(r, v)| {
                (*r as isize - mid_row as isize).abs() > 20 && *v < 200.0
            });

            n_col_checked += 1;
            if bright_near_mid && dark_far { n_col_passed += 1; }
        }

        assert!(n_col_checked > 0, "no valid columns in L-shape test");
        assert!(
            n_col_passed as f64 / n_col_checked as f64 >= 0.5,
            "L-shape: only {n_col_passed}/{n_col_checked} columns had vessel near mid-Y"
        );
    }
}
