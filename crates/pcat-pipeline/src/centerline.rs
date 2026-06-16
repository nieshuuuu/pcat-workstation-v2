use ndarray::Array3;

use crate::interp::trilinear;

/// Clip centerline to proximal segment [start_mm, start_mm + length_mm].
/// Input: centerline as voxel coords [z,y,x], spacing [sz,sy,sx].
/// Output: clipped subset of points.
pub fn clip_by_arclength(
    centerline: &[[f64; 3]],
    spacing: [f64; 3],
    start_mm: f64,
    length_mm: f64,
) -> Vec<[f64; 3]> {
    if centerline.len() < 2 {
        return centerline.to_vec();
    }

    let end_mm = start_mm + length_mm;

    // Compute cumulative arc-length in mm from point 0
    let mut cum_arc = Vec::with_capacity(centerline.len());
    cum_arc.push(0.0);
    for i in 1..centerline.len() {
        let dz = (centerline[i][0] - centerline[i - 1][0]) * spacing[0];
        let dy = (centerline[i][1] - centerline[i - 1][1]) * spacing[1];
        let dx = (centerline[i][2] - centerline[i - 1][2]) * spacing[2];
        let d = (dz * dz + dy * dy + dx * dx).sqrt();
        cum_arc.push(cum_arc[i - 1] + d);
    }

    // Retain points within [start_mm, end_mm].
    // Also interpolate boundary points for exact clipping.
    let mut result = Vec::new();

    for i in 0..centerline.len() {
        let s = cum_arc[i];

        // Check if we need to interpolate a start boundary point
        if i > 0 && cum_arc[i - 1] < start_mm && s >= start_mm {
            let seg_len = s - cum_arc[i - 1];
            if seg_len > 1e-12 {
                let t = (start_mm - cum_arc[i - 1]) / seg_len;
                result.push([
                    centerline[i - 1][0] + t * (centerline[i][0] - centerline[i - 1][0]),
                    centerline[i - 1][1] + t * (centerline[i][1] - centerline[i - 1][1]),
                    centerline[i - 1][2] + t * (centerline[i][2] - centerline[i - 1][2]),
                ]);
            }
        }

        // Include points within range
        if s >= start_mm && s <= end_mm {
            result.push(centerline[i]);
        }

        // Check if we need to interpolate an end boundary point
        if i > 0 && cum_arc[i - 1] <= end_mm && s > end_mm {
            let seg_len = s - cum_arc[i - 1];
            if seg_len > 1e-12 {
                let t = (end_mm - cum_arc[i - 1]) / seg_len;
                result.push([
                    centerline[i - 1][0] + t * (centerline[i][0] - centerline[i - 1][0]),
                    centerline[i - 1][1] + t * (centerline[i][1] - centerline[i - 1][1]),
                    centerline[i - 1][2] + t * (centerline[i][2] - centerline[i - 1][2]),
                ]);
            }
            break; // Past the end, no need to continue
        }
    }

    result
}

/// Estimate vessel radius at each centerline point.
/// For each point, samples radially in 8 directions on the cross-section plane
/// to find the distance to the first non-lumen voxel.
/// Clamps result to [0.5, 8.0] mm.
pub fn estimate_radii(
    volume: &Array3<f32>,
    centerline: &[[f64; 3]], // voxel coords [z,y,x]
    spacing: [f64; 3],
    lumen_range: (f32, f32), // default (150.0, 1200.0)
) -> Vec<f32> {
    use nalgebra::Vector3;

    if centerline.len() < 2 {
        return vec![2.0; centerline.len()]; // fallback
    }

    let n = centerline.len();
    let mut radii = Vec::with_capacity(n);

    // Compute tangent vectors for each point
    let tangents: Vec<Vector3<f64>> = (0..n)
        .map(|i| {
            let (prev, next) = if i == 0 {
                (0, 1)
            } else if i == n - 1 {
                (n - 2, n - 1)
            } else {
                (i - 1, i + 1)
            };
            let t = Vector3::new(
                centerline[next][0] - centerline[prev][0],
                centerline[next][1] - centerline[prev][1],
                centerline[next][2] - centerline[prev][2],
            );
            let norm = t.norm();
            if norm > 1e-12 {
                t / norm
            } else {
                Vector3::new(0.0, 0.0, 1.0)
            }
        })
        .collect();

    let step_mm = 0.1; // radial step size in mm
    let max_r_mm = 8.0;
    let n_directions = 8;

    for i in 0..n {
        let t = tangents[i];

        // Build a local cross-section frame (N, B) perpendicular to tangent
        let world_y = Vector3::new(0.0, 1.0, 0.0);
        let world_x = Vector3::new(1.0, 0.0, 0.0);
        let seed = if t.cross(&world_y).norm() > 0.1 {
            world_y
        } else {
            world_x
        };
        let n_vec = t.cross(&seed).normalize();
        let b_vec = t.cross(&n_vec).normalize();

        let center = Vector3::new(centerline[i][0], centerline[i][1], centerline[i][2]);

        let mut radius_sum = 0.0;
        let mut valid_count = 0;

        for dir_idx in 0..n_directions {
            let angle = 2.0 * std::f64::consts::PI * (dir_idx as f64) / (n_directions as f64);
            let direction = angle.cos() * n_vec + angle.sin() * b_vec;

            let mut r = 0.0;
            let mut found = false;
            while r <= max_r_mm {
                // Convert mm offset to voxel offset
                let sample = center
                    + Vector3::new(
                        direction[0] * r / spacing[0],
                        direction[1] * r / spacing[1],
                        direction[2] * r / spacing[2],
                    );

                let hu = trilinear(volume, sample[0], sample[1], sample[2]);
                if hu.is_nan() || hu < lumen_range.0 || hu > lumen_range.1 {
                    found = true;
                    break;
                }
                r += step_mm;
            }

            let boundary_r = if found { r } else { max_r_mm };
            radius_sum += boundary_r;
            valid_count += 1;
        }

        let avg_radius = if valid_count > 0 {
            radius_sum / valid_count as f64
        } else {
            2.0
        };

        // Clamp to [0.5, 8.0] mm
        radii.push(avg_radius.clamp(0.5, 8.0) as f32);
    }

    radii
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_clip_by_arclength_basic() {
        // Straight centerline along z-axis in voxel coords, spacing = 1mm
        let centerline: Vec<[f64; 3]> = (0..20)
            .map(|i| [i as f64, 0.0, 0.0])
            .collect();
        let spacing = [1.0, 1.0, 1.0];

        let clipped = clip_by_arclength(&centerline, spacing, 5.0, 10.0);
        // Should contain points from arc-length 5 to 15
        assert!(!clipped.is_empty());
        // First point should be at z=5
        assert!((clipped[0][0] - 5.0).abs() < 0.1);
        // Last point should be at z=15
        assert!((clipped.last().unwrap()[0] - 15.0).abs() < 0.1);
    }

    #[test]
    fn test_clip_by_arclength_with_spacing() {
        // Spacing of 2mm per voxel along z
        let centerline: Vec<[f64; 3]> = (0..10)
            .map(|i| [i as f64, 0.0, 0.0])
            .collect();
        let spacing = [2.0, 1.0, 1.0];

        // Total arc = 9 * 2 = 18mm. Clip [0, 10]
        let clipped = clip_by_arclength(&centerline, spacing, 0.0, 10.0);
        assert!(!clipped.is_empty());
        // At 10mm with 2mm spacing, we should reach z=5
        assert!((clipped.last().unwrap()[0] - 5.0).abs() < 0.1);
    }

    #[test]
    fn test_clip_empty_centerline() {
        let clipped = clip_by_arclength(&[], [1.0, 1.0, 1.0], 0.0, 10.0);
        assert!(clipped.is_empty());
    }

    #[test]
    fn test_estimate_radii_uniform_cylinder() {
        // Create a 32^3 volume with a cylinder of HU=300 along z at center (16,16)
        let mut vol = Array3::<f32>::zeros((32, 32, 32));
        let radius_vox = 3.0_f64;
        for z in 0..32 {
            for y in 0..32 {
                for x in 0..32 {
                    let dy = y as f64 - 16.0;
                    let dx = x as f64 - 16.0;
                    if (dy * dy + dx * dx).sqrt() <= radius_vox {
                        vol[[z, y, x]] = 300.0;
                    }
                }
            }
        }

        // Centerline along z at (y=16, x=16)
        let centerline: Vec<[f64; 3]> = (5..27)
            .map(|z| [z as f64, 16.0, 16.0])
            .collect();
        let spacing = [1.0, 1.0, 1.0];

        let radii = estimate_radii(&vol, &centerline, spacing, (150.0, 1200.0));
        assert_eq!(radii.len(), centerline.len());

        // All radii should be approximately 3mm (the cylinder radius)
        for &r in &radii {
            assert!(
                (r - 3.0).abs() < 1.0,
                "radius should be ~3.0, got {}",
                r
            );
        }
    }
}
