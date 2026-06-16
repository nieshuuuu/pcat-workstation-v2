use nalgebra::Vector3;
use ndarray::Array3;

use crate::contour::ContourResult;

/// VOI construction mode.
pub enum VoiMode {
    /// Fixed-width ring: gap_mm outside the boundary, then ring_mm thick shell.
    Crisp { gap_mm: f64, ring_mm: f64 },
    /// Scaled ring: boundary to factor * r_eq.
    #[allow(dead_code)]
    Scaled { factor: f64 },
}

/// Build a perivascular VOI (Volume of Interest) boolean mask.
///
/// For each centerline position, scans nearby voxels and checks whether they
/// fall within the specified shell around the vessel contour.
pub fn build_voi(
    volume_shape: [usize; 3], // [Z, Y, X]
    contours: &ContourResult,
    spacing: [f64; 3],
    mode: VoiMode,
) -> Array3<bool> {
    let [nz, ny, nx] = volume_shape;
    let mut mask = Array3::<bool>::default((nz, ny, nx));

    let n_positions = contours.positions_mm.len();
    if n_positions == 0 {
        return mask;
    }

    let n_angles = if !contours.r_theta.is_empty() {
        contours.r_theta[0].len()
    } else {
        return mask;
    };

    let inv_spacing = [1.0 / spacing[0], 1.0 / spacing[1], 1.0 / spacing[2]];

    // For each centerline position, determine a bounding box of voxels to check
    for pos_idx in 0..n_positions {
        let pos_mm = Vector3::new(
            contours.positions_mm[pos_idx][0],
            contours.positions_mm[pos_idx][1],
            contours.positions_mm[pos_idx][2],
        );
        let n_vec = Vector3::new(
            contours.n_frame[pos_idx][0],
            contours.n_frame[pos_idx][1],
            contours.n_frame[pos_idx][2],
        );
        let b_vec = Vector3::new(
            contours.b_frame[pos_idx][0],
            contours.b_frame[pos_idx][1],
            contours.b_frame[pos_idx][2],
        );

        let r_eq = contours.r_eq[pos_idx];

        // Determine the maximum shell outer radius for this position
        let max_outer_mm = match &mode {
            VoiMode::Crisp { gap_mm, ring_mm } => {
                // Find max boundary radius at this position
                let max_boundary = contours.r_theta[pos_idx]
                    .iter()
                    .cloned()
                    .fold(0.0_f64, f64::max);
                max_boundary + gap_mm + ring_mm + 1.0 // +1mm margin
            }
            VoiMode::Scaled { factor } => factor * r_eq + 1.0,
        };

        // Convert center to voxel coords
        let center_vox = [
            pos_mm[0] * inv_spacing[0],
            pos_mm[1] * inv_spacing[1],
            pos_mm[2] * inv_spacing[2],
        ];

        // Bounding box in voxel coords: scan a cube of max_outer_mm around center
        let margin_vox = [
            (max_outer_mm * inv_spacing[0]).ceil() as i64,
            (max_outer_mm * inv_spacing[1]).ceil() as i64,
            (max_outer_mm * inv_spacing[2]).ceil() as i64,
        ];

        let z_min = (center_vox[0] as i64 - margin_vox[0]).max(0) as usize;
        let z_max = ((center_vox[0] as i64 + margin_vox[0]) as usize).min(nz - 1);
        let y_min = (center_vox[1] as i64 - margin_vox[1]).max(0) as usize;
        let y_max = ((center_vox[1] as i64 + margin_vox[1]) as usize).min(ny - 1);
        let x_min = (center_vox[2] as i64 - margin_vox[2]).max(0) as usize;
        let x_max = ((center_vox[2] as i64 + margin_vox[2]) as usize).min(nx - 1);

        for vz in z_min..=z_max {
            for vy in y_min..=y_max {
                for vx in x_min..=x_max {
                    // Convert voxel to mm
                    let vox_mm = Vector3::new(
                        vz as f64 * spacing[0],
                        vy as f64 * spacing[1],
                        vx as f64 * spacing[2],
                    );

                    // Vector from centerline position to this voxel
                    let delta = vox_mm - pos_mm;

                    // Project onto the N-B plane (remove tangent component)
                    let proj_n = delta.dot(&n_vec);
                    let proj_b = delta.dot(&b_vec);

                    // Distance from centerline axis in mm
                    let r = (proj_n * proj_n + proj_b * proj_b).sqrt();

                    // Convert to polar angle
                    let theta = proj_b.atan2(proj_n);
                    let theta_pos = if theta < 0.0 {
                        theta + 2.0 * std::f64::consts::PI
                    } else {
                        theta
                    };

                    // Find the interpolated boundary radius at this angle
                    let angle_frac =
                        theta_pos / (2.0 * std::f64::consts::PI) * (n_angles as f64);
                    let angle_idx0 = (angle_frac as usize) % n_angles;
                    let angle_idx1 = (angle_idx0 + 1) % n_angles;
                    let t = angle_frac - angle_frac.floor();

                    let boundary_r = contours.r_theta[pos_idx][angle_idx0] * (1.0 - t)
                        + contours.r_theta[pos_idx][angle_idx1] * t;

                    // Check if this voxel is in the VOI shell
                    let in_shell = match &mode {
                        VoiMode::Crisp { gap_mm, ring_mm } => {
                            let inner = boundary_r + gap_mm;
                            let outer = inner + ring_mm;
                            r > inner && r <= outer
                        }
                        VoiMode::Scaled { factor } => {
                            r > boundary_r && r <= factor * r_eq
                        }
                    };

                    if in_shell {
                        mask[[vz, vy, vx]] = true;
                    }
                }
            }
        }
    }

    mask
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_cylinder_contours() -> ContourResult {
        // Simulate contours for a cylinder of radius 3mm along z, centered at (16,16,16) mm
        let n_positions = 10;
        let n_angles = 36;
        let radius = 3.0;

        let r_theta = vec![vec![radius; n_angles]; n_positions];
        let r_eq = vec![radius; n_positions];
        let positions_mm: Vec<[f64; 3]> = (0..n_positions)
            .map(|i| [(8 + i) as f64, 16.0, 16.0])
            .collect();
        // Frame: N = [0,1,0], B = [0,0,1] (perpendicular to z-axis tangent)
        let n_frame = vec![[0.0, 1.0, 0.0]; n_positions];
        let b_frame = vec![[0.0, 0.0, 1.0]; n_positions];
        let arclengths: Vec<f64> = (0..n_positions).map(|i| i as f64).collect();

        ContourResult {
            r_theta,
            r_eq,
            positions_mm,
            n_frame,
            b_frame,
            arclengths,
        }
    }

    #[test]
    fn test_build_voi_crisp() {
        let contours = make_cylinder_contours();
        let shape = [32, 32, 32];
        let spacing = [1.0, 1.0, 1.0];

        let mask = build_voi(
            shape,
            &contours,
            spacing,
            VoiMode::Crisp {
                gap_mm: 0.0,
                ring_mm: 3.0,
            },
        );

        // Count VOI voxels
        let count = mask.iter().filter(|&&v| v).count();
        assert!(count > 0, "VOI should have some voxels");

        // Center of the cylinder (y=16, x=16) at z=12 should NOT be in VOI
        // (it's inside the vessel, not in the shell)
        assert!(
            !mask[[12, 16, 16]],
            "center of vessel should not be in VOI"
        );
    }

    #[test]
    fn test_build_voi_scaled() {
        let contours = make_cylinder_contours();
        let shape = [32, 32, 32];
        let spacing = [1.0, 1.0, 1.0];

        let mask = build_voi(
            shape,
            &contours,
            spacing,
            VoiMode::Scaled { factor: 2.0 },
        );

        let count = mask.iter().filter(|&&v| v).count();
        assert!(count > 0, "scaled VOI should have some voxels");
    }

    #[test]
    fn test_build_voi_empty_contours() {
        let contours = ContourResult {
            r_theta: vec![],
            r_eq: vec![],
            positions_mm: vec![],
            n_frame: vec![],
            b_frame: vec![],
            arclengths: vec![],
        };
        let mask = build_voi([32, 32, 32], &contours, [1.0, 1.0, 1.0], VoiMode::Crisp { gap_mm: 0.0, ring_mm: 3.0 });
        let count = mask.iter().filter(|&&v| v).count();
        assert_eq!(count, 0);
    }
}
