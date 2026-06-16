use nalgebra::Vector3;
use ndarray::Array3;

use crate::cpr::CprFrame;
use crate::interp::trilinear;

/// Result of contour extraction along a vessel centerline.
#[derive(Clone)]
pub struct ContourResult {
    /// Boundary radius at each position and angle, [n_positions][n_angles] in mm.
    pub r_theta: Vec<Vec<f64>>,
    /// Equivalent radius sqrt(area / pi) at each position, in mm.
    pub r_eq: Vec<f64>,
    /// Centerline positions in mm, [z, y, x].
    pub positions_mm: Vec<[f64; 3]>,
    /// Normal vectors (Bishop frame N) at each position.
    pub n_frame: Vec<[f64; 3]>,
    /// Binormal vectors (Bishop frame B) at each position.
    pub b_frame: Vec<[f64; 3]>,
    /// Arc-lengths along centerline in mm.
    #[allow(dead_code)]
    pub arclengths: Vec<f64>,
}

/// Extract vessel boundary contours using polar transform + gradient detection.
///
/// For each centerline position, samples radially at `n_angles` directions
/// from the center outward, detects the lumen boundary via half-max threshold
/// and steepest negative gradient, then smooths the resulting polar contour.
pub fn extract_contours(
    volume: &Array3<f32>,
    centerline: &[[f64; 3]], // voxel coords [z,y,x]
    spacing: [f64; 3],       // [sz,sy,sx]
    n_angles: usize,         // 360
    max_radius_mm: f64,      // 8.0
    sigma_deg: f64,          // 5.0 (smoothing)
) -> ContourResult {
    let n_pts = centerline.len();

    // 1. Convert centerline to mm
    let centerline_mm: Vec<[f64; 3]> = centerline
        .iter()
        .map(|pt| [pt[0] * spacing[0], pt[1] * spacing[1], pt[2] * spacing[2]])
        .collect();

    // 2. Build CprFrame — gives us positions, tangents, normals, binormals via
    //    cubic spline + Bishop frame (same as CPR rendering)
    let frame = CprFrame::from_centerline(&centerline_mm, n_pts);
    let positions: Vec<Vector3<f64>> = frame.positions.iter()
        .map(|p| Vector3::new(p[0], p[1], p[2]))
        .collect();
    let normals = frame.normals;
    let binormals = frame.binormals;
    let arclengths = frame.arclengths;

    // Radial sampling parameters
    let radial_step_mm = 0.2;
    let n_radial = ((max_radius_mm / radial_step_mm) as usize).max(1);

    // Inverse spacing for mm -> voxel
    let inv_spacing = [1.0 / spacing[0], 1.0 / spacing[1], 1.0 / spacing[2]];

    // Sigma in samples (angle index units)
    let sigma_samples = sigma_deg / (360.0 / n_angles as f64);

    let mut r_theta_all = Vec::with_capacity(n_pts);
    let mut r_eq_all = Vec::with_capacity(n_pts);

    for pos_idx in 0..n_pts {
        let pos = positions[pos_idx];
        let n_vec = normals[pos_idx];
        let b_vec = binormals[pos_idx];

        let mut r_theta = vec![max_radius_mm; n_angles];

        for angle_idx in 0..n_angles {
            let theta =
                2.0 * std::f64::consts::PI * (angle_idx as f64) / (n_angles as f64);
            let direction = theta.cos() * n_vec + theta.sin() * b_vec;

            // Sample radial HU profile
            let mut profile = Vec::with_capacity(n_radial);
            for ri in 0..n_radial {
                let r_mm = (ri as f64) * radial_step_mm;
                let sample_mm = pos + r_mm * direction;

                // Convert mm -> voxel
                let vz = sample_mm[0] * inv_spacing[0];
                let vy = sample_mm[1] * inv_spacing[1];
                let vx = sample_mm[2] * inv_spacing[2];

                let hu = trilinear(volume, vz, vy, vx);
                profile.push(if hu.is_nan() { -1024.0 } else { hu as f64 });
            }

            // Detect boundary using half-max threshold + steepest gradient
            r_theta[angle_idx] = detect_boundary(&profile, radial_step_mm, max_radius_mm);
        }

        // Smooth r_theta with circular Gaussian
        smooth_circular(&mut r_theta, sigma_samples);

        // Compute area via polar integration: A = 0.5 * sum(r^2 * dtheta)
        let dtheta = 2.0 * std::f64::consts::PI / (n_angles as f64);
        let area: f64 = 0.5 * r_theta.iter().map(|&r| r * r).sum::<f64>() * dtheta;
        let r_eq = (area / std::f64::consts::PI).sqrt();

        r_theta_all.push(r_theta);
        r_eq_all.push(r_eq);
    }

    // Convert nalgebra vectors to plain arrays
    let positions_mm: Vec<[f64; 3]> = positions.iter().map(|v| [v[0], v[1], v[2]]).collect();
    let n_frame: Vec<[f64; 3]> = normals.iter().map(|v| [v[0], v[1], v[2]]).collect();
    let b_frame: Vec<[f64; 3]> = binormals.iter().map(|v| [v[0], v[1], v[2]]).collect();

    ContourResult {
        r_theta: r_theta_all,
        r_eq: r_eq_all,
        positions_mm,
        n_frame,
        b_frame,
        arclengths,
    }
}

/// Detect the lumen boundary from a radial HU profile.
///
/// Strategy:
/// 1. Find maximum HU in the profile (lumen peak).
/// 2. Compute half-max threshold = max_hu / 2.
/// 3. Walk outward: find where HU drops below half-max.
/// 4. Refine: find the steepest negative gradient near that crossing.
/// 5. If no clear boundary found, return max_radius_mm.
fn detect_boundary(profile: &[f64], step_mm: f64, max_radius_mm: f64) -> f64 {
    if profile.is_empty() {
        return max_radius_mm;
    }

    // Find max HU in the first half of the profile (lumen should be near center)
    let search_len = (profile.len() / 2).max(3).min(profile.len());
    let max_hu = profile[..search_len]
        .iter()
        .cloned()
        .fold(f64::NEG_INFINITY, f64::max);

    if max_hu < 50.0 {
        // No real lumen signal, return max
        return max_radius_mm;
    }

    // Half-max threshold (background ~ -100 for soft tissue)
    let background = -100.0;
    let threshold = (max_hu + background) / 2.0;

    // Walk outward from center to find half-max crossing
    let mut crossing_idx = None;
    for i in 1..profile.len() {
        if profile[i] < threshold {
            crossing_idx = Some(i);
            break;
        }
    }

    let crossing = match crossing_idx {
        Some(idx) => idx,
        None => return max_radius_mm,
    };

    // Search for steepest negative gradient in a window around the crossing
    let window_start = if crossing >= 3 { crossing - 3 } else { 0 };
    let window_end = (crossing + 3).min(profile.len() - 1);

    let mut steepest_idx = crossing;
    let mut steepest_grad = 0.0_f64;

    for i in window_start..window_end {
        let grad = profile[i + 1] - profile[i]; // negative means dropping
        if grad < steepest_grad {
            steepest_grad = grad;
            steepest_idx = i;
        }
    }

    // Boundary is at the steepest gradient point, interpolated to sub-step precision
    let r = (steepest_idx as f64 + 0.5) * step_mm;
    r.clamp(0.5, max_radius_mm)
}

/// Apply 1D Gaussian smoothing with circular wrap-around.
fn smooth_circular(values: &mut [f64], sigma_samples: f64) {
    let n = values.len();
    if n == 0 || sigma_samples < 0.1 {
        return;
    }

    // Build Gaussian kernel
    let kernel_radius = (3.0 * sigma_samples).ceil() as usize;
    let kernel_size = 2 * kernel_radius + 1;
    let mut kernel = Vec::with_capacity(kernel_size);
    let mut kernel_sum = 0.0;
    for k in 0..kernel_size {
        let x = k as f64 - kernel_radius as f64;
        let w = (-0.5 * (x / sigma_samples).powi(2)).exp();
        kernel.push(w);
        kernel_sum += w;
    }
    // Normalize
    for w in &mut kernel {
        *w /= kernel_sum;
    }

    // Apply convolution with wrap-around
    let input = values.to_vec();
    for i in 0..n {
        let mut sum = 0.0;
        for (k, &w) in kernel.iter().enumerate() {
            let j = (i as isize + k as isize - kernel_radius as isize)
                .rem_euclid(n as isize) as usize;
            sum += w * input[j];
        }
        values[i] = sum;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_detect_boundary_clear_edge() {
        // Simulate a radial profile: high HU inside, drops to low outside
        let mut profile = vec![0.0; 40]; // 40 steps at 0.2mm = 8mm
        for i in 0..15 {
            profile[i] = 300.0; // lumen
        }
        for i in 15..20 {
            // gradient region
            profile[i] = 300.0 - (i - 15) as f64 * 80.0;
        }
        // rest stays at 0 (background)

        let r = detect_boundary(&profile, 0.2, 8.0);
        // Boundary should be around 15 * 0.2 = 3.0mm
        assert!(r > 1.0 && r < 5.0, "boundary at r={r}, expected ~3.0");
    }

    #[test]
    fn test_detect_boundary_no_signal() {
        let profile = vec![-100.0; 40];
        let r = detect_boundary(&profile, 0.2, 8.0);
        assert!((r - 8.0).abs() < 0.01, "should return max_radius when no lumen");
    }

    #[test]
    fn test_smooth_circular_identity() {
        // Uniform values should remain unchanged after smoothing
        let mut values = vec![5.0; 360];
        smooth_circular(&mut values, 5.0);
        for &v in &values {
            assert!((v - 5.0).abs() < 1e-10);
        }
    }

    #[test]
    fn test_smooth_circular_spike() {
        // A single spike should be spread out
        let mut values = vec![0.0; 360];
        values[180] = 360.0;
        let original_sum: f64 = values.iter().sum();

        smooth_circular(&mut values, 5.0);

        // Sum should be preserved
        let smoothed_sum: f64 = values.iter().sum();
        assert!(
            (smoothed_sum - original_sum).abs() < 1e-6,
            "sum changed: {} -> {}",
            original_sum,
            smoothed_sum
        );
        // Peak should be reduced
        assert!(values[180] < 360.0);
    }

    #[test]
    fn test_extract_contours_cylinder() {
        // Create a 32^3 volume with a cylinder of HU=300 along z at center (16,16)
        let mut vol = Array3::<f32>::zeros((32, 32, 32));
        let radius_vox = 4.0_f64;
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

        let centerline: Vec<[f64; 3]> = (5..27)
            .map(|z| [z as f64, 16.0, 16.0])
            .collect();
        let spacing = [1.0, 1.0, 1.0];

        let result = extract_contours(&vol, &centerline, spacing, 36, 8.0, 5.0);

        assert_eq!(result.positions_mm.len(), centerline.len());
        assert_eq!(result.r_eq.len(), centerline.len());
        assert_eq!(result.r_theta.len(), centerline.len());

        // Equivalent radius should be roughly the cylinder radius
        for &req in &result.r_eq[2..result.r_eq.len() - 2] {
            assert!(
                req > 1.0 && req < 8.0,
                "r_eq should be roughly cylinder radius, got {}",
                req
            );
        }
    }
}
