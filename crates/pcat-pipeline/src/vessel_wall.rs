//! Vessel wall / lumen geometry on a 2-D cross-section image.
//!
//! Per-ray FWHM boundary detection, anchored on a **globally estimated**
//! lumen HU. The threshold is shared across all rays and derived from a
//! small disk at the image centre (the centreline is, by construction,
//! inside the lumen). This makes the threshold immune to bright adjacent
//! structures — calcified plaque, neighbouring vessels, chambers — which
//! would otherwise corrupt a per-ray peak and collapse the polygon on the
//! affected side.

/// Per-angle vessel radii plus the derived equivalent-circle diameter and
/// boundary polygon.
pub struct VesselGeometry {
    /// Equivalent-circle diameter in mm, `D = 2·√(A/π)` where A is the polar
    /// integral of the per-angle radii.
    pub diameter_mm: f64,
    /// Boundary polygon in image pixel coords `[x, y]`, one point per angle,
    /// starting at θ = 0 (east) and walking counter-clockwise in image space
    /// (y-axis inverted, i.e. top-left origin).
    pub wall: Vec<[f64; 2]>,
    /// Raw radius at each angle in mm (same ordering as `wall`).
    pub r_theta: Vec<f64>,
}

/// Radius of the central disk used to estimate lumen HU, in mm. Picked to
/// sit comfortably inside even a narrow coronary artery (~2 mm radius) so
/// the median is dominated by lumen pixels even when the centreline is
/// placed a little off-axis.
const LUMEN_DISK_RADIUS_MM: f64 = 1.0;

/// Assumed soft-tissue / fat background HU surrounding a coronary artery.
/// Over-estimates mildly when the vessel abuts myocardium (~50 HU); that's
/// the accepted FWHM clinical convention for coronary sizing.
const BACKGROUND_HU: f64 = -100.0;

/// Minimum lumen-peak HU for the measurement to be trusted. Below this, the
/// centre disk is tissue (centreline drift, or no contrast at this slice),
/// and we return a degenerate polygon rather than guessing at a wall.
const MIN_LUMEN_HU: f64 = 50.0;

/// Upper bound on plausible coronary lumen radius. A ray that doesn't find a
/// sub-threshold crossing within this radius has almost certainly walked out
/// of the lumen into an adjacent contrast-filled structure (aortic root at
/// the ostium, cardiac chamber, another vessel). Capping the search here
/// keeps one rogue ray from pulling the polygon's equivalent-circle area up
/// by 10× — the polygon is slightly over-stated on that side, never wildly.
/// 3.5 mm → 7 mm diameter; normal RCA proximal lumen is ~4 mm, clinical max
/// ~6 mm, so the cap sits above real vessels but well below the aorta / any
/// chamber.
const MAX_CORONARY_RADIUS_MM: f64 = 3.5;

/// Compute the lumen boundary + diameter from a cross-section image.
///
/// * `image` — row-major f32 HU, `pixels × pixels`.
/// * `pixels` — image side length.
/// * `width_mm` — physical half-width of the image (image spans
///   `2·width_mm` on each side; matches the `width_mm` used by
///   `CprFrame::render_cross_section`).
/// * `n_angles` — number of radial rays (polygon vertex count).
pub fn compute_vessel_geometry(
    image: &[f32],
    pixels: usize,
    width_mm: f64,
    n_angles: usize,
) -> VesselGeometry {
    let mm_per_pixel = 2.0 * width_mm / pixels as f64;
    let center_px = pixels as f64 / 2.0;

    // 1. Estimate lumen HU from a small disk at image centre (median is
    //    robust to a few stray partial-volume pixels at the disk edge).
    let lumen_hu = estimate_lumen_hu(image, pixels, center_px, mm_per_pixel);
    let lumen_detectable = lumen_hu > MIN_LUMEN_HU;

    let threshold = (lumen_hu + BACKGROUND_HU) / 2.0;
    let step_px = 0.5_f64;
    // Bound the radial scan to a clinically-plausible coronary radius so a
    // ray can't walk through the lumen into an adjacent contrast-filled
    // structure (aorta, chamber). Always at least 0.5 mm so the loop runs.
    let max_radius_px =
        (MAX_CORONARY_RADIUS_MM / mm_per_pixel).min(center_px).max(2.0);
    let n_radial = (max_radius_px / step_px) as usize;

    // 2. Cast a ray per angle, detect boundary at the shared threshold.
    let mut r_theta: Vec<f64> = Vec::with_capacity(n_angles);
    for ai in 0..n_angles {
        if !lumen_detectable {
            r_theta.push(0.0);
            continue;
        }
        let theta = 2.0 * std::f64::consts::PI * (ai as f64) / (n_angles as f64);
        let r_px = detect_ray_boundary(
            image, pixels, center_px, theta, threshold, step_px, n_radial, max_radius_px,
        );
        r_theta.push(r_px * mm_per_pixel);
    }

    // 3. Cartesian polygon (pixel coords, y inverted for top-left origin).
    let wall: Vec<[f64; 2]> = (0..n_angles)
        .map(|ai| {
            let theta = 2.0 * std::f64::consts::PI * (ai as f64) / (n_angles as f64);
            let r_px = r_theta[ai] / mm_per_pixel;
            let x = center_px + r_px * theta.cos();
            let y = center_px - r_px * theta.sin();
            [x, y]
        })
        .collect();

    // 4. Equivalent-circle diameter from polar-integrated area.
    let dtheta = 2.0 * std::f64::consts::PI / (n_angles as f64);
    let area_mm2: f64 = 0.5 * r_theta.iter().map(|&r| r * r).sum::<f64>() * dtheta;
    let diameter_mm = 2.0 * (area_mm2 / std::f64::consts::PI).sqrt();

    VesselGeometry {
        diameter_mm,
        wall,
        r_theta,
    }
}

/// Median HU inside a small disk around the image centre — used as a
/// robust, calcium-insensitive estimate of the lumen HU. Returns
/// `f64::NEG_INFINITY` if no finite samples fall inside the disk.
fn estimate_lumen_hu(image: &[f32], pixels: usize, center_px: f64, mm_per_pixel: f64) -> f64 {
    let disk_r_px = (LUMEN_DISK_RADIUS_MM / mm_per_pixel).max(1.5);
    let disk_r_int = disk_r_px.ceil() as isize;
    let disk_r_sq = disk_r_px * disk_r_px;

    let mut samples: Vec<f64> = Vec::new();
    for dy in -disk_r_int..=disk_r_int {
        for dx in -disk_r_int..=disk_r_int {
            let d2 = (dx * dx + dy * dy) as f64;
            if d2 > disk_r_sq {
                continue;
            }
            let x = center_px + dx as f64;
            let y = center_px + dy as f64;
            let hu = sample_image_bilinear(image, pixels, x, y);
            if hu.is_finite() {
                samples.push(hu as f64);
            }
        }
    }

    if samples.is_empty() {
        return f64::NEG_INFINITY;
    }
    samples.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    samples[samples.len() / 2]
}

/// Cast one ray outward from the centre and return the boundary radius in
/// pixels. "Boundary" = first sample whose HU drops below `threshold`,
/// refined to the steepest negative gradient in a 7-sample window around
/// that crossing.
fn detect_ray_boundary(
    image: &[f32],
    pixels: usize,
    center_px: f64,
    theta: f64,
    threshold: f64,
    step_px: f64,
    n_radial: usize,
    max_radius_px: f64,
) -> f64 {
    let cos_t = theta.cos();
    let sin_t = theta.sin();

    let mut profile: Vec<f64> = Vec::with_capacity(n_radial);
    for ri in 0..n_radial {
        let r_px = (ri as f64) * step_px;
        let x = center_px + r_px * cos_t;
        let y = center_px - r_px * sin_t; // y inverted
        let hu = sample_image_bilinear(image, pixels, x, y);
        profile.push(if hu.is_nan() { -1024.0 } else { hu as f64 });
    }

    // First crossing below threshold, walking outward.
    let mut crossing_idx: Option<usize> = None;
    for i in 1..profile.len() {
        if profile[i] < threshold {
            crossing_idx = Some(i);
            break;
        }
    }
    let crossing = match crossing_idx {
        Some(idx) => idx,
        None => return max_radius_px,
    };

    // Refine to steepest negative gradient near the crossing.
    let window_start = crossing.saturating_sub(3);
    let window_end = (crossing + 3).min(profile.len().saturating_sub(1));
    let mut steepest_idx = crossing;
    let mut steepest_grad = 0.0_f64;
    for i in window_start..window_end {
        let grad = profile[i + 1] - profile[i];
        if grad < steepest_grad {
            steepest_grad = grad;
            steepest_idx = i;
        }
    }

    let r_px = (steepest_idx as f64 + 0.5) * step_px;
    r_px.clamp(0.5, max_radius_px)
}

/// Bilinear interpolation on a row-major f32 image.
fn sample_image_bilinear(image: &[f32], pixels: usize, x: f64, y: f64) -> f32 {
    let max_coord = (pixels as f64) - 1.0;
    if x < 0.0 || y < 0.0 || x > max_coord || y > max_coord {
        return f32::NAN;
    }

    let x = x.clamp(0.0, max_coord);
    let y = y.clamp(0.0, max_coord);

    let x0 = x.floor() as usize;
    let y0 = y.floor() as usize;
    let x1 = (x0 + 1).min(pixels - 1);
    let y1 = (y0 + 1).min(pixels - 1);

    let fx = (x - x0 as f64) as f32;
    let fy = (y - y0 as f64) as f32;

    let v00 = image[y0 * pixels + x0];
    let v10 = image[y0 * pixels + x1];
    let v01 = image[y1 * pixels + x0];
    let v11 = image[y1 * pixels + x1];

    v00 * (1.0 - fx) * (1.0 - fy)
        + v10 * fx * (1.0 - fy)
        + v01 * (1.0 - fx) * fy
        + v11 * fx * fy
}
