use ndarray::Array3;
use serde::{Deserialize, Serialize};

use crate::contour::ContourResult;

/// Radial profile: mean FAI HU at concentric distances from the vessel wall.
#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct RadialProfile {
    /// Bin centers in mm (0.5, 1.5, 2.5, ..., 19.5).
    pub distances_mm: Vec<f64>,
    /// Mean FAI HU per ring (NaN if no voxels).
    pub mean_hu: Vec<f64>,
    /// Std FAI HU per ring.
    pub std_hu: Vec<f64>,
}

/// Per-sector statistics for angular asymmetry analysis.
#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct SectorStats {
    pub label: String,
    pub angle_deg: f64,
    pub hu_mean: f64,
    pub hu_std: f64,
    pub n_voxels: usize,
    pub fai_risk: String,
}

/// Angular asymmetry: mean FAI HU in angular sectors around the vessel.
#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct AngularAsymmetry {
    pub sectors: Vec<SectorStats>,
    /// Per-position heatmap: [n_positions][n_sectors].
    pub per_position_mean: Vec<Vec<f64>>,
}

const SECTOR_LABELS_8: [&str; 8] = [
    "Anterior", "Ant-Right", "Right", "Post-Right",
    "Posterior", "Post-Left", "Left", "Ant-Left",
];

const SECTOR_LABELS_16: [&str; 16] = [
    "N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE",
    "S", "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW",
];

/// FAI (Fat Attenuation Index) statistics for a single vessel.
#[derive(Serialize, Clone, Debug)]
pub struct FaiStats {
    /// Vessel name (e.g. "LAD", "LCx", "RCA").
    pub vessel: String,
    /// Total number of voxels in the perivascular VOI.
    pub n_voi_voxels: usize,
    /// Number of voxels with HU in the fat range within the VOI.
    pub n_fat_voxels: usize,
    /// Fraction of VOI voxels that are fat (n_fat / n_voi).
    pub fat_fraction: f64,
    /// Mean HU of fat voxels.
    pub hu_mean: f64,
    /// Standard deviation of HU of fat voxels.
    pub hu_std: f64,
    /// Median HU of fat voxels.
    pub hu_median: f64,
    /// Risk classification based on mean FAI.
    pub fai_risk: String,
    /// Histogram bin centers (100 bins from -200 to 200).
    pub histogram_bins: Vec<f64>,
    /// Histogram counts for each bin.
    pub histogram_counts: Vec<usize>,
    /// Radial profile analysis (filled by pipeline after initial stats).
    pub radial_profile: Option<RadialProfile>,
    /// Angular asymmetry analysis (filled by pipeline after initial stats).
    pub angular_asymmetry: Option<AngularAsymmetry>,
}

/// Compute PCAT (pericoronary adipose tissue) statistics.
///
/// Collects HU values from the VOI mask, filters to the fat HU range,
/// and computes summary statistics and a histogram.
pub fn compute_pcat_stats(
    volume: &Array3<f32>,
    voi_mask: &Array3<bool>,
    vessel: &str,
    hu_range: (f64, f64), // (-190.0, -30.0)
) -> FaiStats {
    let shape = volume.shape();
    let (nz, ny, nx) = (shape[0], shape[1], shape[2]);

    // 1. Collect all HU values where voi_mask is true
    let mut voi_values = Vec::new();
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                if voi_mask[[z, y, x]] {
                    voi_values.push(volume[[z, y, x]] as f64);
                }
            }
        }
    }

    let n_voi_voxels = voi_values.len();

    // 2. Filter to fat HU range
    let mut fat_values: Vec<f64> = voi_values
        .iter()
        .filter(|&&hu| hu >= hu_range.0 && hu <= hu_range.1)
        .copied()
        .collect();

    let n_fat_voxels = fat_values.len();

    // 3. Compute statistics
    let (hu_mean, hu_std, hu_median) = if fat_values.is_empty() {
        (0.0, 0.0, 0.0)
    } else {
        let mean = fat_values.iter().sum::<f64>() / fat_values.len() as f64;

        let variance = fat_values
            .iter()
            .map(|&v| (v - mean).powi(2))
            .sum::<f64>()
            / fat_values.len() as f64;
        let std = variance.sqrt();

        // Median
        fat_values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let median = if fat_values.len() % 2 == 0 {
            let mid = fat_values.len() / 2;
            (fat_values[mid - 1] + fat_values[mid]) / 2.0
        } else {
            fat_values[fat_values.len() / 2]
        };

        (mean, std, median)
    };

    // 4. Fat fraction
    let fat_fraction = if n_voi_voxels > 0 {
        n_fat_voxels as f64 / n_voi_voxels as f64
    } else {
        0.0
    };

    // 5. Risk classification
    // Based on Oikonomou et al. (2018): FAI > -70.1 HU indicates high risk
    let fai_risk = if hu_mean > -70.1 {
        "HIGH".to_string()
    } else {
        "LOW".to_string()
    };

    // 6. Histogram: 100 bins from -200 to 200
    let n_bins = 100;
    let bin_min = -200.0;
    let bin_max = 200.0;
    let bin_width = (bin_max - bin_min) / n_bins as f64;

    let histogram_bins: Vec<f64> = (0..n_bins)
        .map(|i| bin_min + (i as f64 + 0.5) * bin_width)
        .collect();
    let mut histogram_counts = vec![0usize; n_bins];

    for &hu in &voi_values {
        if hu >= bin_min && hu < bin_max {
            let idx = ((hu - bin_min) / bin_width) as usize;
            let idx = idx.min(n_bins - 1);
            histogram_counts[idx] += 1;
        }
    }

    FaiStats {
        vessel: vessel.to_string(),
        n_voi_voxels,
        n_fat_voxels,
        fat_fraction,
        hu_mean,
        hu_std,
        hu_median,
        fai_risk,
        histogram_bins,
        histogram_counts,
        radial_profile: None,
        angular_asymmetry: None,
    }
}

/// Unified FAI analysis — a SINGLE voxel pass that derives the FAI summary, the
/// radial profile, and the angular asymmetry from the SAME voxels and the SAME
/// contour geometry, so the three dashboard views reconcile by construction:
/// the count-weighted mean of the angular sectors equals the FAI mean exactly,
/// and the radial rings inside the VOI band reconcile with it too.
///
/// Each voxel is assigned to its nearest contour cross-section and projected
/// into that section's (N, B) frame → radial distance `r`, angle `theta`, and
/// `dist_from_wall = r − boundary_r(theta)` using the per-angle contour radius
/// (the identical geometry `voi::build_voi` uses). The VOI is the CRISP-CT
/// shell `dist_from_wall ∈ (gap_mm, gap_mm + ring_mm]`. The FAI is the mean of
/// fat-range voxels in that shell; the angular sectors partition exactly those
/// voxels; the radial profile bins all fat voxels by `dist_from_wall` out to
/// `max_distance_mm` (so it shows the full gradient, and its in-band rings are
/// the same voxels as the VOI).
///
/// Sector angles are anchored to a fixed world reference projected into each
/// section, so a given physical direction maps to the same sector at every
/// position — the section's own N/B frame twists along the vessel (Bishop
/// frame) and would otherwise smear one physical direction across sectors.
/// Sectors are labelled by their angle from that anchor (e.g. "0°", "45°"),
/// not anatomical names, because the anchor is a consistent reference rather
/// than validated patient anatomy.
#[allow(clippy::too_many_arguments)]
pub fn compute_fai_analysis(
    volume: &Array3<f32>,
    contours: &ContourResult,
    spacing: [f64; 3],
    vessel: &str,
    hu_range: (f64, f64),
    gap_mm: f64,
    ring_mm: f64,
    n_sectors: usize,
    max_distance_mm: f64,
    ring_step_mm: f64,
) -> FaiStats {
    use nalgebra::Vector3;
    use std::f64::consts::PI;
    let tau = 2.0 * PI;

    let shape = volume.shape();
    let (nz, ny, nx) = (shape[0], shape[1], shape[2]);
    let n_pos = contours.positions_mm.len();
    let n_angles = contours.r_theta.first().map_or(0, |v| v.len());

    let n_radial = (max_distance_mm / ring_step_mm).max(1.0) as usize;
    let distances_mm: Vec<f64> = (0..n_radial)
        .map(|i| (i as f64 + 0.5) * ring_step_mm)
        .collect();

    // Histogram axis (100 bins, -200..200) — same as compute_pcat_stats.
    let (n_bins, bmin, bmax) = (100usize, -200.0_f64, 200.0_f64);
    let bw = (bmax - bmin) / n_bins as f64;
    let histogram_bins: Vec<f64> = (0..n_bins)
        .map(|i| bmin + (i as f64 + 0.5) * bw)
        .collect();

    // Degenerate guard: no contour → empty but well-formed result.
    if n_pos == 0 || n_angles == 0 {
        return FaiStats {
            vessel: vessel.to_string(),
            n_voi_voxels: 0,
            n_fat_voxels: 0,
            fat_fraction: 0.0,
            hu_mean: 0.0,
            hu_std: 0.0,
            hu_median: 0.0,
            fai_risk: "N/A".to_string(),
            histogram_bins,
            histogram_counts: vec![0; n_bins],
            radial_profile: None,
            angular_asymmetry: None,
        };
    }

    // Section frames as vectors.
    let pos: Vec<Vector3<f64>> = contours
        .positions_mm
        .iter()
        .map(|p| Vector3::new(p[0], p[1], p[2]))
        .collect();
    let nf: Vec<Vector3<f64>> = contours
        .n_frame
        .iter()
        .map(|p| Vector3::new(p[0], p[1], p[2]))
        .collect();
    let bf: Vec<Vector3<f64>> = contours
        .b_frame
        .iter()
        .map(|p| Vector3::new(p[0], p[1], p[2]))
        .collect();

    // Per-section angle anchor: angle (in the N/B basis) of a fixed world
    // reference projected into the section. Subtracting it pins sector 0 to the
    // same world direction at every position. Reference = patient anterior-ish
    // (−Y in the [z,y,x] mm frame); fall back to +Z (slice axis) where the
    // reference is nearly perpendicular to the section.
    let theta_ref: Vec<f64> = (0..n_pos)
        .map(|i| {
            let primary = Vector3::new(0.0, -1.0, 0.0);
            let (mut rn, mut rb) = (primary.dot(&nf[i]), primary.dot(&bf[i]));
            if (rn * rn + rb * rb).sqrt() < 1e-6 {
                let fallback = Vector3::new(1.0, 0.0, 0.0);
                rn = fallback.dot(&nf[i]);
                rb = fallback.dot(&bf[i]);
            }
            rb.atan2(rn)
        })
        .collect();

    // Global bounding box around all sections + margin.
    let max_boundary = contours
        .r_theta
        .iter()
        .flat_map(|v| v.iter().copied())
        .fold(0.0_f64, f64::max);
    let pad = max_boundary + max_distance_mm + 2.0;
    let inv = [1.0 / spacing[0], 1.0 / spacing[1], 1.0 / spacing[2]];
    let (mut zmn, mut zmx) = (f64::MAX, f64::MIN);
    let (mut ymn, mut ymx) = (f64::MAX, f64::MIN);
    let (mut xmn, mut xmx) = (f64::MAX, f64::MIN);
    for p in &pos {
        zmn = zmn.min(p[0]);
        zmx = zmx.max(p[0]);
        ymn = ymn.min(p[1]);
        ymx = ymx.max(p[1]);
        xmn = xmn.min(p[2]);
        xmx = xmx.max(p[2]);
    }
    let iz0 = ((zmn - pad) * inv[0]).floor().max(0.0) as usize;
    let iz1 = (((zmx + pad) * inv[0]).ceil() as usize).min(nz - 1);
    let iy0 = ((ymn - pad) * inv[1]).floor().max(0.0) as usize;
    let iy1 = (((ymx + pad) * inv[1]).ceil() as usize).min(ny - 1);
    let ix0 = ((xmn - pad) * inv[2]).floor().max(0.0) as usize;
    let ix1 = (((xmx + pad) * inv[2]).ceil() as usize).min(nx - 1);

    // Subsample the nearest-section search for very long centerlines.
    let cl_step = if n_pos > 200 { n_pos / 200 } else { 1 };

    // Accumulators.
    let mut bin_sum = vec![0.0f64; n_radial];
    let mut bin_sq = vec![0.0f64; n_radial];
    let mut bin_cnt = vec![0usize; n_radial];
    let mut sec_sum = vec![0.0f64; n_sectors];
    let mut sec_sq = vec![0.0f64; n_sectors];
    let mut sec_cnt = vec![0usize; n_sectors];
    let mut pp_sum = vec![vec![0.0f64; n_sectors]; n_pos];
    let mut pp_cnt = vec![vec![0usize; n_sectors]; n_pos];
    let mut voi_values: Vec<f64> = Vec::new();
    let mut fat_values: Vec<f64> = Vec::new();

    for vz in iz0..=iz1 {
        for vy in iy0..=iy1 {
            for vx in ix0..=ix1 {
                let vm = Vector3::new(
                    vz as f64 * spacing[0],
                    vy as f64 * spacing[1],
                    vx as f64 * spacing[2],
                );
                // Nearest section.
                let (mut best, mut bestd) = (0usize, f64::MAX);
                let mut i = 0;
                while i < n_pos {
                    let d = (vm - pos[i]).norm_squared();
                    if d < bestd {
                        bestd = d;
                        best = i;
                    }
                    i += cl_step;
                }
                let p = best;
                let delta = vm - pos[p];
                let proj_n = delta.dot(&nf[p]);
                let proj_b = delta.dot(&bf[p]);
                let r = (proj_n * proj_n + proj_b * proj_b).sqrt();
                let theta = proj_b.atan2(proj_n);
                let theta_pos = if theta < 0.0 { theta + tau } else { theta };

                // Interpolated per-angle boundary radius (same as build_voi).
                let af = theta_pos / tau * n_angles as f64;
                let a0 = (af as usize) % n_angles;
                let a1 = (a0 + 1) % n_angles;
                let tt = af - af.floor();
                let boundary_r =
                    contours.r_theta[p][a0] * (1.0 - tt) + contours.r_theta[p][a1] * tt;
                let dist_from_wall = r - boundary_r;

                let hu = volume[[vz, vy, vx]] as f64;
                let is_fat = hu >= hu_range.0 && hu <= hu_range.1;

                // Radial profile: every fat voxel, distance from the outer wall.
                if is_fat && dist_from_wall >= 0.0 && dist_from_wall < max_distance_mm {
                    let b = ((dist_from_wall / ring_step_mm) as usize).min(n_radial - 1);
                    bin_sum[b] += hu;
                    bin_sq[b] += hu * hu;
                    bin_cnt[b] += 1;
                }

                // CRISP-CT VOI shell → FAI + angular.
                let in_voi = dist_from_wall > gap_mm && dist_from_wall <= gap_mm + ring_mm;
                if in_voi {
                    voi_values.push(hu);
                    if is_fat {
                        fat_values.push(hu);
                        let anchored = {
                            let a = (theta - theta_ref[p]) % tau;
                            if a < 0.0 { a + tau } else { a }
                        };
                        let s = ((anchored / tau * n_sectors as f64) as usize) % n_sectors;
                        sec_sum[s] += hu;
                        sec_sq[s] += hu * hu;
                        sec_cnt[s] += 1;
                        pp_sum[p][s] += hu;
                        pp_cnt[p][s] += 1;
                    }
                }
            }
        }
    }

    // FAI summary over (in_voi ∧ fat).
    let n_voi_voxels = voi_values.len();
    let n_fat_voxels = fat_values.len();
    fat_values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let (hu_mean, hu_std, hu_median) = if fat_values.is_empty() {
        (0.0, 0.0, 0.0)
    } else {
        let mean = fat_values.iter().sum::<f64>() / n_fat_voxels as f64;
        let var = fat_values.iter().map(|&v| (v - mean).powi(2)).sum::<f64>() / n_fat_voxels as f64;
        let median = if n_fat_voxels % 2 == 0 {
            (fat_values[n_fat_voxels / 2 - 1] + fat_values[n_fat_voxels / 2]) / 2.0
        } else {
            fat_values[n_fat_voxels / 2]
        };
        (mean, var.sqrt(), median)
    };
    let fat_fraction = if n_voi_voxels > 0 {
        n_fat_voxels as f64 / n_voi_voxels as f64
    } else {
        0.0
    };
    let fai_risk = if !fat_values.is_empty() && hu_mean > -70.1 {
        "HIGH"
    } else {
        "LOW"
    }
    .to_string();

    // Histogram over VOI voxels.
    let mut histogram_counts = vec![0usize; n_bins];
    for &hu in &voi_values {
        if hu >= bmin && hu < bmax {
            let idx = (((hu - bmin) / bw) as usize).min(n_bins - 1);
            histogram_counts[idx] += 1;
        }
    }

    // Radial profile.
    let mean_hu: Vec<f64> = (0..n_radial)
        .map(|i| {
            if bin_cnt[i] == 0 {
                f64::NAN
            } else {
                bin_sum[i] / bin_cnt[i] as f64
            }
        })
        .collect();
    let std_hu: Vec<f64> = (0..n_radial)
        .map(|i| {
            if bin_cnt[i] == 0 {
                f64::NAN
            } else {
                let m = bin_sum[i] / bin_cnt[i] as f64;
                (bin_sq[i] / bin_cnt[i] as f64 - m * m).max(0.0).sqrt()
            }
        })
        .collect();
    let radial_profile = RadialProfile {
        distances_mm,
        mean_hu,
        std_hu,
    };

    // Angular asymmetry — sectors labelled by angle from the anchor.
    let sectors: Vec<SectorStats> = (0..n_sectors)
        .map(|s| {
            let angle_deg = s as f64 * 360.0 / n_sectors as f64;
            let (hu_mean_s, hu_std_s) = if sec_cnt[s] == 0 {
                (f64::NAN, f64::NAN)
            } else {
                let m = sec_sum[s] / sec_cnt[s] as f64;
                (m, (sec_sq[s] / sec_cnt[s] as f64 - m * m).max(0.0).sqrt())
            };
            let fai_risk = if hu_mean_s.is_nan() || hu_mean_s <= -70.1 {
                "LOW".to_string()
            } else {
                "HIGH".to_string()
            };
            SectorStats {
                label: format!("{angle_deg:.0}°"),
                angle_deg,
                hu_mean: hu_mean_s,
                hu_std: hu_std_s,
                n_voxels: sec_cnt[s],
                fai_risk,
            }
        })
        .collect();
    let per_position_mean: Vec<Vec<f64>> = (0..n_pos)
        .map(|p| {
            (0..n_sectors)
                .map(|s| {
                    if pp_cnt[p][s] == 0 {
                        f64::NAN
                    } else {
                        pp_sum[p][s] / pp_cnt[p][s] as f64
                    }
                })
                .collect()
        })
        .collect();
    let angular_asymmetry = AngularAsymmetry {
        sectors,
        per_position_mean,
    };

    FaiStats {
        vessel: vessel.to_string(),
        n_voi_voxels,
        n_fat_voxels,
        fat_fraction,
        hu_mean,
        hu_std,
        hu_median,
        fai_risk,
        histogram_bins,
        histogram_counts,
        radial_profile: Some(radial_profile),
        angular_asymmetry: Some(angular_asymmetry),
    }
}

/// **Superseded by [`compute_fai_analysis`]** (kept as a standalone helper; the
/// live pipeline no longer calls it). Radial profile: mean FAI HU at concentric
/// distances from the vessel wall.
///
/// For each voxel near the centerline, distance from the wall =
/// distance_from_centerline − local_radius − `gap_mm`. The `gap_mm` offset skips
/// the vessel wall (intima/media/adventitia) so the profile measures
/// perivascular fat, not the wall itself — matching the FAI VOI's gap. With it,
/// distance 0 sits at the outer wall (where fat begins), so the near-wall bins
/// no longer pick up wall partial-volume that spuriously raises the HU.
/// Voxels are binned into `ring_step_mm` rings, filtered to the FAI HU range,
/// and reduced to mean/std per ring.
pub fn compute_radial_profile(
    volume: &Array3<f32>,
    centerline_vox: &[[f64; 3]],
    radii_mm: &[f32],
    spacing: [f64; 3],
    max_distance_mm: f64,
    ring_step_mm: f64,
    gap_mm: f64,
    hu_range: (f64, f64),
) -> RadialProfile {
    let n_bins = (max_distance_mm / ring_step_mm) as usize;
    let distances_mm: Vec<f64> = (0..n_bins)
        .map(|i| (i as f64 + 1.0) * ring_step_mm)
        .collect();

    // Accumulators per bin: sum, sum_sq, count
    let mut bin_sum = vec![0.0f64; n_bins];
    let mut bin_sum_sq = vec![0.0f64; n_bins];
    let mut bin_count = vec![0usize; n_bins];

    let shape = volume.shape();
    let (nz, ny, nx) = (shape[0], shape[1], shape[2]);

    // Bounding box around centerline with padding
    let max_radius: f64 = if radii_mm.is_empty() {
        2.0
    } else {
        radii_mm.iter().map(|&r| r as f64).fold(0.0f64, f64::max)
    };
    let pad = max_distance_mm + max_radius;
    let mut z_min = f64::MAX;
    let mut z_max = f64::MIN;
    let mut y_min = f64::MAX;
    let mut y_max = f64::MIN;
    let mut x_min = f64::MAX;
    let mut x_max = f64::MIN;

    for pt in centerline_vox {
        z_min = z_min.min(pt[0]);
        z_max = z_max.max(pt[0]);
        y_min = y_min.min(pt[1]);
        y_max = y_max.max(pt[1]);
        x_min = x_min.min(pt[2]);
        x_max = x_max.max(pt[2]);
    }

    // Convert pad from mm to voxel units per axis
    let pad_z = pad / spacing[0];
    let pad_y = pad / spacing[1];
    let pad_x = pad / spacing[2];

    let iz_lo = ((z_min - pad_z).floor().max(0.0)) as usize;
    let iz_hi = ((z_max + pad_z).ceil() as usize).min(nz - 1);
    let iy_lo = ((y_min - pad_y).floor().max(0.0)) as usize;
    let iy_hi = ((y_max + pad_y).ceil() as usize).min(ny - 1);
    let ix_lo = ((x_min - pad_x).floor().max(0.0)) as usize;
    let ix_hi = ((x_max + pad_x).ceil() as usize).min(nx - 1);

    // Subsample centerline for performance if very long (>200 points)
    let cl_step = if centerline_vox.len() > 200 { centerline_vox.len() / 200 } else { 1 };
    let cl_subsample: Vec<(usize, &[f64; 3])> = centerline_vox.iter()
        .enumerate()
        .step_by(cl_step)
        .collect();

    for z in iz_lo..=iz_hi {
        for y in iy_lo..=iy_hi {
            for x in ix_lo..=ix_hi {
                // Distance to nearest centerline point in mm
                let mut min_dist_sq = f64::MAX;
                let mut nearest_idx = 0usize;
                for &(ci, pt) in &cl_subsample {
                    let dz = (z as f64 - pt[0]) * spacing[0];
                    let dy = (y as f64 - pt[1]) * spacing[1];
                    let dx = (x as f64 - pt[2]) * spacing[2];
                    let dsq = dz * dz + dy * dy + dx * dx;
                    if dsq < min_dist_sq {
                        min_dist_sq = dsq;
                        nearest_idx = ci;
                    }
                }
                let dist_from_center = min_dist_sq.sqrt();
                let local_radius = if nearest_idx < radii_mm.len() {
                    radii_mm[nearest_idx] as f64
                } else if !radii_mm.is_empty() {
                    *radii_mm.last().unwrap() as f64
                } else {
                    2.0
                };
                let dist_from_wall = dist_from_center - local_radius - gap_mm;

                if dist_from_wall < 0.0 || dist_from_wall >= max_distance_mm {
                    continue;
                }

                let hu = volume[[z, y, x]] as f64;
                if hu < hu_range.0 || hu > hu_range.1 {
                    continue;
                }

                let bin = (dist_from_wall / ring_step_mm) as usize;
                let bin = bin.min(n_bins - 1);
                bin_sum[bin] += hu;
                bin_sum_sq[bin] += hu * hu;
                bin_count[bin] += 1;
            }
        }
    }

    let mean_hu: Vec<f64> = (0..n_bins)
        .map(|i| {
            if bin_count[i] == 0 {
                f64::NAN
            } else {
                bin_sum[i] / bin_count[i] as f64
            }
        })
        .collect();

    let std_hu: Vec<f64> = (0..n_bins)
        .map(|i| {
            if bin_count[i] == 0 {
                f64::NAN
            } else {
                let mean = bin_sum[i] / bin_count[i] as f64;
                let var = bin_sum_sq[i] / bin_count[i] as f64 - mean * mean;
                var.max(0.0).sqrt()
            }
        })
        .collect();

    RadialProfile {
        distances_mm,
        mean_hu,
        std_hu,
    }
}

/// **Superseded by [`compute_fai_analysis`]** (kept as a standalone helper; the
/// live pipeline no longer calls it). Angular asymmetry: mean FAI HU in angular
/// sectors around the vessel.
///
/// Divides the pericoronary ring into `n_sectors` angular sectors (octants by default),
/// samples cross-sections along the centerline, and computes per-sector statistics.
#[allow(dead_code)]
pub fn compute_angular_asymmetry(
    volume: &Array3<f32>,
    centerline_vox: &[[f64; 3]],
    radii_mm: &[f32],
    spacing: [f64; 3],
    n_sectors: usize,
    hu_range: (f64, f64),
    gap_mm: f64,
    ring_mm: f64,
) -> AngularAsymmetry {
    let n_pts = centerline_vox.len();
    let step = if n_pts > 60 { n_pts / 60 } else { 1 };
    let shape = volume.shape();
    let (nz, ny, nx) = (shape[0], shape[1], shape[2]);

    // Global accumulators per sector
    let mut sector_sum = vec![0.0f64; n_sectors];
    let mut sector_sum_sq = vec![0.0f64; n_sectors];
    let mut sector_count = vec![0usize; n_sectors];
    let mut per_position_mean: Vec<Vec<f64>> = Vec::new();

    let mut idx = 0;
    while idx < n_pts {
        let pt = centerline_vox[idx];
        let radius = if idx < radii_mm.len() {
            radii_mm[idx] as f64
        } else if !radii_mm.is_empty() {
            *radii_mm.last().unwrap() as f64
        } else {
            2.0
        };

        // Tangent via finite differences
        let tangent = if n_pts < 2 {
            [0.0, 0.0, 1.0]
        } else if idx == 0 {
            let next = centerline_vox[1];
            let dz = (next[0] - pt[0]) * spacing[0];
            let dy = (next[1] - pt[1]) * spacing[1];
            let dx = (next[2] - pt[2]) * spacing[2];
            let len = (dz * dz + dy * dy + dx * dx).sqrt().max(1e-12);
            [dz / len, dy / len, dx / len]
        } else if idx >= n_pts - 1 {
            let prev = centerline_vox[n_pts - 2];
            let dz = (pt[0] - prev[0]) * spacing[0];
            let dy = (pt[1] - prev[1]) * spacing[1];
            let dx = (pt[2] - prev[2]) * spacing[2];
            let len = (dz * dz + dy * dy + dx * dx).sqrt().max(1e-12);
            [dz / len, dy / len, dx / len]
        } else {
            let prev = centerline_vox[idx - 1];
            let next = centerline_vox[idx + 1];
            let dz = (next[0] - prev[0]) * spacing[0];
            let dy = (next[1] - prev[1]) * spacing[1];
            let dx = (next[2] - prev[2]) * spacing[2];
            let len = (dz * dz + dy * dy + dx * dx).sqrt().max(1e-12);
            [dz / len, dy / len, dx / len]
        };

        // Compute normal and binormal perpendicular to tangent
        // Pick a reference vector not parallel to tangent
        let ref_vec = if tangent[0].abs() < 0.9 {
            [1.0, 0.0, 0.0]
        } else {
            [0.0, 1.0, 0.0]
        };

        // normal = normalize(ref_vec - (ref_vec . tangent) * tangent)
        let dot = ref_vec[0] * tangent[0] + ref_vec[1] * tangent[1] + ref_vec[2] * tangent[2];
        let n0 = ref_vec[0] - dot * tangent[0];
        let n1 = ref_vec[1] - dot * tangent[1];
        let n2 = ref_vec[2] - dot * tangent[2];
        let nlen = (n0 * n0 + n1 * n1 + n2 * n2).sqrt().max(1e-12);
        let normal = [n0 / nlen, n1 / nlen, n2 / nlen];

        // binormal = tangent x normal
        let binormal = [
            tangent[1] * normal[2] - tangent[2] * normal[1],
            tangent[2] * normal[0] - tangent[0] * normal[2],
            tangent[0] * normal[1] - tangent[1] * normal[0],
        ];

        let mut pos_sector_sum = vec![0.0f64; n_sectors];
        let mut pos_sector_count = vec![0usize; n_sectors];

        for s in 0..n_sectors {
            let angle = (s as f64 + 0.5) * std::f64::consts::TAU / n_sectors as f64;
            let cos_a = angle.cos();
            let sin_a = angle.sin();

            // Sample 5 radial points
            for ri in 0..5 {
                let r = radius + gap_mm + ring_mm * (ri as f64 + 0.5) / 5.0;

                // Offset in mm along normal and binormal
                let off_z = r * (cos_a * normal[0] + sin_a * binormal[0]);
                let off_y = r * (cos_a * normal[1] + sin_a * binormal[1]);
                let off_x = r * (cos_a * normal[2] + sin_a * binormal[2]);

                // Convert mm offset to voxel offset
                let vz = pt[0] + off_z / spacing[0];
                let vy = pt[1] + off_y / spacing[1];
                let vx = pt[2] + off_x / spacing[2];

                // Nearest-neighbor sampling
                let iz = vz.round() as isize;
                let iy = vy.round() as isize;
                let ix = vx.round() as isize;

                if iz < 0 || iy < 0 || ix < 0 {
                    continue;
                }
                let (uz, uy, ux) = (iz as usize, iy as usize, ix as usize);
                if uz >= nz || uy >= ny || ux >= nx {
                    continue;
                }

                let hu = volume[[uz, uy, ux]] as f64;
                if hu < hu_range.0 || hu > hu_range.1 {
                    continue;
                }

                sector_sum[s] += hu;
                sector_sum_sq[s] += hu * hu;
                sector_count[s] += 1;

                pos_sector_sum[s] += hu;
                pos_sector_count[s] += 1;
            }
        }

        // Per-position means
        let pos_means: Vec<f64> = (0..n_sectors)
            .map(|s| {
                if pos_sector_count[s] == 0 {
                    f64::NAN
                } else {
                    pos_sector_sum[s] / pos_sector_count[s] as f64
                }
            })
            .collect();
        per_position_mean.push(pos_means);

        idx += step;
    }

    // Build sector stats
    let labels: Vec<&str> = if n_sectors == 16 {
        SECTOR_LABELS_16.to_vec()
    } else if n_sectors == 8 {
        SECTOR_LABELS_8.to_vec()
    } else {
        (0..n_sectors)
            .map(|i| match i {
                0 => "S0",
                1 => "S1",
                2 => "S2",
                3 => "S3",
                4 => "S4",
                5 => "S5",
                6 => "S6",
                7 => "S7",
                _ => "S?",
            })
            .collect()
    };

    let sectors: Vec<SectorStats> = (0..n_sectors)
        .map(|s| {
            let angle_deg = (s as f64 + 0.5) * 360.0 / n_sectors as f64;
            let (hu_mean, hu_std) = if sector_count[s] == 0 {
                (f64::NAN, f64::NAN)
            } else {
                let mean = sector_sum[s] / sector_count[s] as f64;
                let var = sector_sum_sq[s] / sector_count[s] as f64 - mean * mean;
                (mean, var.max(0.0).sqrt())
            };
            let fai_risk = if hu_mean.is_nan() || hu_mean <= -70.1 {
                "LOW".to_string()
            } else {
                "HIGH".to_string()
            };
            let label = if s < labels.len() {
                labels[s].to_string()
            } else {
                format!("S{s}")
            };
            SectorStats {
                label,
                angle_deg,
                hu_mean,
                hu_std,
                n_voxels: sector_count[s],
                fai_risk,
            }
        })
        .collect();

    AngularAsymmetry {
        sectors,
        per_position_mean,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::Array3;

    #[test]
    fn test_compute_pcat_stats_basic() {
        // Create a small volume and mask
        let mut vol = Array3::<f32>::zeros((10, 10, 10));
        let mut mask = Array3::<bool>::default((10, 10, 10));

        // Fill some voxels with fat-range HU values
        for z in 2..8 {
            for y in 2..8 {
                for x in 2..8 {
                    vol[[z, y, x]] = -80.0; // typical fat HU
                    mask[[z, y, x]] = true;
                }
            }
        }

        let stats = compute_pcat_stats(&vol, &mask, "LAD", (-190.0, -30.0));

        assert_eq!(stats.vessel, "LAD");
        assert_eq!(stats.n_voi_voxels, 6 * 6 * 6);
        assert_eq!(stats.n_fat_voxels, 6 * 6 * 6); // all are in fat range
        assert!((stats.fat_fraction - 1.0).abs() < 1e-10);
        assert!((stats.hu_mean - (-80.0)).abs() < 0.1);
        assert!(stats.hu_std < 0.01); // all same value
        assert!((stats.hu_median - (-80.0)).abs() < 0.1);
        assert_eq!(stats.fai_risk, "LOW"); // -80 < -70.1
    }

    #[test]
    fn test_compute_pcat_stats_high_risk() {
        let mut vol = Array3::<f32>::zeros((10, 10, 10));
        let mut mask = Array3::<bool>::default((10, 10, 10));

        // HU = -60 is above -70.1 threshold
        for z in 0..5 {
            for y in 0..5 {
                for x in 0..5 {
                    vol[[z, y, x]] = -60.0;
                    mask[[z, y, x]] = true;
                }
            }
        }

        let stats = compute_pcat_stats(&vol, &mask, "RCA", (-190.0, -30.0));
        assert_eq!(stats.fai_risk, "HIGH"); // -60 > -70.1
        assert!((stats.hu_mean - (-60.0)).abs() < 0.1);
    }

    #[test]
    fn test_compute_pcat_stats_empty_voi() {
        let vol = Array3::<f32>::zeros((10, 10, 10));
        let mask = Array3::<bool>::default((10, 10, 10));

        let stats = compute_pcat_stats(&vol, &mask, "LCx", (-190.0, -30.0));

        assert_eq!(stats.n_voi_voxels, 0);
        assert_eq!(stats.n_fat_voxels, 0);
        assert!((stats.fat_fraction - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_compute_pcat_stats_mixed() {
        let mut vol = Array3::<f32>::zeros((10, 10, 10));
        let mut mask = Array3::<bool>::default((10, 10, 10));

        // Mix of fat and non-fat in the VOI
        for x in 0..10 {
            vol[[5, 5, x]] = if x < 5 { -80.0 } else { 50.0 }; // fat vs non-fat
            mask[[5, 5, x]] = true;
        }

        let stats = compute_pcat_stats(&vol, &mask, "LAD", (-190.0, -30.0));
        assert_eq!(stats.n_voi_voxels, 10);
        assert_eq!(stats.n_fat_voxels, 5); // only the -80 HU voxels
        assert!((stats.fat_fraction - 0.5).abs() < 0.01);
    }

    #[test]
    fn test_histogram_bins() {
        let vol = Array3::<f32>::zeros((10, 10, 10));
        let mask = Array3::<bool>::default((10, 10, 10));

        let stats = compute_pcat_stats(&vol, &mask, "LAD", (-190.0, -30.0));

        assert_eq!(stats.histogram_bins.len(), 100);
        assert_eq!(stats.histogram_counts.len(), 100);
        // First bin center should be at -200 + 0.5 * 4 = -198
        assert!((stats.histogram_bins[0] - (-198.0)).abs() < 0.01);
    }

    /// Synthetic straight-cylinder contour (radius 3 mm along z), matching
    /// voi.rs's test fixture — used to exercise `compute_fai_analysis`.
    fn cylinder_contours() -> ContourResult {
        let (n_positions, n_angles, radius) = (16usize, 36usize, 3.0_f64);
        ContourResult {
            r_theta: vec![vec![radius; n_angles]; n_positions],
            r_eq: vec![radius; n_positions],
            positions_mm: (0..n_positions).map(|i| [(8 + i) as f64, 16.0, 16.0]).collect(),
            n_frame: vec![[0.0, 1.0, 0.0]; n_positions],
            b_frame: vec![[0.0, 0.0, 1.0]; n_positions],
            arclengths: (0..n_positions).map(|i| i as f64).collect(),
        }
    }

    #[test]
    fn fai_analysis_views_reconcile() {
        // Volume whose HU varies spatially but stays inside the fat window
        // everywhere, so every VOI voxel is counted (none dropped by the HU
        // filter). The angular sectors then partition EXACTLY the FAI fat
        // voxels — this is the invariant the dashboard depends on.
        let (nz, ny, nx) = (32usize, 32usize, 32usize);
        let mut vol = Array3::<f32>::zeros((nz, ny, nx));
        for ((_z, y, x), v) in vol.indexed_iter_mut() {
            // ∈ [-98, -82] HU across the volume → all within the fat range,
            // but varying by (y, x) so different sectors get different means.
            *v = (-90.0 + 0.4 * (x as f64 - 16.0) + 0.3 * (y as f64 - 16.0)) as f32;
        }
        let contours = cylinder_contours();
        let spacing = [1.0, 1.0, 1.0];

        let stats = compute_fai_analysis(
            &vol, &contours, spacing, "RCA",
            (-190.0, -30.0), 1.0, 3.0, 8, 20.0, 1.0,
        );

        let ang = stats.angular_asymmetry.expect("angular present");

        // Sectors partition exactly the FAI fat voxels (no voxel dropped or
        // double-counted).
        let sector_total: usize = ang.sectors.iter().map(|s| s.n_voxels).sum();
        assert_eq!(
            sector_total, stats.n_fat_voxels,
            "angular sectors must partition the FAI fat voxels exactly"
        );
        assert!(stats.n_fat_voxels > 100, "expected a populated VOI");

        // Count-weighted mean of the sectors == FAI mean (the consistency
        // guarantee that was the whole point of the unified pass).
        let wsum: f64 = ang
            .sectors
            .iter()
            .filter(|s| s.n_voxels > 0)
            .map(|s| s.hu_mean * s.n_voxels as f64)
            .sum();
        let weighted = wsum / sector_total as f64;
        assert!(
            (weighted - stats.hu_mean).abs() < 1e-6,
            "angular aggregate {weighted} must equal FAI mean {}",
            stats.hu_mean
        );

        // At least two sectors must actually differ, else the test is vacuous.
        let means: Vec<f64> = ang.sectors.iter().filter(|s| s.n_voxels > 0).map(|s| s.hu_mean).collect();
        let spread = means.iter().cloned().fold(f64::MIN, f64::max)
            - means.iter().cloned().fold(f64::MAX, f64::min);
        assert!(spread > 0.5, "sectors should vary (spread {spread})");

        // Radial profile: the rings inside the VOI band (1–4 mm) are populated.
        let radial = stats.radial_profile.expect("radial present");
        let in_band: usize = radial
            .distances_mm
            .iter()
            .zip(&radial.mean_hu)
            .filter(|(&d, m)| d > 1.0 && d <= 4.0 && !m.is_nan())
            .count();
        assert!(in_band >= 2, "VOI-band radial rings should be populated");
    }

    #[test]
    fn fai_analysis_empty_contour_is_safe() {
        let vol = Array3::<f32>::zeros((8, 8, 8));
        let empty = ContourResult {
            r_theta: vec![],
            r_eq: vec![],
            positions_mm: vec![],
            n_frame: vec![],
            b_frame: vec![],
            arclengths: vec![],
        };
        let stats = compute_fai_analysis(
            &vol, &empty, [1.0, 1.0, 1.0], "LAD",
            (-190.0, -30.0), 1.0, 3.0, 8, 20.0, 1.0,
        );
        assert_eq!(stats.n_fat_voxels, 0);
        assert_eq!(stats.fai_risk, "N/A");
    }
}
