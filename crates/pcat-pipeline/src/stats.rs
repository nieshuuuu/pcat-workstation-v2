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
}

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

    // Validate the contour contract at the boundary (golden rule 7): the four
    // per-section vectors must share one length and a uniform angle count — the
    // voxel sweep indexes nf[p]/bf[p]/r_theta[p] by section with no further
    // bounds checks. A ragged ContourResult is a producer bug, so fail loud.
    assert_eq!(
        contours.n_frame.len(),
        n_pos,
        "compute_fai_analysis: n_frame len {} != positions {n_pos}",
        contours.n_frame.len()
    );
    assert_eq!(
        contours.b_frame.len(),
        n_pos,
        "compute_fai_analysis: b_frame len {} != positions {n_pos}",
        contours.b_frame.len()
    );
    assert_eq!(
        contours.r_theta.len(),
        n_pos,
        "compute_fai_analysis: r_theta len {} != positions {n_pos}",
        contours.r_theta.len()
    );
    assert!(
        contours.r_theta.iter().all(|v| v.len() == n_angles),
        "compute_fai_analysis: ragged contour — every section must have {n_angles} angles"
    );

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

    // Per-z-slice accumulator for the parallel voxel sweep. Slices are summed
    // back in index order (collect → serial fold), so the FAI mean is identical
    // run-to-run — no nondeterministic reduction order leaks into the biomarker.
    use rayon::prelude::*;
    struct Acc {
        bin_sum: Vec<f64>,
        bin_sq: Vec<f64>,
        bin_cnt: Vec<usize>,
        sec_sum: Vec<f64>,
        sec_sq: Vec<f64>,
        sec_cnt: Vec<usize>,
        voi_values: Vec<f64>,
        fat_values: Vec<f64>,
    }
    impl Acc {
        fn zeros(n_radial: usize, n_sectors: usize) -> Self {
            Acc {
                bin_sum: vec![0.0; n_radial],
                bin_sq: vec![0.0; n_radial],
                bin_cnt: vec![0; n_radial],
                sec_sum: vec![0.0; n_sectors],
                sec_sq: vec![0.0; n_sectors],
                sec_cnt: vec![0; n_sectors],
                voi_values: Vec::new(),
                fat_values: Vec::new(),
            }
        }
        fn merge(mut self, other: Acc) -> Acc {
            for i in 0..self.bin_sum.len() {
                self.bin_sum[i] += other.bin_sum[i];
                self.bin_sq[i] += other.bin_sq[i];
                self.bin_cnt[i] += other.bin_cnt[i];
            }
            for s in 0..self.sec_sum.len() {
                self.sec_sum[s] += other.sec_sum[s];
                self.sec_sq[s] += other.sec_sq[s];
                self.sec_cnt[s] += other.sec_cnt[s];
            }
            self.voi_values.extend(other.voi_values);
            self.fat_values.extend(other.fat_values);
            self
        }
    }

    // One voxel pass over the bounding box, parallelised across z-slices. Each
    // voxel is assigned to its nearest contour section by an exact full scan;
    // the slices run concurrently, so there is no need to subsample the
    // centerline for speed (the old serial `cl_step` did, at the cost of
    // accuracy on long vessels). Projecting into that section's (N, B) frame
    // gives r/theta/dist_from_wall, exactly as `build_voi`.
    let partials: Vec<Acc> = (iz0..iz1 + 1)
        .into_par_iter()
        .map(|vz| {
            let mut acc = Acc::zeros(n_radial, n_sectors);
            for vy in iy0..=iy1 {
                for vx in ix0..=ix1 {
                    let vm = Vector3::new(
                        vz as f64 * spacing[0],
                        vy as f64 * spacing[1],
                        vx as f64 * spacing[2],
                    );
                    // Nearest section (exact full scan).
                    let (mut best, mut bestd) = (0usize, f64::MAX);
                    for i in 0..n_pos {
                        let d = (vm - pos[i]).norm_squared();
                        if d < bestd {
                            bestd = d;
                            best = i;
                        }
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
                        acc.bin_sum[b] += hu;
                        acc.bin_sq[b] += hu * hu;
                        acc.bin_cnt[b] += 1;
                    }

                    // CRISP-CT VOI shell → FAI + angular.
                    let in_voi = dist_from_wall > gap_mm && dist_from_wall <= gap_mm + ring_mm;
                    if in_voi {
                        acc.voi_values.push(hu);
                        if is_fat {
                            acc.fat_values.push(hu);
                            let anchored = {
                                let aa = (theta - theta_ref[p]) % tau;
                                if aa < 0.0 { aa + tau } else { aa }
                            };
                            let s = ((anchored / tau * n_sectors as f64) as usize) % n_sectors;
                            acc.sec_sum[s] += hu;
                            acc.sec_sq[s] += hu * hu;
                            acc.sec_cnt[s] += 1;
                        }
                    }
                }
            }
            acc
        })
        .collect();

    let Acc {
        bin_sum,
        bin_sq,
        bin_cnt,
        sec_sum,
        sec_sq,
        sec_cnt,
        voi_values,
        mut fat_values,
    } = partials
        .into_iter()
        .fold(Acc::zeros(n_radial, n_sectors), Acc::merge);

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
    let angular_asymmetry = AngularAsymmetry { sectors };

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

    /// L-shaped (curved) contour: a +z limb that bends into a +y limb, so the
    /// reconciliation invariants are exercised on non-straight geometry — the
    /// case the straight cylinder cannot reach (each voxel's nearest section,
    /// hence its VOI membership and sector, depends on the bend).
    fn curved_contours() -> ContourResult {
        let (n_angles, radius) = (36usize, 3.0_f64);
        let mut positions_mm = Vec::new();
        let mut n_frame = Vec::new();
        let mut b_frame = Vec::new();
        // Limb 1 runs along +z (frame N=+y, B=+x).
        for i in 0..8 {
            positions_mm.push([(8 + i) as f64, 16.0, 16.0]);
            n_frame.push([0.0, 1.0, 0.0]);
            b_frame.push([0.0, 0.0, 1.0]);
        }
        // Limb 2 runs along +y (frame N=+z, B=+x).
        for i in 1..8 {
            positions_mm.push([15.0, (16 + i) as f64, 16.0]);
            n_frame.push([1.0, 0.0, 0.0]);
            b_frame.push([0.0, 0.0, 1.0]);
        }
        let n_positions = positions_mm.len();
        ContourResult {
            r_theta: vec![vec![radius; n_angles]; n_positions],
            r_eq: vec![radius; n_positions],
            positions_mm,
            n_frame,
            b_frame,
            arclengths: (0..n_positions).map(|i| i as f64).collect(),
        }
    }

    #[test]
    fn fai_analysis_reconciles_on_curved_vessel() {
        let (nz, ny, nx) = (40usize, 40usize, 40usize);
        let mut vol = Array3::<f32>::zeros((nz, ny, nx));
        for ((_z, y, x), v) in vol.indexed_iter_mut() {
            *v = (-90.0 + 0.4 * (x as f64 - 16.0) + 0.3 * (y as f64 - 16.0)) as f32;
        }
        let contours = curved_contours();
        let stats = compute_fai_analysis(
            &vol, &contours, [1.0, 1.0, 1.0], "RCA",
            (-190.0, -30.0), 1.0, 3.0, 8, 20.0, 1.0,
        );
        let ang = stats.angular_asymmetry.expect("angular present");

        // Sectors still partition exactly the FAI fat voxels on a bent vessel.
        let sector_total: usize = ang.sectors.iter().map(|s| s.n_voxels).sum();
        assert_eq!(
            sector_total, stats.n_fat_voxels,
            "sectors must partition the FAI fat voxels on a curved vessel"
        );
        assert!(stats.n_fat_voxels > 50, "expected a populated curved VOI");

        // Count-weighted sector mean == FAI mean, the dashboard's consistency
        // guarantee, must hold on the bend too.
        let wsum: f64 = ang
            .sectors
            .iter()
            .filter(|s| s.n_voxels > 0)
            .map(|s| s.hu_mean * s.n_voxels as f64)
            .sum();
        assert!(
            (wsum / sector_total as f64 - stats.hu_mean).abs() < 1e-6,
            "angular aggregate must equal FAI mean on a curved vessel"
        );
    }

    #[test]
    fn fai_analysis_is_deterministic() {
        // The parallel voxel sweep must yield a bit-identical FAI mean run to
        // run — the collect→serial-fold reduction order does not depend on
        // thread scheduling.
        let (nz, ny, nx) = (32usize, 32usize, 32usize);
        let mut vol = Array3::<f32>::zeros((nz, ny, nx));
        for ((_z, _y, x), v) in vol.indexed_iter_mut() {
            *v = (-90.0 + 0.4 * (x as f64 - 16.0)) as f32;
        }
        let contours = cylinder_contours();
        let run = || {
            compute_fai_analysis(
                &vol, &contours, [1.0, 1.0, 1.0], "LAD",
                (-190.0, -30.0), 1.0, 3.0, 8, 20.0, 1.0,
            )
        };
        let a = run();
        let b = run();
        assert_eq!(a.n_fat_voxels, b.n_fat_voxels);
        assert_eq!(
            a.hu_mean.to_bits(),
            b.hu_mean.to_bits(),
            "FAI mean must be bit-identical across runs"
        );
    }
}
