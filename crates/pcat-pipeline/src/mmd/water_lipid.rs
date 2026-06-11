//! Noise-aware water/lipid (2-material) decomposition by GLS projection onto
//! the water→lipid line — ported 1:1 from `wl-noise-aware-mmd`
//! (`src/wl_decompose.jl` + `examples/decompose_57955439_calfree.jl`).
//!
//! Unlike the 3-material direct/PWSQS solvers in this module, this runs on the
//! WHOLE volume with no ROI, self-calibrates every slot from the patient's own
//! anatomy (no external phantom), and reports a per-voxel water fraction f_w
//! plus its standard deviation σ_f. The lipid fraction is `1 − f_w`.
//!
//! The estimator core (`gls_fw`, `sigma_cov`) is pure and deterministic; the
//! calibration scan is the only volume-wide pass. Slices are decomposed on
//! demand by the caller (cheap per-voxel solve) — there is no cached result
//! volume, so switching the lipid anchor costs nothing.

use ndarray::Array3;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use super::materials::{Material, MaterialLibrary};

/// Lipid endpoint choice. Both anchors share Σ, ρ, the noise lines and the
/// gate; only `HU_l` (hence the line direction `g`) differs.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum WlAnchor {
    /// Fat-referenced: `HU_l` = the patient's measured subcutaneous-fat ROI.
    /// `f_w = 0` at the patient's adipose, `1` at water. Robust, no LAC table.
    Adipose,
    /// Absolute: `HU_l` = pure-lipid `1000·(μ_lipid/μ_water − 1)` from the
    /// app's `MaterialLibrary`. `f_l = 1` ⇒ pure lipid; patient adipose < 1.
    Theoretical,
}

/// Self-measured calibration for the GLS water/lipid solver. Every field is
/// derived from the patient's own anatomy in the loaded dual-energy volume.
/// Small and serializable — this is the single source of truth; per-slice maps
/// are derived from it on demand.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WlCalibration {
    pub low_kev: f64,
    pub high_kev: f64,
    /// Water endpoint HU at (low, high). Rides the vendor water calibration: (0, 0).
    pub hu_w: [f64; 2],
    /// Lipid endpoint HU at (low, high) — measured subcutaneous fat ROI mean.
    pub hu_l_adipose: [f64; 2],
    /// Lipid endpoint HU at (low, high) — theoretical pure lipid from the LAC table.
    pub hu_l_theo: [f64; 2],
    /// Muscle reference HU at (low, high) — for the noise line only, not an endpoint.
    pub hu_mus: [f64; 2],
    /// Noise line σ(HU)=a·HU+b at the low energy.
    pub ab_low: [f64; 2],
    /// Noise line σ(HU)=a·HU+b at the high energy.
    pub ab_high: [f64; 2],
    /// In-ROI HU std at the fat ROI, (low, high) — the low-HU noise anchor.
    pub sigma_fat: [f64; 2],
    /// In-ROI HU std at the muscle ROI, (low, high) — the mid-HU noise anchor.
    pub sigma_mus: [f64; 2],
    /// Low/high-energy noise correlation, measured from the fat ROI residuals.
    pub rho: f64,
    /// Soft-tissue gate on the low-energy HU: [lo, hi]. Outside ⇒ NaN.
    pub gate: [f64; 2],
    pub n_fat: usize,
    pub n_mus: usize,
    /// Fraction of volume voxels that pass the gate (water/lipid mixtures).
    pub fraction_kept: f64,
    /// Volume dimensions (nz, ny, nx) so the frontend can size the slice viewer.
    pub dims: [usize; 3],
}

impl WlCalibration {
    /// Lipid endpoint HU at the requested anchor.
    pub fn hu_l(&self, anchor: WlAnchor) -> [f64; 2] {
        match anchor {
            WlAnchor::Adipose => self.hu_l_adipose,
            WlAnchor::Theoretical => self.hu_l_theo,
        }
    }
}

/// 2×2 noise covariance Σ(HU) from per-energy σ-lines and correlation ρ.
///
/// `σ_E = max(max(a_E·HU_E + b_E, 0), 1e-6)`; off-diagonal
/// `c = clamp(ρ, −0.999, 0.999)·σ_lo·σ_hi`. The 1e-6 floor keeps Σ
/// positive-definite even where a σ-line extrapolates ≤ 0; the ρ clamp keeps it
/// non-singular. Mirrors `sigma_cov` in `wl_decompose.jl` (pv = 0).
#[inline]
pub fn sigma_cov(hu_low: f64, hu_high: f64, ab_low: [f64; 2], ab_high: [f64; 2], rho: f64) -> [[f64; 2]; 2] {
    let s_lo = (ab_low[0] * hu_low + ab_low[1]).max(0.0).max(1e-6);
    let s_hi = (ab_high[0] * hu_high + ab_high[1]).max(0.0).max(1e-6);
    let c = rho.clamp(-0.999, 0.999) * s_lo * s_hi;
    [[s_lo * s_lo, c], [c, s_hi * s_hi]]
}

/// GLS projection of (HU_low, HU_high) onto the water→lipid line.
///
/// Returns `(f_w, sigma_f)`: the water-fraction estimate and its standard
/// deviation. `f̂_w = (gᵀΣ⁻¹d)/(gᵀΣ⁻¹g)`, `σ_f = √(1/(gᵀΣ⁻¹g))`, with
/// `d = HU − HU_l`, `g = HU_w − HU_l`. The 2×2 inverse is hand-coded. Mirrors
/// `gls_fw` in `wl_decompose.jl`.
#[inline]
pub fn gls_fw(hu: [f64; 2], hu_w: [f64; 2], hu_l: [f64; 2], sigma: [[f64; 2]; 2]) -> (f64, f64) {
    let d = [hu[0] - hu_l[0], hu[1] - hu_l[1]];
    let g = [hu_w[0] - hu_l[0], hu_w[1] - hu_l[1]];

    // Σ⁻¹ for the symmetric 2×2 [[a, b], [b, e]].
    let (a, b, e) = (sigma[0][0], sigma[0][1], sigma[1][1]);
    let det = a * e - b * b;
    // det > 0 by construction (PD), but guard against a degenerate Σ.
    let inv_det = if det.abs() < 1e-300 { 0.0 } else { 1.0 / det };
    // Σ⁻¹ = inv_det * [[e, -b], [-b, a]].
    let si = [[e * inv_det, -b * inv_det], [-b * inv_det, a * inv_det]];

    // Σ⁻¹·g and Σ⁻¹·d.
    let si_g = [si[0][0] * g[0] + si[0][1] * g[1], si[1][0] * g[0] + si[1][1] * g[1]];
    let si_d = [si[0][0] * d[0] + si[0][1] * d[1], si[1][0] * d[0] + si[1][1] * d[1]];

    let g_si_g = g[0] * si_g[0] + g[1] * si_g[1];
    let g_si_d = g[0] * si_d[0] + g[1] * si_d[1];

    let fhat = g_si_d / g_si_g;
    let var = 1.0 / g_si_g;
    (fhat, var.max(0.0).sqrt())
}

/// Theoretical pure-lipid endpoint HU at (low, high) from the LAC table:
/// `HU_l(E) = 1000·(μ_lipid(E)/μ_water(E) − 1)`.
fn theoretical_lipid_hu(lib: &MaterialLibrary) -> [f64; 2] {
    let lo = 1000.0 * (lib.lac_low(Material::Lipid) / lib.lac_low(Material::Water) - 1.0);
    let hi = 1000.0 * (lib.lac_high(Material::Lipid) / lib.lac_high(Material::Water) - 1.0);
    [lo, hi]
}

/// 1-voxel 4-neighbour erosion of a 2D boolean mask (ny rows × nx cols,
/// row-major). A cell survives only if it and its 4 edge-neighbours are all
/// true; border cells are dropped. Mirrors `erode1` in the Julia script.
fn erode1(mask: &[bool], ny: usize, nx: usize) -> Vec<bool> {
    let mut out = vec![false; ny * nx];
    if ny < 3 || nx < 3 {
        return out;
    }
    for j in 1..ny - 1 {
        for i in 1..nx - 1 {
            let idx = j * nx + i;
            out[idx] = mask[idx]
                && mask[idx - 1]
                && mask[idx + 1]
                && mask[idx - nx]
                && mask[idx + nx];
        }
    }
    out
}

/// Per-energy HU values of one ROI on one slice — kept grouped by slice so the
/// noise std can be detrended per slice (`roi_sigma`).
struct SliceRoi {
    fat_lo: Vec<f32>,
    fat_hi: Vec<f32>,
    mus_lo: Vec<f32>,
    mus_hi: Vec<f32>,
}

/// Pooled per-slice-detrended residual std: subtract each slice's own mean,
/// pool the residuals across slices, take the unbiased std. Because residuals
/// sum to zero within each slice, the pooled mean is exactly zero, so
/// `σ = √(Σ residual² / (N − 1))`. Mirrors `roi_sigma` in the Julia script.
fn roi_sigma(groups: &[&[f32]]) -> f64 {
    let mut sum_sq = 0.0_f64;
    let mut n = 0_usize;
    for g in groups {
        if g.is_empty() {
            continue;
        }
        let mean = g.iter().map(|&v| v as f64).sum::<f64>() / g.len() as f64;
        for &v in g.iter() {
            let r = v as f64 - mean;
            sum_sq += r * r;
        }
        n += g.len();
    }
    if n < 2 {
        return 0.0;
    }
    (sum_sq / (n - 1) as f64).sqrt()
}

/// Self-calibrate every GLS slot from the patient's own anatomy across the
/// whole dual-energy volume (cal-free). Returns an error if the volume is not
/// dual-energy or the fat/muscle ROIs are too small to anchor the line.
pub fn self_calibrate(
    low: &Array3<f32>,
    high: &Array3<f32>,
    lib: &MaterialLibrary,
) -> Result<WlCalibration, String> {
    if low.shape() != high.shape() {
        return Err(format!(
            "self_calibrate: low {:?} and high {:?} shapes differ",
            low.shape(),
            high.shape()
        ));
    }
    if (lib.low_kev - lib.high_kev).abs() < 1e-6 {
        return Err("self_calibrate: low and high energies are equal — not a dual-energy volume".into());
    }

    let (nz, ny, nx) = (low.shape()[0], low.shape()[1], low.shape()[2]);
    let lo = low.as_slice().ok_or("low energy volume not contiguous")?;
    let hi = high.as_slice().ok_or("high energy volume not contiguous")?;
    let plane = ny * nx;

    // Scan each slice in parallel: build fat/muscle masks (HU windows + body
    // mask), erode 1 voxel, collect surviving voxel HU per channel.
    let per_slice: Vec<SliceRoi> = (0..nz)
        .into_par_iter()
        .map(|z| {
            let base = z * plane;
            let slo = &lo[base..base + plane];
            let shi = &hi[base..base + plane];

            let mut fat_mask = vec![false; plane];
            let mut mus_mask = vec![false; plane];
            for k in 0..plane {
                let l = slo[k];
                let h = shi[k];
                let body = l > -250.0;
                fat_mask[k] = body && (-130.0..=-70.0).contains(&l) && (-110.0..=-50.0).contains(&h);
                mus_mask[k] = body && (25.0..=75.0).contains(&l) && (20.0..=70.0).contains(&h);
            }
            let fat_e = erode1(&fat_mask, ny, nx);
            let mus_e = erode1(&mus_mask, ny, nx);

            let mut roi = SliceRoi {
                fat_lo: Vec::new(),
                fat_hi: Vec::new(),
                mus_lo: Vec::new(),
                mus_hi: Vec::new(),
            };
            for k in 0..plane {
                if fat_e[k] {
                    roi.fat_lo.push(slo[k]);
                    roi.fat_hi.push(shi[k]);
                }
                if mus_e[k] {
                    roi.mus_lo.push(slo[k]);
                    roi.mus_hi.push(shi[k]);
                }
            }
            roi
        })
        .collect();

    // Aggregate.
    let n_fat: usize = per_slice.iter().map(|s| s.fat_lo.len()).sum();
    let n_mus: usize = per_slice.iter().map(|s| s.mus_lo.len()).sum();
    if n_fat < 100 {
        return Err(format!(
            "self_calibrate: fat ROI too small ({n_fat} voxels) — need clean subcutaneous fat in the volume"
        ));
    }
    if n_mus < 100 {
        return Err(format!(
            "self_calibrate: muscle ROI too small ({n_mus} voxels) — need muscle tissue for the noise line"
        ));
    }

    let mean = |sel: &dyn Fn(&SliceRoi) -> &Vec<f32>| -> f64 {
        let mut sum = 0.0_f64;
        let mut n = 0_usize;
        for s in &per_slice {
            let v = sel(s);
            sum += v.iter().map(|&x| x as f64).sum::<f64>();
            n += v.len();
        }
        if n == 0 {
            0.0
        } else {
            sum / n as f64
        }
    };

    let hu_l = [mean(&|s| &s.fat_lo), mean(&|s| &s.fat_hi)];
    let hu_mus = [mean(&|s| &s.mus_lo), mean(&|s| &s.mus_hi)];

    // Per-slice-detrended residual std per channel/ROI.
    let fat_lo_groups: Vec<&[f32]> = per_slice.iter().map(|s| s.fat_lo.as_slice()).collect();
    let fat_hi_groups: Vec<&[f32]> = per_slice.iter().map(|s| s.fat_hi.as_slice()).collect();
    let mus_lo_groups: Vec<&[f32]> = per_slice.iter().map(|s| s.mus_lo.as_slice()).collect();
    let mus_hi_groups: Vec<&[f32]> = per_slice.iter().map(|s| s.mus_hi.as_slice()).collect();
    let sigma_fat = [roi_sigma(&fat_lo_groups), roi_sigma(&fat_hi_groups)];
    let sigma_mus = [roi_sigma(&mus_lo_groups), roi_sigma(&mus_hi_groups)];

    // 2-point noise line σ(HU)=a·HU+b through (HU_fat, σ_fat) and (HU_mus, σ_mus).
    let abline = |hu1: f64, s1: f64, hu2: f64, s2: f64| -> [f64; 2] {
        let denom = hu2 - hu1;
        if denom.abs() < 1e-9 {
            // Degenerate (fat and muscle landed at the same HU): flat line.
            return [0.0, s1];
        }
        let a = (s2 - s1) / denom;
        [a, s1 - a * hu1]
    };
    let ab_low = abline(hu_l[0], sigma_fat[0], hu_mus[0], sigma_mus[0]);
    let ab_high = abline(hu_l[1], sigma_fat[1], hu_mus[1], sigma_mus[1]);

    // ρ from fat-ROI residuals about the global fat mean (HU_l). Because the
    // residual means are zero by construction, Pearson = Σ(d_lo·d_hi)/√(Σd_lo²·Σd_hi²).
    let (mut s_lh, mut s_ll, mut s_hh) = (0.0_f64, 0.0_f64, 0.0_f64);
    for s in &per_slice {
        for k in 0..s.fat_lo.len() {
            let dl = s.fat_lo[k] as f64 - hu_l[0];
            let dh = s.fat_hi[k] as f64 - hu_l[1];
            s_lh += dl * dh;
            s_ll += dl * dl;
            s_hh += dh * dh;
        }
    }
    let rho = if s_ll > 0.0 && s_hh > 0.0 {
        s_lh / (s_ll.sqrt() * s_hh.sqrt())
    } else {
        0.0
    };

    // Soft-tissue gate on the low-energy HU (uses the adipose endpoint).
    let gate = [hu_l[0] - 40.0, 150.0];

    // Fraction of volume voxels that pass the gate.
    let kept: usize = lo
        .par_iter()
        .filter(|&&l| l as f64 >= gate[0] && l as f64 <= gate[1])
        .count();
    let fraction_kept = kept as f64 / lo.len() as f64;

    Ok(WlCalibration {
        low_kev: lib.low_kev,
        high_kev: lib.high_kev,
        hu_w: [0.0, 0.0],
        hu_l_adipose: hu_l,
        hu_l_theo: theoretical_lipid_hu(lib),
        hu_mus,
        ab_low,
        ab_high,
        sigma_fat,
        sigma_mus,
        rho,
        gate,
        n_fat,
        n_mus,
        fraction_kept,
        dims: [nz, ny, nx],
    })
}

/// Decompose a single axial slice (row-major `ny·nx` HU at low and high energy)
/// into `(f_w, σ_f)` flat arrays. Voxels failing the soft-tissue gate are
/// `f32::NAN` in both outputs. The per-voxel solve is `sigma_cov` then
/// `gls_fw` at the chosen lipid anchor — identical to `decompose_volume` in the
/// Julia, applied one slice at a time.
pub fn decompose_slice(
    low_slice: &[f32],
    high_slice: &[f32],
    calib: &WlCalibration,
    anchor: WlAnchor,
) -> (Vec<f32>, Vec<f32>) {
    assert_eq!(low_slice.len(), high_slice.len(), "slice length mismatch");
    let hu_l = calib.hu_l(anchor);
    let n = low_slice.len();
    let mut fw = vec![f32::NAN; n];
    let mut sf = vec![f32::NAN; n];

    for k in 0..n {
        let l = low_slice[k] as f64;
        let h = high_slice[k] as f64;
        if l < calib.gate[0] || l > calib.gate[1] {
            continue; // gated out → NaN
        }
        let sigma = sigma_cov(l, h, calib.ab_low, calib.ab_high, calib.rho);
        let (f, s) = gls_fw([l, h], calib.hu_w, hu_l, sigma);
        fw[k] = f as f32;
        sf[k] = s as f32;
    }
    (fw, sf)
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::Array3;

    fn diag_sigma(s_lo: f64, s_hi: f64) -> [[f64; 2]; 2] {
        [[s_lo * s_lo, 0.0], [0.0, s_hi * s_hi]]
    }

    #[test]
    fn gls_recovers_water_endpoint() {
        // At HU = HU_w the water fraction must be 1.
        let hu_w = [0.0, 0.0];
        let hu_l = [-104.0, -81.0];
        let sigma = diag_sigma(11.0, 10.0);
        let (fw, sf) = gls_fw(hu_w, hu_w, hu_l, sigma);
        assert!((fw - 1.0).abs() < 1e-9, "expected f_w=1 at water, got {fw}");
        assert!(sf > 0.0, "sigma_f must be positive, got {sf}");
    }

    #[test]
    fn gls_recovers_lipid_endpoint() {
        // At HU = HU_l the water fraction must be 0.
        let hu_w = [0.0, 0.0];
        let hu_l = [-104.0, -81.0];
        let sigma = diag_sigma(11.0, 10.0);
        let (fw, _) = gls_fw(hu_l, hu_w, hu_l, sigma);
        assert!(fw.abs() < 1e-9, "expected f_w=0 at lipid, got {fw}");
    }

    #[test]
    fn gls_recovers_midpoint_on_line() {
        // A point exactly on the line at fraction t must recover f_w=t,
        // independent of Σ (GLS is exact for points on the line).
        let hu_w = [0.0, 0.0];
        let hu_l = [-104.0, -81.0];
        let t = 0.6;
        let hu = [
            hu_l[0] + t * (hu_w[0] - hu_l[0]),
            hu_l[1] + t * (hu_w[1] - hu_l[1]),
        ];
        let sigma = sigma_cov(hu[0], hu[1], [-0.012, 10.7], [-0.010, 10.2], 0.43);
        let (fw, _) = gls_fw(hu, hu_w, hu_l, sigma);
        assert!((fw - t).abs() < 1e-9, "expected f_w={t} on the line, got {fw}");
    }

    #[test]
    fn sigma_cov_is_symmetric_and_pd() {
        let s = sigma_cov(-100.0, -80.0, [-0.012, 10.7], [-0.010, 10.2], 0.9);
        assert_eq!(s[0][1], s[1][0], "Σ must be symmetric");
        let det = s[0][0] * s[1][1] - s[0][1] * s[1][0];
        assert!(det > 0.0, "Σ must be positive-definite, det={det}");
        // |ρ| is clamped below 1, so the correlation cannot make Σ singular.
        let s_extreme = sigma_cov(-100.0, -80.0, [-0.012, 10.7], [-0.010, 10.2], 1.0);
        let det_e = s_extreme[0][0] * s_extreme[1][1] - s_extreme[0][1] * s_extreme[1][0];
        assert!(det_e > 0.0, "Σ must stay PD at ρ→1, det={det_e}");
    }

    #[test]
    fn sigma_cov_floors_negative_line() {
        // A σ-line extrapolating below zero is floored to a tiny positive value.
        let s = sigma_cov(5000.0, 5000.0, [-0.012, 10.7], [-0.010, 10.2], 0.0);
        assert!(s[0][0] > 0.0 && s[1][1] > 0.0, "variances must stay positive");
    }

    #[test]
    fn self_calibrate_recovers_planted_fat_and_positive_slope() {
        // Build a synthetic volume with three plateaus: fat (~ -100 / -80),
        // muscle (~ 50 / 45) and water (0 / 0), each with a little spread so the
        // noise std is non-zero and the muscle has more spread than fat (so the
        // 2-point line slope is well defined).
        let (nz, ny, nx) = (4usize, 30usize, 30usize);
        let mut low = Array3::<f32>::zeros((nz, ny, nx));
        let mut high = Array3::<f32>::zeros((nz, ny, nx));
        // Deterministic pseudo-jitter so the test is reproducible.
        let jit = |k: usize, amp: f32| ((k % 7) as f32 - 3.0) / 3.0 * amp;
        for z in 0..nz {
            for y in 0..ny {
                for x in 0..nx {
                    let k = z * ny * nx + y * nx + x;
                    let (l, h) = if x < 10 {
                        (-100.0 + jit(k, 3.0), -80.0 + jit(k + 1, 3.0)) // fat
                    } else if x < 20 {
                        (50.0 + jit(k, 8.0), 45.0 + jit(k + 1, 8.0)) // muscle (wider)
                    } else {
                        (0.0 + jit(k, 2.0), 0.0 + jit(k + 1, 2.0)) // water
                    };
                    low[[z, y, x]] = l;
                    high[[z, y, x]] = h;
                }
            }
        }
        let lib = MaterialLibrary::naeotom_70_150();
        let calib = self_calibrate(&low, &high, &lib).expect("calibration should succeed");

        // Fat endpoint recovered near the planted plateau (erosion trims edges).
        assert!(
            (calib.hu_l_adipose[0] - (-100.0)).abs() < 6.0,
            "fat low endpoint off: {}",
            calib.hu_l_adipose[0]
        );
        assert!(calib.n_fat > 100 && calib.n_mus > 100);
        // Muscle has more spread than fat ⇒ positive noise-line slope (σ rises
        // toward muscle HU as |HU| decreases... here muscle HU > fat HU and
        // σ_mus > σ_fat, so slope a > 0).
        assert!(calib.sigma_mus[0] > calib.sigma_fat[0], "muscle should be noisier in this phantom");
        assert!(calib.ab_low[0] > 0.0, "expected positive low-energy noise slope, got {}", calib.ab_low[0]);
        // ρ is finite and in range.
        assert!(calib.rho.abs() <= 1.0, "rho out of range: {}", calib.rho);
        // Theoretical lipid endpoint is negative (lipid attenuates less than water).
        assert!(calib.hu_l_theo[0] < 0.0 && calib.hu_l_theo[1] < 0.0);
    }

    #[test]
    fn self_calibrate_errors_on_equal_energies() {
        let low = Array3::<f32>::zeros((2, 4, 4));
        let high = Array3::<f32>::zeros((2, 4, 4));
        let mut lib = MaterialLibrary::naeotom_70_150();
        lib.high_kev = 70.0;
        assert!(self_calibrate(&low, &high, &lib).is_err());
    }

    #[test]
    fn decompose_slice_gates_and_solves() {
        // 1×3 slice: a gas voxel (gated), a fat voxel (f_w≈0), a water voxel (f_w≈1).
        let calib = WlCalibration {
            low_kev: 70.0,
            high_kev: 150.0,
            hu_w: [0.0, 0.0],
            hu_l_adipose: [-104.0, -81.0],
            hu_l_theo: [-112.0, -81.0],
            hu_mus: [50.0, 45.0],
            ab_low: [-0.012, 10.7],
            ab_high: [-0.010, 10.2],
            sigma_fat: [11.0, 10.0],
            sigma_mus: [12.0, 11.0],
            rho: 0.43,
            gate: [-144.0, 150.0],
            n_fat: 1000,
            n_mus: 1000,
            fraction_kept: 0.3,
            dims: [1, 1, 3],
        };
        let low = [-900.0_f32, -104.0, 0.0];
        let high = [-900.0_f32, -81.0, 0.0];
        let (fw, sf) = decompose_slice(&low, &high, &calib, WlAnchor::Adipose);
        assert!(fw[0].is_nan() && sf[0].is_nan(), "gas voxel must be gated to NaN");
        assert!(fw[1].abs() < 1e-4, "fat voxel f_w≈0, got {}", fw[1]);
        assert!((fw[2] - 1.0).abs() < 1e-4, "water voxel f_w≈1, got {}", fw[2]);
        assert!(sf[1] > 0.0 && sf[2] > 0.0);
    }
}
