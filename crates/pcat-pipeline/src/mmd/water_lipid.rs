//! Noise-aware water/lipid (2-material) decomposition by GLS projection onto
//! the water→lipid line — ported 1:1 from `wl-noise-aware-mmd`
//! (`src/wl_decompose.jl` + `examples/decompose_57955439_calfree.jl`).
//!
//! This is the sole material-decomposition solver (the ill-conditioned
//! 3-material direct/PWSQS solvers were dropped). It runs on the WHOLE volume
//! with no ROI, self-calibrates every slot from the patient's own anatomy (no
//! external phantom), and reports a per-voxel water fraction f_w plus its
//! standard deviation σ_f. The lipid fraction is `1 − f_w`.
//!
//! The estimator core (`gls_fw`, `sigma_cov`) is pure and deterministic; the
//! calibration scan is the only volume-wide pass. The whole-volume viewer
//! decomposes slices on demand from the calibration; `decompose_volume_gls`
//! packs a masked-ROI result into the shared `MmdResult` for the pericoronary
//! tool. No result volume is cached, so switching the lipid anchor is free.

use ndarray::Array3;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use super::result::MmdResult;

/// Lipid endpoint choice. Both anchors share Σ, ρ, the noise lines and the
/// gate; only `HU_l` (hence the line direction `g`) differs.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum WlAnchor {
    /// Fat-referenced: `HU_l` = the patient's measured subcutaneous-fat ROI.
    /// `f_w = 0` at the patient's adipose, `1` at water. Robust, no table.
    Adipose,
    /// Absolute: `HU_l` = pure-triglyceride lipid from the baked NIST curve
    /// (`THEO_LIPID_HU`). `f_l = 1` ⇒ pure lipid; patient adipose reads < 1.
    Theoretical,
}

/// BasisSimulator NIST pure-triglyceride lipid endpoint `HU_l(E) =
/// 1000·(μ_lipid(E)/μ_water(E) − 1)` vs energy (keV). Baked from
/// wl-noise-aware-mmd `theoretical_endpoints` (src/wl_analysis.jl) so the
/// workstation reproduces the reference theoretical anchor without a Julia /
/// BasisSimulator dependency. Pure triglyceride is not adipose tissue, so this
/// is its own curve — distinct from the patient-measured adipose endpoint
/// (`hu_l_adipose`). Verified: 70→−111.69, 150→−81.21 match the reference TOML.
const THEO_LIPID_HU: [(f64, f64); 7] = [
    (40.0, -212.7195),
    (60.0, -129.1293),
    (70.0, -111.6949),
    (80.0, -101.5788),
    (100.0, -90.4788),
    (120.0, -85.0432),
    (150.0, -81.2130),
];

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
    /// Lipid endpoint HU at (low, high) — measured subcutaneous-fat ROI mean.
    pub hu_l_adipose: [f64; 2],
    /// Lipid endpoint HU at (low, high) — theoretical pure triglyceride (baked NIST curve).
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

/// Linear interpolation of the theoretical pure-lipid HU endpoint at `energy`,
/// clamped to the table range.
fn theo_lipid_hu_at(energy: f64) -> f64 {
    let first = THEO_LIPID_HU[0];
    let last = THEO_LIPID_HU[THEO_LIPID_HU.len() - 1];
    if energy <= first.0 {
        return first.1;
    }
    if energy >= last.0 {
        return last.1;
    }
    for w in THEO_LIPID_HU.windows(2) {
        let (e0, v0) = w[0];
        let (e1, v1) = w[1];
        if energy >= e0 && energy <= e1 {
            let t = (energy - e0) / (e1 - e0);
            return v0 * (1.0 - t) + v1 * t;
        }
    }
    last.1
}

/// Theoretical pure-lipid endpoint HU at (low, high) keV from the baked NIST
/// triglyceride curve.
fn theoretical_lipid_hu(low_kev: f64, high_kev: f64) -> [f64; 2] {
    [theo_lipid_hu_at(low_kev), theo_lipid_hu_at(high_kev)]
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
/// whole dual-energy volume (cal-free). `low_kev`/`high_kev` are the VMI
/// energies (used for the theoretical anchor + reporting). Returns an error if
/// the volume is not dual-energy or the fat/muscle ROIs are too small to anchor
/// the line.
pub fn self_calibrate(
    low: &Array3<f32>,
    high: &Array3<f32>,
    low_kev: f64,
    high_kev: f64,
) -> Result<WlCalibration, String> {
    if low.shape() != high.shape() {
        return Err(format!(
            "self_calibrate: low {:?} and high {:?} shapes differ",
            low.shape(),
            high.shape()
        ));
    }
    if (low_kev - high_kev).abs() < 1e-6 {
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
        low_kev,
        high_kev,
        hu_w: [0.0, 0.0],
        hu_l_adipose: hu_l,
        hu_l_theo: theoretical_lipid_hu(low_kev, high_kev),
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

/// Decompose the masked voxels of a whole volume into water/lipid fractions by
/// GLS, packed into the shared `MmdResult` so the existing pericoronary ROI
/// tool (overlay, surface plot, CSV) renders f_w/f_l directly. Iodine/calcium
/// are zero (2-material model). `f_w` is NOT clamped — muscle/fibrous voxels
/// read f_w > 1 (the binary model's contamination signature). The drawn ROI is
/// the selection, so the soft-tissue gate is NOT applied here. Unmasked voxels
/// are zero.
pub fn decompose_volume_gls(
    low_energy: &Array3<f32>,
    high_energy: &Array3<f32>,
    mask: &Array3<bool>,
    calib: &WlCalibration,
    anchor: WlAnchor,
) -> MmdResult {
    let shape = low_energy.shape();
    assert_eq!(shape, high_energy.shape(), "low/high energy shape mismatch");
    assert_eq!(shape, mask.shape(), "mask shape mismatch");
    let dim = (shape[0], shape[1], shape[2]);
    let hu_l = calib.hu_l(anchor);

    // mg/mL from the canonical densities in materials.rs (g/cm^3 × 1000).
    let rho_w = (super::materials::DENSITY_WATER * 1000.0) as f32;
    let rho_l = (super::materials::DENSITY_LIPID * 1000.0) as f32;

    let lo = low_energy.as_slice().expect("low_energy not contiguous");
    let hi = high_energy.as_slice().expect("high_energy not contiguous");
    let mk = mask.as_slice().expect("mask not contiguous");

    let results: Vec<[f32; 5]> = lo
        .par_iter()
        .zip(hi.par_iter())
        .zip(mk.par_iter())
        .map(|((&l, &h), &m)| {
            if !m {
                return [0.0; 5];
            }
            let sigma = sigma_cov(l as f64, h as f64, calib.ab_low, calib.ab_high, calib.rho);
            let (fw, _sf) = gls_fw([l as f64, h as f64], calib.hu_w, hu_l, sigma);
            let fw = fw as f32;
            let fl = 1.0 - fw;
            [fw, fl, fw * rho_w, fl * rho_l, fw * rho_w + fl * rho_l]
        })
        .collect();

    let n = results.len();
    let mut wf = vec![0.0f32; n];
    let mut lf = vec![0.0f32; n];
    let mut wm = vec![0.0f32; n];
    let mut lm = vec![0.0f32; n];
    let mut td = vec![0.0f32; n];
    for (i, r) in results.iter().enumerate() {
        wf[i] = r[0];
        lf[i] = r[1];
        wm[i] = r[2];
        lm[i] = r[3];
        td[i] = r[4];
    }

    MmdResult {
        water_frac: Array3::from_shape_vec(dim, wf).unwrap(),
        lipid_frac: Array3::from_shape_vec(dim, lf).unwrap(),
        iodine_frac: Array3::<f32>::zeros(dim),
        calcium_frac: Array3::<f32>::zeros(dim),
        water_mass: Array3::from_shape_vec(dim, wm).unwrap(),
        lipid_mass: Array3::from_shape_vec(dim, lm).unwrap(),
        iodine_mass: Array3::<f32>::zeros(dim),
        calcium_mass: Array3::<f32>::zeros(dim),
        total_density: Array3::from_shape_vec(dim, td).unwrap(),
        mask: mask.clone(),
        iterations: 1,
        converged: true,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::{array, Array3};

    fn diag_sigma(s_lo: f64, s_hi: f64) -> [[f64; 2]; 2] {
        [[s_lo * s_lo, 0.0], [0.0, s_hi * s_hi]]
    }

    #[test]
    fn gls_recovers_water_endpoint() {
        let hu_w = [0.0, 0.0];
        let hu_l = [-104.0, -81.0];
        let sigma = diag_sigma(11.0, 10.0);
        let (fw, sf) = gls_fw(hu_w, hu_w, hu_l, sigma);
        assert!((fw - 1.0).abs() < 1e-9, "expected f_w=1 at water, got {fw}");
        assert!(sf > 0.0, "sigma_f must be positive, got {sf}");
    }

    #[test]
    fn gls_recovers_lipid_endpoint() {
        let hu_w = [0.0, 0.0];
        let hu_l = [-104.0, -81.0];
        let sigma = diag_sigma(11.0, 10.0);
        let (fw, _) = gls_fw(hu_l, hu_w, hu_l, sigma);
        assert!(fw.abs() < 1e-9, "expected f_w=0 at lipid, got {fw}");
    }

    #[test]
    fn gls_recovers_midpoint_on_line() {
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
        let s_extreme = sigma_cov(-100.0, -80.0, [-0.012, 10.7], [-0.010, 10.2], 1.0);
        let det_e = s_extreme[0][0] * s_extreme[1][1] - s_extreme[0][1] * s_extreme[1][0];
        assert!(det_e > 0.0, "Σ must stay PD at ρ→1, det={det_e}");
    }

    #[test]
    fn sigma_cov_floors_negative_line() {
        let s = sigma_cov(5000.0, 5000.0, [-0.012, 10.7], [-0.010, 10.2], 0.0);
        assert!(s[0][0] > 0.0 && s[1][1] > 0.0, "variances must stay positive");
    }

    #[test]
    fn theoretical_endpoint_matches_reference_70_150() {
        // Baked BasisSimulator NIST triglyceride curve must reproduce the
        // reference TOML's theoretical pure-lipid endpoint at 70/150 keV.
        let hu = theoretical_lipid_hu(70.0, 150.0);
        assert!((hu[0] - (-111.6949)).abs() < 1e-3, "theo lipid 70keV: {}", hu[0]);
        assert!((hu[1] - (-81.2130)).abs() < 1e-3, "theo lipid 150keV: {}", hu[1]);
        // Interpolation between table points (e.g. 110 keV) lands between neighbours.
        let mid = theo_lipid_hu_at(110.0);
        assert!(mid < -85.0432 && mid > -90.4788, "110 keV interp out of bracket: {mid}");
    }

    #[test]
    fn self_calibrate_recovers_planted_fat_and_positive_slope() {
        let (nz, ny, nx) = (4usize, 30usize, 30usize);
        let mut low = Array3::<f32>::zeros((nz, ny, nx));
        let mut high = Array3::<f32>::zeros((nz, ny, nx));
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
        let calib = self_calibrate(&low, &high, 70.0, 150.0).expect("calibration should succeed");
        assert!((calib.hu_l_adipose[0] - (-100.0)).abs() < 6.0, "fat low endpoint off: {}", calib.hu_l_adipose[0]);
        assert!(calib.n_fat > 100 && calib.n_mus > 100);
        assert!(calib.sigma_mus[0] > calib.sigma_fat[0], "muscle should be noisier in this phantom");
        assert!(calib.ab_low[0] > 0.0, "expected positive low-energy noise slope, got {}", calib.ab_low[0]);
        assert!(calib.rho.abs() <= 1.0, "rho out of range: {}", calib.rho);
        // Theoretical endpoint is the baked NIST value at 70/150, not the planted fat.
        assert!((calib.hu_l_theo[0] - (-111.6949)).abs() < 1e-3);
        assert!((calib.hu_l_theo[1] - (-81.2130)).abs() < 1e-3);
    }

    #[test]
    fn self_calibrate_errors_on_equal_energies() {
        let low = Array3::<f32>::zeros((2, 4, 4));
        let high = Array3::<f32>::zeros((2, 4, 4));
        assert!(self_calibrate(&low, &high, 70.0, 70.0).is_err());
    }

    fn demo_calib() -> WlCalibration {
        WlCalibration {
            low_kev: 70.0,
            high_kev: 150.0,
            hu_w: [0.0, 0.0],
            hu_l_adipose: [-104.0, -81.0],
            hu_l_theo: [-111.69, -81.21],
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
        }
    }

    #[test]
    fn decompose_slice_gates_and_solves() {
        let calib = demo_calib();
        let low = [-900.0_f32, -104.0, 0.0];
        let high = [-900.0_f32, -81.0, 0.0];
        let (fw, sf) = decompose_slice(&low, &high, &calib, WlAnchor::Adipose);
        assert!(fw[0].is_nan() && sf[0].is_nan(), "gas voxel must be gated to NaN");
        assert!(fw[1].abs() < 1e-4, "fat voxel f_w≈0, got {}", fw[1]);
        assert!((fw[2] - 1.0).abs() < 1e-4, "water voxel f_w≈1, got {}", fw[2]);
        assert!(sf[1] > 0.0 && sf[2] > 0.0);
    }

    #[test]
    fn decompose_volume_gls_fills_masked_fractions() {
        // 1×1×2: one masked water voxel (f_w≈1), one unmasked (zero).
        let calib = demo_calib();
        let low = array![[[0.0_f32, 0.0]]];
        let high = array![[[0.0_f32, 0.0]]];
        let mask = array![[[true, false]]];
        let r = decompose_volume_gls(&low, &high, &mask, &calib, WlAnchor::Adipose);
        assert!((r.water_frac[[0, 0, 0]] - 1.0).abs() < 1e-3, "masked water f_w≈1");
        assert!(r.lipid_frac[[0, 0, 0]].abs() < 1e-3, "masked water f_l≈0");
        assert_eq!(r.water_frac[[0, 0, 1]], 0.0, "unmasked voxel is zero");
        assert_eq!(r.iodine_frac[[0, 0, 0]], 0.0);
        assert!(r.converged);
    }
}
