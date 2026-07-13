//! Noise-aware water/lipid/PROTEIN (3-material) decomposition by a quadratic
//! calibration surface — ported 1:1 from `wlp-decomposition`
//! (`apply_wlp_model.jl`, itself the frozen export of `wlp_decomposition.jl`).
//!
//! This is the 3-material successor to the 2-material GLS in [`super::water_lipid`].
//! The decode is a quadratic surface `f = poly2(HU_lo, HU_hi)` fit against
//! known-composition rods with the three PURE endpoints pinned (equality-
//! constrained LS) so pure water/lipid/protein decode to exactly 100% of that
//! material, normalised to the water/lipid/protein simplex; ROI
//! means that scatter outside the simplex under noise are projected back with a
//! Mahalanobis (Σ⁻¹) MLE onto the triangle. Delivered maps add the coupled
//! σ_f-weighted edge-preserving Huber-TV (never Gaussian — that wipes the
//! fat↔muscle boundary, the PVAT signal).
//!
//! # Provenance & transfer caveat (READ BEFORE USING ON REAL DATA)
//!
//! Unlike the 2-material GLS, this surface is **NOT self-calibrating**. Its
//! coefficients are frozen to ONE acquisition chain: BasisSimulator 80/140 kVp
//! EICT → Cong water/iodine → FBP → VMI at 70/150 keV, on a stadium QRM-thorax.
//! The endpoints are the same theoretical NIST anchors the GLS "Theoretical"
//! anchor uses (water 0, lipid −111.7/−81.2, protein 270.6/290.1), so on the
//! water/lipid axis it broadly agrees with the validated GLS; the QUADRATIC
//! terms encode the sim's beam-hardening curvature and the PROTEIN axis, neither
//! of which is validated on real scanner data. The held-out CCC (0.982/0.990/
//! 0.989) is a SIM number. On real 70/150 VMI it is a reasonable first estimate
//! but its accuracy must be re-earned; a new scanner/chain should refit the
//! surface from known-fraction ROIs in `wlp-decomposition` (`apply_wlp_model.jl`
//! `fit_surface` → a fresh `wlp_model_<pair>.toml`), then re-bake [`WlpModel::baked`]
//! from it. The struct is `Deserialize`, so a refit model can also be loaded at
//! runtime once a second scanner actually needs it (YAGNI until then). The baked
//! noise (σ_lo ≈ 9.7, σ_hi ≈ 2.3 HU) is likewise the sim's; it only weights the
//! TV pool and the MLE metric, not the point decode.

use rayon::prelude::*;
use serde::{Deserialize, Serialize};

/// Default TV smoothing strength for the delivered map. The data weight is
/// MEDIAN-NORMALIZED (see `decompose_slice_wlp`) so λ is noise-scale-invariant:
/// it no longer needs re-tuning when the dose/scanner (hence the absolute noise)
/// changes. On real 70/150 data λ≈0.5–0.8 keeps rod boundaries crisp while
/// cleaning flat tissue; λ≳1.2 over-smooths (boundaries bleed, TV staircasing).
/// Edge preservation comes from the Huber `eps`. (Pre-normalization this was 8,
/// but that broke once noise self-calibration shrank the weight ~8×.)
pub const TV_LAMBDA_DEFAULT: f64 = 0.7;
pub const TV_ITERS_DEFAULT: usize = 150;
pub const TV_EPS_DEFAULT: f64 = 0.05;

/// The 6-term quadratic surface basis `[1, h_lo, h_hi, h_lo², h_hi², h_lo·h_hi]`.
/// Identical to `poly2` in `apply_wlp_model.jl`.
#[inline]
fn poly2(h_lo: f64, h_hi: f64) -> [f64; 6] {
    [1.0, h_lo, h_hi, h_lo * h_lo, h_hi * h_hi, h_lo * h_hi]
}

/// ∂poly2/∂h_lo and ∂poly2/∂h_hi — for the σ_f delta-method noise propagation.
#[inline]
fn dpoly_lo(h_lo: f64, h_hi: f64) -> [f64; 6] {
    [0.0, 1.0, 0.0, 2.0 * h_lo, 0.0, h_hi]
}
#[inline]
fn dpoly_hi(h_lo: f64, h_hi: f64) -> [f64; 6] {
    [0.0, 0.0, 1.0, 0.0, 2.0 * h_hi, h_lo]
}

#[inline]
fn dot6(c: &[f64; 6], x: &[f64; 6]) -> f64 {
    let mut s = 0.0;
    for i in 0..6 {
        s += c[i] * x[i];
    }
    s
}

/// `σ(HU) = c₀·HU² + c₁·HU + c₂` — the per-energy convex-quadratic noise smile.
#[inline]
fn quad_sigma(c: &[f64; 3], hu: f64) -> f64 {
    c[0] * hu * hu + c[1] * hu + c[2]
}

/// The frozen (or refit) water/lipid/protein calibration surface + noise model.
///
/// Serializable so it can be stored in `AppState` / a session exactly like
/// [`super::water_lipid::WlCalibration`]. The default [`WlpModel::baked`] is the
/// sim-validated surface; [`WlpModel::from_calibration`] refits it from measured
/// ROIs for a new scanner.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WlpModel {
    pub low_kev: f64,
    pub high_kev: f64,
    /// poly2 → f_w coefficients.
    pub cw: [f64; 6],
    /// poly2 → f_l coefficients.
    pub cl: [f64; 6],
    /// poly2 → f_p coefficients.
    pub cp: [f64; 6],
    /// Affine (linear) f_l coefficients `[c0, c_lo, c_hi]` — commutes with the
    /// recon PSF, kept for mass-conservation integrals (unused in the point map).
    pub cl_aff: [f64; 3],
    /// Endpoint HU at (low, high): water, lipid, protein.
    pub hu_w: [f64; 2],
    pub hu_l: [f64; 2],
    pub hu_p: [f64; 2],
    /// σ(HU) convex-quadratic coefficients, low / high energy.
    pub sigma_lo: [f64; 3],
    pub sigma_hi: [f64; 3],
    /// Inter-energy noise correlation.
    pub rho: f64,
    /// 2×2 noise covariance in HU space (for the MLE projection metric).
    pub sigma_hu: [[f64; 2]; 2],
    /// `mean(measured − linear-mix theory)` over the cal rods — de-biases the
    /// exact G-inverse used for out-of-triangle ROI means.
    pub bias_hu: [f64; 2],
    /// Soft-tissue gate on the HIGH-energy HU: [lo, hi]. Outside ⇒ NaN.
    pub gate: [f64; 2],
    /// Whether this surface is the frozen sim model (true) or a real-data refit.
    pub is_baked: bool,
    /// Whether the noise (`sigma_lo`/`sigma_hi`/`rho`/`sigma_hu`) was self-measured
    /// from the loaded volume (true) or is the sim's baked noise (false). The
    /// sim's σ₁₅₀≈2.3 is ~10× too low vs real NAEOTOM (~22 HU), so the σ_f readout
    /// + TV weight are only honest once this is true. The SURFACE stays frozen
    /// either way (self-calibrating it needs known W/L/P rods).
    #[serde(default)]
    pub noise_measured: bool,
    /// [nz, ny, nx] of the dual-energy grid this model was registered against,
    /// so the frontend can size the slice viewer. The surface itself is
    /// volume-independent; this is filled in at registration (baked ⇒ [0,0,0]).
    #[serde(default)]
    pub dims: [usize; 3],
}

impl WlpModel {
    /// The sim-validated 70/150 keV surface, baked from `wlp_model_70_150.toml`.
    /// These numbers are the single source of truth; re-export from
    /// `wlp-decomposition/wlp_decomposition.jl` if the sim chain changes.
    pub fn baked() -> Self {
        WlpModel {
            low_kev: 70.0,
            high_kev: 150.0,
            // Endpoint-pinned (barycentric-indicator at water/lipid/protein) so
            // pure materials decode to exactly 100% — see fit note in the notebook.
            cw: [
                1.0,
                0.03155876909952963,
                -0.03356858676164307,
                1.9516893554874072e-5,
                9.66840988347038e-5,
                -0.00011932712142043603,
            ],
            cl: [
                -6.416748650785238e-17,
                -0.02438851460716694,
                0.023030338830520646,
                -1.5686512041428118e-6,
                -5.177952563700287e-5,
                5.593330573449891e-5,
            ],
            cp: [
                1.8145058998874857e-17,
                -0.007170254492362232,
                0.010538247931122019,
                -1.7948242350714233e-5,
                -4.490457319768504e-5,
                6.339381568590412e-5,
            ],
            cl_aff: [
                0.009616319311225403,
                -0.023923284134564383,
                0.021944191288115065,
            ],
            hu_w: [0.0, 0.0],
            hu_l: [-111.69487998422547, -81.21299890609298],
            hu_p: [270.5833989652018, 290.09786893836923],
            sigma_lo: [0.0, 0.0021938436543037414, 9.72910391285317],
            sigma_hi: [3.5298319690056125e-5, -0.00026561525919214845, 2.229225748394469],
            rho: 0.7874343156500762,
            sigma_hu: [
                [94.1948266304702, 17.314344949204337],
                [17.314344949204337, 5.132824391416491],
            ],
            bias_hu: [-2.2625705360740223, 0.5131352796373362],
            gate: [-300.0, 250.0],
            is_baked: true,
            noise_measured: false,
            dims: [0, 0, 0],
        }
    }

    /// Replace the baked (sim) noise with real noise self-measured from the
    /// loaded volume. Real PCCT noise is ~material-uniform (path-dominated), so a
    /// single σ per energy over soft tissue is representative — hence flat σ(HU)
    /// (constant term only), not the sim's convex smile. The surface is untouched.
    pub fn set_measured_noise(&mut self, sigma_lo: f64, sigma_hi: f64, rho: f64) {
        let (s_lo, s_hi) = (sigma_lo.max(1e-3), sigma_hi.max(1e-3));
        let r = rho.clamp(-0.999, 0.999);
        self.sigma_lo = [0.0, 0.0, s_lo];
        self.sigma_hi = [0.0, 0.0, s_hi];
        self.rho = r;
        self.sigma_hu = [[s_lo * s_lo, r * s_lo * s_hi], [r * s_lo * s_hi, s_hi * s_hi]];
        self.noise_measured = true;
    }

    /// Endpoint matrix `G = [[hu_l0−hu_w0, hu_p0−hu_w0], [hu_l1−hu_w1, hu_p1−hu_w1]]`,
    /// mapping (f_l, f_p) → measured HU offset.
    #[inline]
    fn g_matrix(&self) -> [[f64; 2]; 2] {
        [
            [self.hu_l[0] - self.hu_w[0], self.hu_p[0] - self.hu_w[0]],
            [self.hu_l[1] - self.hu_w[1], self.hu_p[1] - self.hu_w[1]],
        ]
    }

    /// Mahalanobis metric `A = Gᵀ Σ⁻¹ G` on (f_l, f_p) for the simplex projection.
    fn a_metric(&self) -> [[f64; 2]; 2] {
        let g = self.g_matrix();
        let s = self.sigma_hu;
        let det = s[0][0] * s[1][1] - s[0][1] * s[1][0];
        let inv_det = if det.abs() < 1e-300 { 0.0 } else { 1.0 / det };
        // Σ⁻¹ = inv_det · [[s11, −s01], [−s10, s00]].
        let si = [
            [s[1][1] * inv_det, -s[0][1] * inv_det],
            [-s[1][0] * inv_det, s[0][0] * inv_det],
        ];
        // M = Σ⁻¹ G.
        let m = [
            [si[0][0] * g[0][0] + si[0][1] * g[1][0], si[0][0] * g[0][1] + si[0][1] * g[1][1]],
            [si[1][0] * g[0][0] + si[1][1] * g[1][0], si[1][0] * g[0][1] + si[1][1] * g[1][1]],
        ];
        // A = Gᵀ M.
        [
            [g[0][0] * m[0][0] + g[1][0] * m[1][0], g[0][0] * m[0][1] + g[1][0] * m[1][1]],
            [g[0][1] * m[0][0] + g[1][1] * m[1][0], g[0][1] * m[0][1] + g[1][1] * m[1][1]],
        ]
    }

    /// Raw per-voxel decode: three surfaces, normalised to sum 1. Mirrors
    /// `decode` in `apply_wlp_model.jl`. May land outside the simplex under
    /// noise (that is expected; do NOT rectify per-voxel before averaging — it
    /// would bias the ROI mean by Jensen). Returns (f_w, f_l, f_p).
    #[inline]
    pub fn decode(&self, hu_lo: f64, hu_hi: f64) -> [f64; 3] {
        let b = poly2(hu_lo, hu_hi);
        let x = dot6(&self.cw, &b);
        let y = dot6(&self.cl, &b);
        let z = dot6(&self.cp, &b);
        let s = x + y + z;
        [x / s, y / s, z / s]
    }

    /// Feasible decode for a POOLED ROI mean: interior → the validated surface;
    /// exterior → the bias-corrected G-inverse projected onto the simplex in the
    /// Mahalanobis metric. Mirrors `decode_feas` in `apply_wlp_model.jl`.
    pub fn decode_feas(&self, hu_lo: f64, hu_hi: f64) -> [f64; 3] {
        let b = poly2(hu_lo, hu_hi);
        let r = [dot6(&self.cw, &b), dot6(&self.cl, &b), dot6(&self.cp, &b)];
        if r.iter().all(|&v| v >= -1e-9) {
            let s = r[0] + r[1] + r[2];
            return [r[0] / s, r[1] / s, r[2] / s];
        }
        // θ = G \ ([hu] − bias); then MLE-project (f_l, f_p) onto the simplex.
        let g = self.g_matrix();
        let v = [hu_lo - self.bias_hu[0], hu_hi - self.bias_hu[1]];
        let det = g[0][0] * g[1][1] - g[0][1] * g[1][0];
        let (fl0, fp0) = if det.abs() < 1e-300 {
            (0.0, 0.0)
        } else {
            (
                (v[0] * g[1][1] - g[0][1] * v[1]) / det,
                (g[0][0] * v[1] - v[0] * g[1][0]) / det,
            )
        };
        let (fl, fp) = self.proj_simplex(fl0, fp0);
        [1.0 - fl - fp, fl, fp]
    }

    /// Mahalanobis-closest point on the water/lipid/protein simplex (triangle
    /// with corners water=(0,0), lipid=(1,0), protein=(0,1) in (f_l, f_p)).
    /// Mirrors `proj_simplex` in `apply_wlp_model.jl`.
    fn proj_simplex(&self, fl: f64, fp: f64) -> (f64, f64) {
        if fl >= -1e-9 && fp >= -1e-9 && (fl + fp) <= 1.0 + 1e-9 {
            return (fl, fp);
        }
        let a = self.a_metric();
        let theta = [fl, fp];
        let verts = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]];
        let edges = [(0usize, 1usize), (0, 2), (1, 2)];
        // A·v helper.
        let av = |v: [f64; 2]| [a[0][0] * v[0] + a[0][1] * v[1], a[1][0] * v[0] + a[1][1] * v[1]];
        let dot = |x: [f64; 2], y: [f64; 2]| x[0] * y[0] + x[1] * y[1];
        let mut best = verts[0];
        let mut bd = f64::INFINITY;
        for &(pi, qi) in &edges {
            let p = verts[pi];
            let q = verts[qi];
            let u = [q[0] - p[0], q[1] - p[1]];
            let tm = [theta[0] - p[0], theta[1] - p[1]];
            let denom = dot(u, av(u));
            let t = if denom.abs() < 1e-300 {
                0.0
            } else {
                (dot(u, av(tm)) / denom).clamp(0.0, 1.0)
            };
            let c = [p[0] + t * u[0], p[1] + t * u[1]];
            let dv = [c[0] - theta[0], c[1] - theta[1]];
            let d = dot(dv, av(dv));
            if d < bd {
                bd = d;
                best = c;
            }
        }
        (best[0], best[1])
    }

    /// σ_f data weight `1/(σ²_fl + σ²_fp)` from the delta method: propagate the
    /// per-energy σ smiles through the lipid+protein surface gradients with the
    /// inter-energy ρ. Mirrors `sigma_f_weight` in `apply_wlp_model.jl`. Used to
    /// weight the coupled TV pool (strong axis held, weak axis smoothed).
    #[inline]
    pub fn sigma_f_weight(&self, hu_lo: f64, hu_hi: f64) -> f64 {
        let s4 = quad_sigma(&self.sigma_lo, hu_lo);
        let s7 = quad_sigma(&self.sigma_hi, hu_hi);
        let dlo = dpoly_lo(hu_lo, hu_hi);
        let dhi = dpoly_hi(hu_lo, hu_hi);
        let g4l = dot6(&self.cl, &dlo);
        let g7l = dot6(&self.cl, &dhi);
        let g4p = dot6(&self.cp, &dlo);
        let g7p = dot6(&self.cp, &dhi);
        let vl = g4l * g4l * s4 * s4 + g7l * g7l * s7 * s7 + 2.0 * self.rho * g4l * g7l * s4 * s7;
        let vp = g4p * g4p * s4 * s4 + g7p * g7p * s7 * s7 + 2.0 * self.rho * g4p * g7p * s4 * s7;
        1.0 / (vl + vp).max(1e-6)
    }

    /// Per-voxel σ_f (std of the recovered lipid fraction) for the readout —
    /// √(variance of f_l via the same delta method). Reported alongside f_w/f_l/f_p.
    #[inline]
    pub fn sigma_fl(&self, hu_lo: f64, hu_hi: f64) -> f64 {
        let s4 = quad_sigma(&self.sigma_lo, hu_lo);
        let s7 = quad_sigma(&self.sigma_hi, hu_hi);
        let dlo = dpoly_lo(hu_lo, hu_hi);
        let dhi = dpoly_hi(hu_lo, hu_hi);
        let g4l = dot6(&self.cl, &dlo);
        let g7l = dot6(&self.cl, &dhi);
        let vl = g4l * g4l * s4 * s4 + g7l * g7l * s7 * s7 + 2.0 * self.rho * g4l * g7l * s4 * s7;
        vl.max(0.0).sqrt()
    }
}

/// Coupled edge-preserving Huber-TV on (f_l, f_p), σ_f-weighted. `mask` is the
/// gated soft-tissue mask (row-major ny·nx); `w` is the per-voxel data weight
/// (`sigma_f_weight`), `NaN`/non-finite ⇒ unit weight. Never Gaussian (that
/// wipes the fat↔muscle boundary). Mirrors `tv_coupled` in `apply_wlp_model.jl`.
/// Returns (f_l, f_p) with `NaN` outside the mask.
#[allow(clippy::too_many_arguments)]
pub fn tv_coupled(
    yl: &[f64],
    yp: &[f64],
    mask: &[bool],
    ny: usize,
    nx: usize,
    lambda: f64,
    iters: usize,
    eps: f64,
    w: &[f64],
) -> (Vec<f64>, Vec<f64>) {
    let n = ny * nx;
    let mut fl: Vec<f64> = (0..n).map(|k| if mask[k] { yl[k] } else { 0.0 }).collect();
    let mut fp: Vec<f64> = (0..n).map(|k| if mask[k] { yp[k] } else { 0.0 }).collect();
    let mut fl2 = fl.clone();
    let mut fp2 = fp.clone();
    // (di, dj) 4-neighbours.
    let neigh: [(isize, isize); 4] = [(1, 0), (-1, 0), (0, 1), (0, -1)];
    // Jacobi iterations: each output voxel reads only the PREVIOUS iterate, so
    // the sweep is embarrassingly parallel over voxels (rayon).
    for _ in 0..iters {
        let (fl_prev, fp_prev) = (&fl, &fp);
        fl2.par_iter_mut()
            .zip(fp2.par_iter_mut())
            .enumerate()
            .for_each(|(k, (out_l, out_p))| {
                if !mask[k] {
                    *out_l = fl_prev[k];
                    *out_p = fp_prev[k];
                    return;
                }
                let i = k % nx;
                let j = k / nx;
                let wij = if w[k].is_finite() { w[k] } else { 1.0 };
                let mut rl = wij * yl[k];
                let mut rp = wij * yp[k];
                let mut den = wij;
                for &(di, dj) in &neigh {
                    let ni = i as isize + di;
                    let nj = j as isize + dj;
                    if ni < 0 || nj < 0 || ni >= nx as isize || nj >= ny as isize {
                        continue;
                    }
                    let nk = (nj as usize) * nx + (ni as usize);
                    if !mask[nk] {
                        continue;
                    }
                    let dl = fl_prev[nk] - fl_prev[k];
                    let dp = fp_prev[nk] - fp_prev[k];
                    // Huber weight: small |∇| ⇒ strong smoothing (c≈λ/eps),
                    // large |∇| (a real fat↔muscle edge) ⇒ c collapses, edge kept.
                    let c = lambda / (dl * dl + dp * dp).sqrt().max(eps);
                    rl += c * fl_prev[nk];
                    rp += c * fp_prev[nk];
                    den += c;
                }
                // Project onto the non-negative simplex (a,b ≥ 0, a+b ≤ 1).
                let mut a = (rl / den).max(0.0);
                let mut b = (rp / den).max(0.0);
                let s = a + b;
                if s > 1.0 {
                    a /= s;
                    b /= s;
                }
                *out_l = a;
                *out_p = b;
            });
        std::mem::swap(&mut fl, &mut fl2);
        std::mem::swap(&mut fp, &mut fp2);
    }
    for k in 0..n {
        if !mask[k] {
            fl[k] = f64::NAN;
            fp[k] = f64::NAN;
        }
    }
    (fl, fp)
}

/// Self-measure the real per-energy noise `(σ_lo, σ_hi, ρ)` from the loaded
/// dual-energy volume over gated soft tissue, so the σ_f readout + TV weight
/// reflect the actual scanner rather than the sim.
///
/// Estimator: a robust MAD of the 4-neighbour Laplacian response (a high-pass
/// that kills smooth anatomy, leaving noise) over the central-third slices —
/// edge-insensitive (median), matches bcmmd-gaussian's `local_noise_map` to
/// ~10%. `σ = 1.4826·median(|∇²|)/√20` (√20 = the Laplacian's noise gain). ρ
/// from the correlation of the two channels' responses. Returns `None` if there
/// is too little soft tissue to measure (caller keeps the baked σ).
pub fn measure_noise(
    low: &ndarray::Array3<f32>,
    high: &ndarray::Array3<f32>,
    gate: [f64; 2],
) -> Option<(f64, f64, f64)> {
    let (nz, ny, nx) = (low.shape()[0], low.shape()[1], low.shape()[2]);
    if ny < 3 || nx < 3 || nz == 0 {
        return None;
    }
    let lo = low.as_slice()?;
    let hi = high.as_slice()?;
    let plane = ny * nx;
    let z0 = nz / 3;
    let z1 = (2 * nz / 3).max(z0 + 1).min(nz);
    let mut abs_lo = Vec::<f64>::new();
    let mut abs_hi = Vec::<f64>::new();
    let (mut sll, mut shh, mut slh) = (0.0_f64, 0.0_f64, 0.0_f64);
    for z in z0..z1 {
        let b = z * plane;
        for j in 1..ny - 1 {
            for i in 1..nx - 1 {
                let k = b + j * nx + i;
                let h = hi[k] as f64;
                if h <= gate[0] || h >= gate[1] {
                    continue;
                }
                let lap = |s: &[f32]| {
                    4.0 * s[k] as f64 - s[k - 1] as f64 - s[k + 1] as f64 - s[k - nx] as f64 - s[k + nx] as f64
                };
                let (ll, lh) = (lap(lo), lap(hi));
                abs_lo.push(ll.abs());
                abs_hi.push(lh.abs());
                sll += ll * ll;
                shh += lh * lh;
                slh += ll * lh;
            }
        }
    }
    if abs_lo.len() < 2000 {
        return None;
    }
    let median = |v: &mut Vec<f64>| {
        v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        v[v.len() / 2]
    };
    let gain = 20.0_f64.sqrt(); // std of the 4-neighbour Laplacian for white noise = σ·√20
    let s_lo = 1.4826 * median(&mut abs_lo) / gain;
    let s_hi = 1.4826 * median(&mut abs_hi) / gain;
    let rho = if sll > 0.0 && shh > 0.0 {
        slh / (sll.sqrt() * shh.sqrt())
    } else {
        0.0
    };
    Some((s_lo, s_hi, rho))
}

/// The four delivered per-voxel maps of one axial slice, row-major `ny·nx`.
/// Gated voxels are `NaN` in every plane. `fw + fl + fp = 1` where finite.
pub struct WlpMaps {
    pub fw: Vec<f32>,
    pub fl: Vec<f32>,
    pub fp: Vec<f32>,
    /// Per-voxel σ of the recovered lipid fraction (delta method) — the
    /// uncertainty proxy shown in the viewer. NOTE: from the model's baked
    /// (sim) noise, not the patient's — see the module transfer caveat.
    pub sf: Vec<f32>,
}

/// Decompose one axial slice into the delivered water/lipid/protein maps — the
/// exact `deliver()` pipeline of `wlp_decomposition.jl`: per-voxel poly2 decode
/// over the gated soft tissue, then the σ_f-weighted coupled Huber-TV pool.
///
/// Gating is on the HIGH-energy channel (`model.gate`, −300..250 HU at 150 keV
/// for the baked model) — matching the notebook, which gates on the high VMI.
/// The CT shown by the caller is the low-energy image; gate channel and display
/// channel are independent.
///
/// `tv_lambda` is the TV smoothing strength: `0` ⇒ raw per-voxel maps (noisy,
/// the honest point estimate); a positive value applies the coupled Huber-TV at
/// that strength ([`TV_LAMBDA_DEFAULT`] ≈ WL-map smoothness). The 3-material
/// per-voxel decode is intrinsically ~10× noisier than the 2-material line, so a
/// legible delivered map needs strong smoothing.
pub fn decompose_slice_wlp(
    low_slice: &[f32],
    high_slice: &[f32],
    model: &WlpModel,
    ny: usize,
    nx: usize,
    tv_lambda: f64,
) -> WlpMaps {
    assert_eq!(low_slice.len(), high_slice.len(), "slice length mismatch");
    assert_eq!(low_slice.len(), ny * nx, "slice length ≠ ny·nx");
    let n = ny * nx;

    let mask: Vec<bool> = (0..n)
        .map(|k| {
            let h = high_slice[k] as f64;
            h > model.gate[0] && h < model.gate[1]
        })
        .collect();

    let mut fl_raw = vec![0.0_f64; n];
    let mut fp_raw = vec![0.0_f64; n];
    let mut w = vec![1.0_f64; n];
    let mut sf = vec![f32::NAN; n];
    for k in 0..n {
        if !mask[k] {
            continue;
        }
        let l = low_slice[k] as f64;
        let h = high_slice[k] as f64;
        let d = model.decode(l, h);
        fl_raw[k] = d[1];
        fp_raw[k] = d[2];
        w[k] = model.sigma_f_weight(l, h);
        sf[k] = model.sigma_fl(l, h) as f32;
    }

    // Median-normalize the data weight over soft tissue so the TV strength λ is
    // independent of the absolute noise level (self-calibrating the noise shifts
    // the weight scale ~8×; without this λ would need re-tuning per dose). The
    // relative structure — noisier/weaker-axis voxels get less weight — is kept.
    {
        let mut wm: Vec<f64> = (0..n).filter(|&k| mask[k]).map(|k| w[k]).collect();
        if !wm.is_empty() {
            wm.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            let med = wm[wm.len() / 2];
            if med > 0.0 {
                for wk in w.iter_mut() {
                    *wk /= med;
                }
            }
        }
    }

    let (fl_d, fp_d) = if tv_lambda > 0.0 {
        tv_coupled(&fl_raw, &fp_raw, &mask, ny, nx, tv_lambda, TV_ITERS_DEFAULT, TV_EPS_DEFAULT, &w)
    } else {
        let fl: Vec<f64> = (0..n).map(|k| if mask[k] { fl_raw[k] } else { f64::NAN }).collect();
        let fp: Vec<f64> = (0..n).map(|k| if mask[k] { fp_raw[k] } else { f64::NAN }).collect();
        (fl, fp)
    };

    let mut fw = vec![f32::NAN; n];
    let mut fl = vec![f32::NAN; n];
    let mut fp = vec![f32::NAN; n];
    for k in 0..n {
        if mask[k] && fl_d[k].is_finite() && fp_d[k].is_finite() {
            fl[k] = fl_d[k] as f32;
            fp[k] = fp_d[k] as f32;
            fw[k] = (1.0 - fl_d[k] - fp_d[k]) as f32;
        }
    }
    WlpMaps { fw, fl, fp, sf }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn decompose_slice_gates_and_solves() {
        let m = WlpModel::baked();
        // 3 voxels: lung (gated on the high channel), water, real-adipose.
        let low = [-800.0_f32, 0.0, -104.0];
        let high = [-700.0_f32, 0.0, -81.0];
        let maps = decompose_slice_wlp(&low, &high, &m, 1, 3, 0.0); // 0 = raw, no TV
        assert!(
            maps.fw[0].is_nan() && maps.fl[0].is_nan() && maps.fp[0].is_nan(),
            "lung must gate to NaN"
        );
        // Water voxel: f_w dominates; fractions sum to 1 where finite.
        assert!(
            maps.fw[1] > maps.fl[1] && maps.fw[1] > maps.fp[1],
            "water: fw should dominate, got {:?}",
            (maps.fw[1], maps.fl[1], maps.fp[1])
        );
        for k in 1..3 {
            let s = maps.fw[k] + maps.fl[k] + maps.fp[k];
            assert!((s - 1.0).abs() < 1e-4, "sum≠1 at {k}: {s}");
        }
    }

    #[test]
    fn decode_matches_julia_golden() {
        // GOLDEN: the Rust port must reproduce `apply_wlp_model.jl`'s `decode` /
        // `decode_feas` bit-for-bit (proof the extraction is lossless). Values
        // captured from the verified Julia reference at these HU pairs. The
        // Endpoints are PINNED so pure materials decode to exactly 100% (water
        // (0,0)→f_w=1). Real adipose (−104,−81)→f_l≈0.785 (vs 0.69 for the old
        // unconstrained fit — pinning also straightens the near-pure regime where
        // real fat lives). Values captured from the verified Julia reference.
        let m = WlpModel::baked();
        let cases: [(f64, f64, [f64; 3], [f64; 3]); 4] = [
            (0.0, 0.0,
             [1.0, -6.41674865078524e-17, 1.81450589988749e-17],
             [1.0, -6.41674865078524e-17, 1.81450589988749e-17]),
            (-104.0, -81.0,
             [0.277170963640265, 0.785438242252227, -0.062609205892492],
             [0.0, 1.0, 0.0]),
            (50.0, 45.0,
             [1.04344356153402, -0.165985712507582, 0.122542150973565],
             [0.853655710372383, 0.0, 0.146344289627617]),
            (-60.0, -30.0,
             [0.0560191440695096, 0.820831944428277, 0.123148911502213],
             [0.0560191440695096, 0.820831944428277, 0.123148911502213]),
        ];
        for (lo, hi, want_dec, want_feas) in cases {
            let d = m.decode(lo, hi);
            let f = m.decode_feas(lo, hi);
            for k in 0..3 {
                assert!((d[k] - want_dec[k]).abs() < 1e-9, "decode({lo},{hi})[{k}]={} want {}", d[k], want_dec[k]);
                assert!((f[k] - want_feas[k]).abs() < 1e-9, "feas({lo},{hi})[{k}]={} want {}", f[k], want_feas[k]);
            }
        }
    }

    #[test]
    fn baked_endpoints_have_correct_argmax() {
        let m = WlpModel::baked();
        let argmax = |f: [f64; 3]| (0..3).max_by(|&a, &b| f[a].total_cmp(&f[b])).unwrap();
        assert_eq!(argmax(m.decode(m.hu_w[0], m.hu_w[1])), 0, "water endpoint");
        assert_eq!(argmax(m.decode(m.hu_l[0], m.hu_l[1])), 1, "lipid endpoint");
        assert_eq!(argmax(m.decode(m.hu_p[0], m.hu_p[1])), 2, "protein endpoint");
    }

    #[test]
    fn decode_normalises_to_one() {
        let m = WlpModel::baked();
        for &(a, b) in &[(-100.0, -80.0), (0.0, 0.0), (50.0, 45.0), (-50.0, -20.0)] {
            let f = m.decode(a, b);
            assert!((f[0] + f[1] + f[2] - 1.0).abs() < 1e-9, "sum≠1 at ({a},{b}): {f:?}");
        }
    }

    #[test]
    fn decode_feas_is_in_simplex() {
        let m = WlpModel::baked();
        // A deliberately off-triangle HU pair (well below the water/lipid edge).
        let f = m.decode_feas(-90.0, -60.0);
        assert!(f.iter().all(|&v| v >= -1e-9), "feasible decode left simplex: {f:?}");
        assert!((f[0] + f[1] + f[2] - 1.0).abs() < 1e-9, "sum≠1: {f:?}");
    }

    #[test]
    fn sigma_weight_positive() {
        let m = WlpModel::baked();
        assert!(m.sigma_f_weight(-100.0, -80.0) > 0.0);
        assert!(m.sigma_fl(-100.0, -80.0) > 0.0);
    }

    #[test]
    fn model_json_round_trips() {
        // The model crosses the Tauri IPC + session boundaries as JSON; a
        // round-trip must preserve every coefficient (else the frontend decode
        // silently diverges from the backend).
        let m = WlpModel::baked();
        let json = serde_json::to_string(&m).unwrap();
        let back: WlpModel = serde_json::from_str(&json).unwrap();
        let (a, b) = (back.decode(-104.0, -81.0), m.decode(-104.0, -81.0));
        for k in 0..3 {
            assert!((a[k] - b[k]).abs() < 1e-12, "decode diverged after round-trip: {a:?} vs {b:?}");
        }
        assert_eq!(back.gate, m.gate);
        assert_eq!(back.sigma_hu, m.sigma_hu);
        // A pre-dims session (no `dims` field) must still deserialize (serde default).
        let stripped = json.replace(",\"dims\":[0,0,0]", "");
        let old: WlpModel = serde_json::from_str(&stripped).unwrap();
        assert_eq!(old.dims, [0, 0, 0]);
    }

    #[test]
    fn measure_noise_recovers_planted_sigma() {
        // Plant known white noise (σ_lo=20, σ_hi=22, ρ≈0) on a smooth soft-tissue
        // background; the MAD-Laplacian estimator must recover it within ~15%.
        let (nz, ny, nx) = (9usize, 40usize, 40usize);
        let mut low = ndarray::Array3::<f32>::zeros((nz, ny, nx));
        let mut high = ndarray::Array3::<f32>::zeros((nz, ny, nx));
        // deterministic N(0,1) via Irwin-Hall (12 LCG uniforms − 6; mean 6, var 1).
        let gauss = |seed: u64| -> f64 {
            let mut s = seed.wrapping_mul(0x9E37_79B9_7F4A_7C15).wrapping_add(1);
            let mut sum = 0.0;
            for _ in 0..12 {
                s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                sum += ((s >> 40) as f64) / ((1u64 << 24) as f64);
            }
            sum - 6.0
        };
        for z in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let k = (z * ny * nx + j * nx + i) as u64;
                    let bg = 20.0 + 0.5 * (i as f64); // smooth ramp (killed by the Laplacian)
                    low[[z, j, i]] = (bg + 20.0 * gauss(2 * k)) as f32;
                    high[[z, j, i]] = (bg + 22.0 * gauss(2 * k + 1)) as f32;
                }
            }
        }
        let (s_lo, s_hi, _rho) = measure_noise(&low, &high, [-300.0, 250.0]).expect("should measure");
        assert!((s_lo - 20.0).abs() < 3.5, "σ_lo off: {s_lo}");
        assert!((s_hi - 22.0).abs() < 3.5, "σ_hi off: {s_hi}");
    }

    #[test]
    fn set_measured_noise_is_flat_and_flagged() {
        let mut m = WlpModel::baked();
        assert!(!m.noise_measured);
        m.set_measured_noise(22.0, 21.0, 0.71);
        assert!(m.noise_measured);
        assert_eq!(m.sigma_lo, [0.0, 0.0, 22.0]); // flat: material-uniform
        assert!((m.sigma_hu[0][0] - 484.0).abs() < 1e-9);
        assert!((m.sigma_fl(-100.0, -80.0) > 0.3), "real σ_fl should be ~2× the sim's ~0.2");
    }

    #[test]
    fn tv_denoises_but_preserves_edge_at_normalized_weight() {
        // `decompose_slice_wlp` median-normalizes the weight, so tv_coupled runs
        // with weight ≈ 1. In that regime the default λ must clean a flat region
        // (weak λ leaves the noise) AND keep a real edge. Build a noisy patch
        // (true f_l=0.6) split by a real edge (0.6 | 0.1) at unit weight.
        let (ny, nx) = (24usize, 24usize);
        let n = ny * nx;
        let mut yl = vec![0.0_f64; n];
        let yp = vec![0.05_f64; n];
        let mask = vec![true; n];
        let w = vec![1.0_f64; n]; // median-normalized weight
        let jit = |k: usize| ((k * 2654435761usize >> 13) & 0xff) as f64 / 255.0 - 0.5;
        for k in 0..n {
            let i = k % nx;
            let base = if i < nx / 2 { 0.60 } else { 0.10 };
            yl[k] = (base + 0.5 * jit(k)).clamp(0.0, 1.0);
        }
        let std_of = |v: &[f64], lo: usize, hi: usize| {
            let c: Vec<f64> = (0..n).filter(|&k| { let i = k % nx; i >= lo && i < hi }).map(|k| v[k]).collect();
            let m = c.iter().sum::<f64>() / c.len() as f64;
            (c.iter().map(|x| (x - m).powi(2)).sum::<f64>() / c.len() as f64).sqrt()
        };
        let mean_of = |v: &[f64], lo: usize, hi: usize| {
            let c: Vec<f64> = (0..n).filter(|&k| { let i = k % nx; i >= lo && i < hi }).map(|k| v[k]).collect();
            c.iter().sum::<f64>() / c.len() as f64
        };
        let noisy = std_of(&yl, 0, nx / 2);
        let weak = tv_coupled(&yl, &yp, &mask, ny, nx, 0.02, TV_ITERS_DEFAULT, TV_EPS_DEFAULT, &w);
        let strong = tv_coupled(&yl, &yp, &mask, ny, nx, TV_LAMBDA_DEFAULT, TV_ITERS_DEFAULT, TV_EPS_DEFAULT, &w);
        assert!(std_of(&weak.0, 0, nx / 2) > 0.5 * noisy, "weak λ should leave most noise");
        assert!(std_of(&strong.0, 0, nx / 2) < 0.35 * noisy, "default λ should clean the flat region");
        // The real edge survives (halves stay well separated).
        let sep = mean_of(&strong.0, 2, nx / 2 - 2) - mean_of(&strong.0, nx / 2 + 2, nx - 2);
        assert!(sep > 0.35, "edge preserved: separation {sep}");
    }

    #[test]
    fn tv_coupled_preserves_uniform_region() {
        // A uniform f_l=0.8, f_p=0.1 patch must survive TV unchanged (no edges).
        let (ny, nx) = (5, 5);
        let n = ny * nx;
        let yl = vec![0.8; n];
        let yp = vec![0.1; n];
        let mask = vec![true; n];
        let w = vec![1.0; n];
        let (fl, fp) = tv_coupled(&yl, &yp, &mask, ny, nx, 0.05, 10, 0.04, &w);
        // Interior voxel unchanged.
        let c = 2 * nx + 2;
        assert!((fl[c] - 0.8).abs() < 1e-6, "f_l drifted: {}", fl[c]);
        assert!((fp[c] - 0.1).abs() < 1e-6, "f_p drifted: {}", fp[c]);
    }
}
