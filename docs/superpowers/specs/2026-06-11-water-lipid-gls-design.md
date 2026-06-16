# Whole-Volume Water/Lipid GLS Decomposition — Design

**Date:** 2026-06-11
**Goal:** Port the noise-aware GLS water/lipid quantification from
`wl-noise-aware-mmd` (the method that produced
`data/naeotom/analysis_57955439/fw_fl_maps_theolipid_57955439.png`) into
pcat-workstation-v2 as a **whole-volume axial viewer**, replacing the previous
material-decomposition experience (a tiny CPR ring-ROI you "can't see anything"
in).

## Why the old view fails

`run_mmd_on_roi` builds a 10 mm pericoronary ring mask from CPR contours and
runs a 3-material Cramer's-rule solve only inside that mask, rendered on a CPR
cross-section. There is no whole-volume picture, and the 3-material direct
inversion is ill-conditioned at soft-tissue contrast. The reference figure is
the opposite: every soft-tissue voxel of an **axial** slab, water/lipid only,
jet 0..1 over the CT.

## Algorithm (ported 1:1 from `src/wl_decompose.jl` + `examples/decompose_57955439_calfree.jl`)

Two-material noise-aware GLS projection onto the water→lipid line.

- **GLS estimate** (`gls_fw`): with `d = HU − HU_l`, `g = HU_w − HU_l`,
  `f̂_w = (gᵀΣ⁻¹d)/(gᵀΣ⁻¹g)`, `σ_f² = 1/(gᵀΣ⁻¹g)`. 2×2 inverse hand-coded.
- **Noise covariance** (`sigma_cov`): per-energy `σ_E = a_E·HU + b_E` (floored
  ≥0, then ≥1e-6), off-diagonal `c = clamp(ρ,−0.999,0.999)·σ_lo·σ_hi`.
- **Self-calibration (cal-free)** — every slot measured from the patient's own
  anatomy over the whole volume (no external phantom):
  - Body mask `HU_low > −250`.
  - Fat ROI `−130≤HU_lo≤−70 ∧ −110≤HU_hi≤−50 ∧ body`; muscle ROI
    `25≤HU_lo≤75 ∧ 20≤HU_hi≤70 ∧ body`. 1-voxel 4-neighbour erosion per slice.
  - `HU_l` (adipose endpoint) = mean fat ROI per channel; `HU_w = (0,0)`.
  - Noise line per channel = 2-point line through `(HU_fat, σ_fat)` and
    `(HU_mus, σ_mus)`, where `σ` is the pooled per-slice-detrended residual std
    (`roi_sigma`).
  - `ρ` = Pearson of fat-ROI residuals `(HU_lo−HU_l_lo, HU_hi−HU_l_hi)`.
- **Gate**: keep `HU_low ≥ HU_l_adipose_lo − 40 ∧ HU_low ≤ 150`; else NaN
  (drops gas/lung/iodine/calcium/bone). Gate uses the adipose endpoint for
  **both** anchors (matches the Julia).
- **Anchors**: `Adipose` (measured fat, fat-referenced scale) and `Theoretical`
  (pure-lipid `HU_l(E)=1000·(μ_lipid(E)/μ_water(E)−1)` from the app's existing
  `MaterialLibrary` LAC table — SSoT within the app). Both share Σ, ρ, gate,
  noise lines; only `HU_l` (hence `g`) changes.

## Architecture

SSoT: store only the **calibration** (a few scalars). Slices are **derived on
demand** — the per-voxel GLS solve is cheap, so `get_wl_slice(z, anchor)`
decomposes just that slice. No giant cached f_w/σ_f volume (honours the
no-cache preference; anchor switch is free).

### Backend (Rust)

- `crates/pcat-pipeline/src/mmd/water_lipid.rs` (new):
  `WlCalibration` (serializable), `WlAnchor`, `gls_fw`, `sigma_cov`,
  `self_calibrate(low, high, lib) -> Result<WlCalibration>`,
  `decompose_slice(low_slice, high_slice, &calib, anchor) -> (Vec<f32> fw, Vec<f32> sigma_f)`.
  Parallel calibration scan over slices (rayon); fail loud if fat/muscle ROI < 100 vox.
- `mmd/mod.rs`: export the new items.
- `src-tauri/src/state.rs`: add `wl_calibration: Option<WlCalibration>`.
- `src-tauri/src/commands/water_lipid.rs` (new): `run_water_lipid` (self-calibrate
  the loaded `dual_energy`, store calibration, return summary) and
  `get_wl_slice(z, anchor)` (framed binary: `{z,ny,nx}` + i16 CT + f32 f_w + f32 σ_f).
- register in `commands/mod.rs` + `lib.rs`.

### Frontend (Svelte)

- `src/lib/api.ts`: `runWaterLipid()`, `getWlSlice(z, anchor)` (framed decode),
  `WlCalibration`/`WlSlice` types.
- `src/components/WaterLipidView.svelte` (new): "Run", slice slider, material
  toggle (f_w / f_l / σ_f), anchor toggle (adipose / theoretical), CT window
  −160..240, jet 0..1 over grayscale with NaN transparent, anterior-up; colorbar
  + calibration readout (endpoints, ρ, noise lines, % kept).
- `src/App.svelte`: third tab "Water/Lipid".

## Testing

Rust unit tests: `gls_fw` recovers f_w=1 at the water endpoint and f_w=0 at the
lipid endpoint; synthetic 60/40 mixture on the line recovers 0.6; `sigma_cov`
symmetric/PD; `self_calibrate` on a synthetic fat+muscle+water volume returns
HU_l≈ planted fat mean and a positive noise slope. `cargo test -p pcat-pipeline`
+ `cargo check` for the Tauri crate; `npm run build` for the frontend.

## Out of scope

The non-contrast single-energy lever-rule + Huber-TV MAP path and the Δ
(theo−adipose) diverging panel from the Julia script — the deliverable is the
whole-volume f_w/f_l jet viewer. f_l is `1−f_w` client-side.
