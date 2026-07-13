//! Whole-volume water/lipid/PROTEIN (3-material) decomposition commands.
//!
//! The 3-material successor to `commands::water_lipid`. Where the 2-material GLS
//! self-calibrates the water→lipid line from the patient, this uses the FROZEN
//! poly2 calibration surface [`pcat_pipeline::mmd::WlpModel`] (baked from the
//! `wlp-decomposition` sim at 70/150 keV). `run_water_lipid_protein` validates
//! that the loaded pair is 70/150 and registers the baked model;
//! `get_wlp_slice` derives one axial slice's f_w/f_l/f_p/σ_f maps on demand.
//!
//! The surface is sim-calibrated — see the `WlpModel` module docs for the
//! sim→real transfer caveat. The 2-material GLS commands are untouched and still
//! back the pericoronary ROI overlay.

use std::sync::{Arc, Mutex};

use ndarray::Array3;
use serde::Serialize;
use tauri::ipc::Response;

use pcat_pipeline::mmd::{self, WlpModel};

use crate::commands::framed::encode_frame;
use crate::state::AppState;

/// The baked surface is fit at these energies; applying it to a different pair
/// is physically wrong (see `WlpModel` transfer caveat). Tolerance is generous
/// because MonoPlus keV are parsed from folder names (integer keV).
const WLP_LOW_KEV: f64 = 70.0;
const WLP_HIGH_KEV: f64 = 150.0;
const WLP_KEV_TOL: f64 = 1.0;

/// Register the baked WLP poly2 surface for the loaded dual-energy volume.
///
/// Fails loud if no dual-energy volume is loaded or the pair is not 70/150 keV
/// (the only energies the baked surface is valid for). Returns the model so the
/// frontend can show the endpoints + the sim-calibration caveat and size the
/// slice viewer. Runs instantly — the surface is a constant, not a scan.
#[tauri::command]
pub async fn run_water_lipid_protein(
    state: tauri::State<'_, Mutex<AppState>>,
) -> Result<WlpModel, String> {
    let (low_kev, high_kev, dims) = {
        let guard = state.lock().map_err(|e| format!("lock poisoned: {e}"))?;
        let de = guard.dual_energy.as_ref().ok_or_else(|| {
            "no dual-energy volume loaded — load a two-keV (MonoPlus 70/150) series first".to_string()
        })?;
        (
            de.low_energy_kev,
            de.high_energy_kev,
            [de.low.shape()[0], de.low.shape()[1], de.low.shape()[2]],
        )
    };

    if (low_kev - WLP_LOW_KEV).abs() > WLP_KEV_TOL || (high_kev - WLP_HIGH_KEV).abs() > WLP_KEV_TOL {
        return Err(format!(
            "the baked water/lipid/protein surface is calibrated for {WLP_LOW_KEV:.0}/{WLP_HIGH_KEV:.0} keV, \
             but the loaded pair is {low_kev:.0}/{high_kev:.0} keV. Refit the surface for this pair in \
             wlp-decomposition (apply_wlp_model.jl `fit_surface`) before using it here."
        ));
    }

    let mut model = WlpModel::baked();
    model.dims = dims; // record the grid this registration applies to (viewer sizing)
    {
        let mut guard = state.lock().map_err(|e| format!("lock poisoned: {e}"))?;
        guard.wlp_model = Some(model.clone());
    }
    Ok(model)
}

/// Restore a previously-saved WLP model into backend state (session reload) so
/// `get_wlp_slice` can derive maps without re-running. Mirrors
/// `restore_wl_calibration`.
#[tauri::command]
pub async fn restore_wlp_model(
    model: WlpModel,
    state: tauri::State<'_, Mutex<AppState>>,
) -> Result<(), String> {
    let mut guard = state.lock().map_err(|e| format!("lock poisoned: {e}"))?;
    guard.wlp_model = Some(model);
    Ok(())
}

/// Metadata header for a framed WLP slice response.
#[derive(Serialize)]
struct WlpSliceMeta {
    z: usize,
    ny: usize,
    nx: usize,
}

/// Derive one axial slice's CT + f_w + f_l + f_p + σ_f maps from the stored WLP
/// model.
///
/// Framed binary layout: `{z, ny, nx}` JSON header, then payload = `ny·nx` i16
/// CT HU (low energy, LE) ++ `ny·nx` f32 f_w ++ f_l ++ f_p ++ σ_f. Gated voxels
/// are `NaN` in every fraction plane. Unlike the 2-material slice, ALL THREE
/// fractions are sent (f_l is no longer 1 − f_w).
/// `smoothing` is the TV strength (λ): 0 ⇒ raw per-voxel maps (noisy but the
/// honest point estimate); ~8 (the frontend default) ⇒ WL-map smoothness. The
/// 3-material decode is ~10× noisier than the 2-material line, so a legible map
/// needs strong smoothing — this is a live knob because the ideal strength is
/// data-dependent.
#[tauri::command]
pub async fn get_wlp_slice(
    z: usize,
    smoothing: f64,
    state: tauri::State<'_, Mutex<AppState>>,
) -> Result<Response, String> {
    let (low, high, model) = {
        let guard = state.lock().map_err(|e| format!("lock poisoned: {e}"))?;
        let model = guard
            .wlp_model
            .clone()
            .ok_or_else(|| "no water/lipid/protein model — run the decomposition first".to_string())?;
        let de = guard
            .dual_energy
            .as_ref()
            .ok_or_else(|| "no dual-energy volume loaded".to_string())?;
        (Arc::clone(&de.low), Arc::clone(&de.high), model)
    };

    let framed = tokio::task::spawn_blocking(move || build_wlp_slice_frame(z, smoothing, &low, &high, &model))
        .await
        .map_err(|e| format!("get_wlp_slice task failed: {e}"))??;

    Ok(Response::new(framed))
}

/// Extract slice `z`, decompose it (poly2 + coupled TV at strength `smoothing`),
/// and pack the framed bytes.
fn build_wlp_slice_frame(
    z: usize,
    smoothing: f64,
    low: &Array3<f32>,
    high: &Array3<f32>,
    model: &WlpModel,
) -> Result<Vec<u8>, String> {
    let (nz, ny, nx) = (low.shape()[0], low.shape()[1], low.shape()[2]);
    if z >= nz {
        return Err(format!("slice index z={z} out of range (0..{nz})"));
    }
    let plane = ny * nx;
    let base = z * plane;

    let lo = low.as_slice().ok_or("low energy volume not contiguous")?;
    let hi = high.as_slice().ok_or("high energy volume not contiguous")?;
    let low_slice = &lo[base..base + plane];
    let high_slice = &hi[base..base + plane];

    let maps = mmd::decompose_slice_wlp(low_slice, high_slice, model, ny, nx, smoothing.max(0.0));

    // CT as i16 HU (low energy, clamped to the loader's range).
    let mut ct = Vec::<i16>::with_capacity(plane);
    for &v in low_slice {
        ct.push(v.round().clamp(-1024.0, 3071.0) as i16);
    }

    let mut payload = Vec::<u8>::with_capacity(plane * 2 + plane * 4 * 4);
    payload.extend_from_slice(bytemuck::cast_slice(&ct));
    payload.extend_from_slice(bytemuck::cast_slice(&maps.fw));
    payload.extend_from_slice(bytemuck::cast_slice(&maps.fl));
    payload.extend_from_slice(bytemuck::cast_slice(&maps.fp));
    payload.extend_from_slice(bytemuck::cast_slice(&maps.sf));

    encode_frame(&WlpSliceMeta { z, ny, nx }, &payload)
}
