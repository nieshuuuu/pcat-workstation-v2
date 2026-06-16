//! Whole-volume noise-aware GLS water/lipid decomposition commands.
//!
//! `run_water_lipid` self-calibrates the loaded dual-energy volume (every slot
//! measured from the patient's own fat/muscle — no external phantom) and stores
//! the calibration. `get_wl_slice` derives one axial slice's f_w/σ_f maps on
//! demand from that calibration — there is no cached result volume, so an
//! anchor switch is free. The estimator is `pcat_pipeline::mmd::water_lipid`.

use std::sync::{Arc, Mutex};

use ndarray::Array3;
use tauri::ipc::Response;
use serde::Serialize;

use pcat_pipeline::mmd::{self, WlAnchor, WlCalibration};

use crate::commands::framed::encode_frame;
use crate::state::AppState;

/// Parse the anchor string the frontend sends into the solver enum.
fn parse_anchor(anchor: &str) -> Result<WlAnchor, String> {
    match anchor {
        "adipose" => Ok(WlAnchor::Adipose),
        "theoretical" => Ok(WlAnchor::Theoretical),
        other => Err(format!("unknown anchor '{other}': expected 'adipose' or 'theoretical'")),
    }
}

/// Self-calibrate and store the GLS water/lipid calibration for the loaded
/// dual-energy volume. Returns the measured calibration (endpoints, noise
/// lines, ρ, gate, dims) so the frontend can show the readout and size the
/// slice viewer. Run once per dual-energy load; slices are derived afterwards.
#[tauri::command]
pub async fn run_water_lipid(
    state: tauri::State<'_, Mutex<AppState>>,
) -> Result<WlCalibration, String> {
    // Pull the dual-energy handles out under the lock, then release it for the
    // CPU-heavy calibration scan.
    let (low, high, low_kev, high_kev) = {
        let guard = state.lock().map_err(|e| format!("lock poisoned: {e}"))?;
        let de = guard
            .dual_energy
            .as_ref()
            .ok_or_else(|| "no dual-energy volume loaded — load a two-keV (e.g. MonoPlus 70/150) series first".to_string())?;
        (
            Arc::clone(&de.low),
            Arc::clone(&de.high),
            de.low_energy_kev,
            de.high_energy_kev,
        )
    };

    let calib = tokio::task::spawn_blocking(move || {
        mmd::self_calibrate(&low, &high, low_kev, high_kev)
    })
    .await
    .map_err(|e| format!("run_water_lipid task failed: {e}"))??;

    {
        let mut guard = state.lock().map_err(|e| format!("lock poisoned: {e}"))?;
        guard.wl_calibration = Some(calib.clone());
    }

    Ok(calib)
}

/// Restore a previously-saved water/lipid calibration into backend state. Used
/// by the session loader so the whole-volume WL maps reappear on reopen without
/// re-running the self-calibration scan. Slices are still derived on demand by
/// `get_wl_slice`, which requires the dual-energy volume to be loaded (it is,
/// after a patient load that has a keV pair).
#[tauri::command]
pub async fn restore_wl_calibration(
    calib: WlCalibration,
    state: tauri::State<'_, Mutex<AppState>>,
) -> Result<(), String> {
    let mut guard = state.lock().map_err(|e| format!("lock poisoned: {e}"))?;
    guard.wl_calibration = Some(calib);
    Ok(())
}

/// Metadata header for a framed water/lipid slice response.
#[derive(Serialize)]
struct WlSliceMeta {
    z: usize,
    ny: usize,
    nx: usize,
    anchor: String,
}

/// Derive one axial slice's CT + f_w + σ_f maps from the stored calibration.
///
/// Framed binary layout: `{z, ny, nx, anchor}` JSON header, then payload =
/// `ny·nx` i16 CT HU (low energy, little-endian) ++ `ny·nx` f32 f_w ++ `ny·nx`
/// f32 σ_f. Gated-out voxels are `NaN` in f_w/σ_f. f_l = 1 − f_w is computed on
/// the frontend.
#[tauri::command]
pub async fn get_wl_slice(
    z: usize,
    anchor: String,
    state: tauri::State<'_, Mutex<AppState>>,
) -> Result<Response, String> {
    let anchor_enum = parse_anchor(&anchor)?;

    let (low, high, calib) = {
        let guard = state.lock().map_err(|e| format!("lock poisoned: {e}"))?;
        let calib = guard
            .wl_calibration
            .clone()
            .ok_or_else(|| "no water/lipid calibration — run the decomposition first".to_string())?;
        let de = guard
            .dual_energy
            .as_ref()
            .ok_or_else(|| "no dual-energy volume loaded".to_string())?;
        (Arc::clone(&de.low), Arc::clone(&de.high), calib)
    };

    let framed = tokio::task::spawn_blocking(move || build_slice_frame(z, anchor_enum, &anchor, &low, &high, &calib))
        .await
        .map_err(|e| format!("get_wl_slice task failed: {e}"))??;

    Ok(Response::new(framed))
}

/// Extract slice `z`, decompose it, and pack the framed response bytes.
fn build_slice_frame(
    z: usize,
    anchor_enum: WlAnchor,
    anchor_label: &str,
    low: &Array3<f32>,
    high: &Array3<f32>,
    calib: &WlCalibration,
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

    let (fw, sf) = mmd::decompose_slice(low_slice, high_slice, calib, anchor_enum);

    // CT as i16 HU (clamped to the loader's range), then f_w, then σ_f.
    let mut ct = Vec::<i16>::with_capacity(plane);
    for &v in low_slice {
        ct.push(v.round().clamp(-1024.0, 3071.0) as i16);
    }

    let mut payload = Vec::<u8>::with_capacity(plane * 2 + plane * 4 + plane * 4);
    payload.extend_from_slice(bytemuck::cast_slice(&ct));
    payload.extend_from_slice(bytemuck::cast_slice(&fw));
    payload.extend_from_slice(bytemuck::cast_slice(&sf));

    let meta = WlSliceMeta {
        z,
        ny,
        nx,
        anchor: anchor_label.to_string(),
    };
    encode_frame(&meta, &payload)
}
