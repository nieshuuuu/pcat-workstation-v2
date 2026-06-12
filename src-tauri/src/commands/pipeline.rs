use std::collections::HashMap;
use std::sync::Mutex;

use serde::{Deserialize, Serialize};
use tauri::Emitter;

use pcat_pipeline::centerline;
use pcat_pipeline::contour;
use pcat_pipeline::stats::FaiStats;
use crate::state::{AnalysisResults, AppState, Vessel, VesselResult};

/// Seed points and segment specification for a single vessel.
#[derive(Deserialize, Clone)]
pub struct VesselSeeds {
    /// Ostium (origin) point in mm [z, y, x].
    pub ostium_mm: [f64; 3],
    /// Waypoints along the vessel in mm [z, y, x].
    pub waypoints_mm: Vec<[f64; 3]>,
    /// Start of the proximal segment in mm from ostium.
    pub segment_start_mm: f64,
    /// Length of the proximal segment in mm.
    pub segment_length_mm: f64,
}

/// Progress event payload emitted during pipeline execution.
#[derive(Serialize, Clone)]
struct PipelineProgress {
    vessel: String,
    stage: String,
    progress: f64,
}

/// Parse a vessel name string into the Vessel enum.
fn parse_vessel(name: &str) -> Option<Vessel> {
    match name.to_uppercase().as_str() {
        "LAD" => Some(Vessel::LAD),
        "LCX" => Some(Vessel::LCx),
        "RCA" => Some(Vessel::RCA),
        _ => None,
    }
}

/// Build a dense centerline from ostium + waypoints using linear interpolation
/// at approximately `step_mm` spacing.
fn build_dense_centerline(
    ostium_mm: &[f64; 3],
    waypoints_mm: &[[f64; 3]],
    spacing: [f64; 3],
    origin: [f64; 3],
    direction: &[f64; 9],
    step_mm: f64,
) -> Vec<[f64; 3]> {
    // Collect all control points in order: ostium, then waypoints
    let mut control_pts: Vec<[f64; 3]> = Vec::with_capacity(1 + waypoints_mm.len());
    control_pts.push(*ostium_mm);
    control_pts.extend_from_slice(waypoints_mm);

    if control_pts.len() < 2 {
        return control_pts;
    }

    // Convert mm → voxel through the IOP-aware helper so non-axial acquisitions
    // (rare but possible) don't silently mis-locate the centerline.
    let inv_spacing = [1.0 / spacing[0], 1.0 / spacing[1], 1.0 / spacing[2]];
    let control_vox: Vec<[f64; 3]> = control_pts
        .iter()
        .map(|pt| pcat_pipeline::types::patient_to_voxel(*pt, origin, inv_spacing, direction))
        .collect();

    // Densely interpolate between consecutive control points at step_mm intervals
    let mut dense = Vec::new();
    dense.push(control_vox[0]);

    for seg in 0..control_vox.len() - 1 {
        let p0 = control_vox[seg];
        let p1 = control_vox[seg + 1];

        // Segment length in mm
        let dz = (p1[0] - p0[0]) * spacing[0];
        let dy = (p1[1] - p0[1]) * spacing[1];
        let dx = (p1[2] - p0[2]) * spacing[2];
        let seg_len_mm = (dz * dz + dy * dy + dx * dx).sqrt();

        let n_steps = (seg_len_mm / step_mm).ceil() as usize;
        if n_steps == 0 {
            continue;
        }

        for k in 1..=n_steps {
            let t = k as f64 / n_steps as f64;
            dense.push([
                p0[0] + t * (p1[0] - p0[0]),
                p0[1] + t * (p1[1] - p0[1]),
                p0[2] + t * (p1[2] - p0[2]),
            ]);
        }
    }

    dense
}

/// Run the full PCAT analysis pipeline for one or more vessels.
///
/// For each vessel:
/// 1. Build dense centerline from seed points (at ~0.5mm spacing)
/// 2. Clip to the proximal segment
/// 3. Estimate radii
/// 4. Extract contours (polar transform + gradient detection)
/// 5. Build perivascular VOI mask
/// 6. Compute FAI statistics
///
/// Emits `pipeline-progress` events via the Tauri event system.
#[tauri::command]
pub async fn run_pipeline(
    seeds: HashMap<String, VesselSeeds>,
    app: tauri::AppHandle,
    state: tauri::State<'_, Mutex<AppState>>,
) -> Result<HashMap<String, FaiStats>, String> {
    if seeds.is_empty() {
        return Err("no vessel seeds provided".into());
    }

    // Extract volume data under the lock, then release immediately
    let (volume_data, spacing, origin, direction) = {
        let guard = state.lock().map_err(|e| format!("lock poisoned: {e}"))?;
        let vol = guard
            .volume
            .as_ref()
            .ok_or_else(|| "no volume loaded".to_string())?;
        (vol.data.clone(), vol.spacing, vol.origin, vol.direction)
    };

    let seeds_clone = seeds.clone();
    let app_clone = app.clone();

    // Run heavy computation on a blocking thread
    let results = tokio::task::spawn_blocking(move || {
        let mut all_stats: HashMap<String, FaiStats> = HashMap::new();
        let total_vessels = seeds_clone.len();

        for (vessel_idx, (vessel_name, vessel_seeds)) in seeds_clone.iter().enumerate() {
            let vessel_base_progress = vessel_idx as f64 / total_vessels as f64;
            let vessel_progress_weight = 1.0 / total_vessels as f64;

            // Helper to emit progress
            let emit_progress = |stage: &str, stage_frac: f64| {
                let progress = vessel_base_progress + vessel_progress_weight * stage_frac;
                let _ = app_clone.emit(
                    "pipeline-progress",
                    PipelineProgress {
                        vessel: vessel_name.clone(),
                        stage: stage.to_string(),
                        progress,
                    },
                );
            };

            // --- Stage 1: Build dense centerline ---
            emit_progress("centerline", 0.0);

            let dense_centerline = build_dense_centerline(
                &vessel_seeds.ostium_mm,
                &vessel_seeds.waypoints_mm,
                spacing,
                origin,
                &direction,
                0.5, // 0.5mm step
            );

            if dense_centerline.len() < 2 {
                all_stats.insert(
                    vessel_name.clone(),
                    FaiStats {
                        vessel: vessel_name.clone(),
                        n_voi_voxels: 0,
                        n_fat_voxels: 0,
                        fat_fraction: 0.0,
                        hu_mean: 0.0,
                        hu_std: 0.0,
                        hu_median: 0.0,
                        fai_risk: "N/A".to_string(),
                        histogram_bins: vec![],
                        histogram_counts: vec![],
                        radial_profile: None,
                        angular_asymmetry: None,
                    },
                );
                continue;
            }

            // --- Stage 2: Clip to proximal segment ---
            emit_progress("clipping", 0.1);

            let clipped = centerline::clip_by_arclength(
                &dense_centerline,
                spacing,
                vessel_seeds.segment_start_mm,
                vessel_seeds.segment_length_mm,
            );

            if clipped.len() < 2 {
                all_stats.insert(
                    vessel_name.clone(),
                    FaiStats {
                        vessel: vessel_name.clone(),
                        n_voi_voxels: 0,
                        n_fat_voxels: 0,
                        fat_fraction: 0.0,
                        hu_mean: 0.0,
                        hu_std: 0.0,
                        hu_median: 0.0,
                        fai_risk: "N/A".to_string(),
                        histogram_bins: vec![],
                        histogram_counts: vec![],
                        radial_profile: None,
                        angular_asymmetry: None,
                    },
                );
                continue;
            }

            // --- Stage 4: Extract contours ---
            emit_progress("contours", 0.3);

            let contours = contour::extract_contours(
                &volume_data,
                &clipped,
                spacing,
                360,  // n_angles
                8.0,  // max_radius_mm
                5.0,  // sigma_deg
            );

            // --- Stage 5–8: FAI analysis ---
            // One voxel pass derives the FAI summary, the radial profile, and
            // the angular asymmetry from the SAME voxels, all referenced to the
            // half-max contour (per-angle) with a 1 mm wall gap. The three
            // dashboard views therefore reconcile by construction: the angular
            // sectors and the in-band radial rings average back to the FAI mean.
            emit_progress("stats", 0.85);

            let stats = pcat_pipeline::stats::compute_fai_analysis(
                &volume_data,
                &contours,
                spacing,
                vessel_name,
                (-190.0, -30.0), // FAI fat HU window
                1.0,             // gap_mm (CRISP-CT: skip the vessel wall)
                3.0,             // ring_mm
                8,               // n_sectors
                20.0,            // max_distance_mm (radial-profile extent)
                1.0,             // ring_step_mm
            );

            emit_progress("done", 1.0);
            all_stats.insert(vessel_name.clone(), stats);
        }

        all_stats
    })
    .await
    .map_err(|e| format!("pipeline task failed: {e}"))?;

    // Store results in AppState
    {
        let mut guard = state.lock().map_err(|e| format!("lock poisoned: {e}"))?;
        let mut vessel_results = HashMap::new();
        for (name, stats) in &results {
            if let Some(vessel) = parse_vessel(name) {
                vessel_results.insert(
                    vessel,
                    VesselResult {
                        fai_mean_hu: stats.hu_mean,
                        fai_risk: stats.fai_risk.clone(),
                        fat_fraction: stats.fat_fraction,
                        n_voi_voxels: stats.n_voi_voxels,
                        n_fat_voxels: stats.n_fat_voxels,
                        hu_std: stats.hu_std,
                        hu_median: stats.hu_median,
                        histogram_bins: stats.histogram_bins.clone(),
                        histogram_counts: stats.histogram_counts.clone(),
                        radial_profile: stats.radial_profile.clone(),
                        angular_asymmetry: stats.angular_asymmetry.clone(),
                    },
                );
            }
        }
        guard.analysis_results = Some(AnalysisResults {
            vessels: vessel_results,
        });
    }

    Ok(results)
}
