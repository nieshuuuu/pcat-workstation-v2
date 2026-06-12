use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};

use tauri::{AppHandle, Emitter, Manager, State};
use tauri_plugin_dialog::DialogExt;
use tauri::ipc::Response;

use ndarray::Array3;
use pcat_pipeline::dicom_scan::{self, SeriesDescriptor};
use pcat_pipeline::dicom_load::{self, LoadedVolume as PipelineLoadedVolume, VolumeMetadata as PipelineVolumeMetadata};
use pcat_pipeline::types::LoadedVolume as StateLoadedVolume;
use crate::state::AppState;
use crate::volume_cache::CachedVolume;
use crate::commands::framed::encode_frame;

const MAX_RECENT: usize = 10;

fn recent_paths_file(app: &tauri::AppHandle) -> PathBuf {
    let dir = app.path().app_data_dir().expect("app data dir");
    dir.join("recent_dicoms.json")
}

fn load_recent_list(app: &tauri::AppHandle) -> Vec<String> {
    let path = recent_paths_file(app);
    if let Ok(data) = std::fs::read_to_string(&path) {
        serde_json::from_str(&data).unwrap_or_default()
    } else {
        Vec::new()
    }
}

fn save_recent_list(app: &tauri::AppHandle, paths: &[String]) {
    let path = recent_paths_file(app);
    if let Some(parent) = path.parent() {
        let _ = std::fs::create_dir_all(parent);
    }
    let _ = std::fs::write(&path, serde_json::to_string_pretty(paths).unwrap_or_default());
}

fn push_recent(app: &tauri::AppHandle, new_path: &str) {
    let mut list = load_recent_list(app);
    list.retain(|p| p != new_path);
    list.insert(0, new_path.to_string());
    list.truncate(MAX_RECENT);
    save_recent_list(app, &list);
}

/// Opens a native folder-picker dialog. Returns the selected path as a string,
/// or `None` if the user cancelled.
#[tauri::command]
pub async fn open_dicom_dialog(app: tauri::AppHandle) -> Result<Option<String>, String> {
    let path = tokio::task::spawn_blocking(move || {
        app.dialog().file().blocking_pick_folder()
    })
    .await
    .map_err(|e| format!("dialog task failed: {e}"))?;

    Ok(path.map(|p| p.to_string()))
}

/// Return the list of recently opened DICOM folder paths.
#[tauri::command]
pub async fn get_recent_dicoms(app: tauri::AppHandle) -> Result<Vec<String>, String> {
    Ok(load_recent_list(&app))
}

/// Sanitize a path string into a safe filename component.
fn sanitize_for_filename(s: &str) -> String {
    s.chars()
        .map(|c| if c.is_ascii_alphanumeric() || c == '-' || c == '_' || c == '.' { c } else { '_' })
        .collect::<String>()
        .replace("..", "_")
}

/// Seed-file key derived from the **full** DICOM folder path, so two patients
/// whose folders end in the same last component (e.g. `.../001/DICOM` and
/// `.../002/DICOM`) never collide. The last path component is prefixed for
/// human browseability — the sanitized full path after it guarantees uniqueness.
fn seeds_filename(dicom_path: &str) -> String {
    let short = Path::new(dicom_path)
        .file_name()
        .map(|f| f.to_string_lossy().to_string())
        .unwrap_or_default();
    let full = sanitize_for_filename(dicom_path);
    if short.is_empty() {
        format!("{}.json", full)
    } else {
        format!("{}__{}.json", sanitize_for_filename(&short), full)
    }
}

/// Save seeds JSON keyed by the full DICOM folder path.
#[tauri::command]
pub async fn save_seeds(app: tauri::AppHandle, seeds_json: String, dicom_path: String) -> Result<String, String> {
    let dir = app.path().app_data_dir().expect("app data dir").join("seeds");
    let _ = std::fs::create_dir_all(&dir);
    let path = dir.join(seeds_filename(&dicom_path));
    std::fs::write(&path, &seeds_json).map_err(|e| format!("write failed: {e}"))?;
    Ok(path.to_string_lossy().to_string())
}

/// Load seeds JSON for the given DICOM folder path.
///
/// Keyed by the full sanitized path so last-component collisions cannot
/// return another patient's seeds. Note: files saved before the collision
/// fix lived under a last-component-only name; those are intentionally
/// not consulted here, because a legacy file may contain seeds from a
/// *different* folder that happened to share the same final component —
/// serving those would reproduce the very bug this change fixes.
#[tauri::command]
pub async fn load_seeds(app: tauri::AppHandle, dicom_path: String) -> Result<Option<String>, String> {
    let dir = app.path().app_data_dir().expect("app data dir").join("seeds");
    let path = dir.join(seeds_filename(&dicom_path));
    if path.exists() {
        let data = std::fs::read_to_string(&path).map_err(|e| format!("read failed: {e}"))?;
        Ok(Some(data))
    } else {
        Ok(None)
    }
}

// ---------------------------------------------------------------------------
// Patient browser: list patient folders + status from saved annotations
// ---------------------------------------------------------------------------

/// Per-patient progress summary for the patient browser.
#[derive(serde::Serialize)]
pub struct PatientInfo {
    /// Folder name (e.g. "57955439"), used as a stable patient ID.
    pub id: String,
    /// Absolute path to the patient's DICOM folder.
    pub path: String,
    /// `not_started` | `in_progress` | `complete`.
    pub status: String,
    /// Number of cross-sections marked finalized in saved annotations (0 if none).
    pub finalized_count: usize,
    /// Whether MMD has been run and stored in saved annotations.
    pub has_mmd: bool,
}

/// Decide whether a directory looks like a patient folder.
///
/// Heuristic: name is non-hidden, contains at least one regular file (we don't
/// scan deeply for DICOM headers — that would be slow over SMB).
fn looks_like_patient_dir(entry: &std::fs::DirEntry) -> bool {
    if !entry.file_type().map(|t| t.is_dir()).unwrap_or(false) {
        return false;
    }
    let name = entry.file_name();
    let name = name.to_string_lossy();
    if name.starts_with('.') || name.starts_with('_') {
        return false;
    }
    // Quick check: directory must contain at least one regular file (DICOM slice).
    if let Ok(mut iter) = std::fs::read_dir(entry.path()) {
        iter.any(|e| e.ok().and_then(|e| e.file_type().ok()).is_some_and(|t| t.is_file()))
    } else {
        false
    }
}

/// True if any seeds JSON file has been saved under the given patient folder.
/// Seeds are saved per-series (`seeds/{series}__{sanitized-full-path}.json`),
/// so we scan the seeds directory for any filename that contains the
/// sanitized patient path — any hit means the user placed + saved seeds on
/// at least one series of this patient.
fn patient_has_seeds(app: &tauri::AppHandle, patient_path: &str) -> bool {
    let seeds_dir = app.path().app_data_dir().expect("app data dir").join("seeds");
    let Ok(read) = std::fs::read_dir(&seeds_dir) else { return false };
    let needle = sanitize_for_filename(patient_path);
    if needle.is_empty() { return false }
    for entry in read.flatten() {
        let name = entry.file_name().to_string_lossy().to_string();
        if name.contains(&needle) {
            return true;
        }
    }
    false
}

/// Reads saved annotation JSON for a given patient folder name and returns
/// (finalized_count, has_mmd). Returns (0, false) if no save exists.
fn read_annotation_summary(app: &tauri::AppHandle, folder_name: &str) -> (usize, bool) {
    let dir = app.path().app_data_dir().expect("app data dir").join("annotations");
    // Filename is `sanitize(folder_name).json` — match the convention in
    // `save_annotations` / `load_annotations`.
    let file = dir.join(format!("{}.json", sanitize_for_filename(folder_name)));
    let Ok(data) = std::fs::read_to_string(&file) else {
        return (0, false);
    };
    let Ok(json) = serde_json::from_str::<serde_json::Value>(&data) else {
        return (0, false);
    };
    let finalized_count = json
        .get("finalized")
        .and_then(|v| v.as_object())
        .map(|m| m.values().filter(|v| v.as_bool().unwrap_or(false)).count())
        .unwrap_or(0);
    let has_mmd = json
        .get("mmd_method")
        .map(|v| !v.is_null())
        .unwrap_or(false);
    (finalized_count, has_mmd)
}

/// Walk `root_dir` and return a sorted list of patient folders with status badges.
///
/// Status: `complete` if MMD has been run AND ≥1 contour finalized,
/// `in_progress` if any annotation file exists, `not_started` otherwise.
#[tauri::command]
pub async fn list_patients(
    app: tauri::AppHandle,
    root_dir: String,
) -> Result<Vec<PatientInfo>, String> {
    let root = PathBuf::from(&root_dir);
    if !root.is_dir() {
        return Err(format!("not a directory: {root_dir}"));
    }

    // Collect candidate patient directories on a blocking thread (SMB walks
    // can be slow).
    let root_for_walk = root.clone();
    let entries = tokio::task::spawn_blocking(move || -> Result<Vec<(String, PathBuf)>, String> {
        let mut out = Vec::new();
        let read = std::fs::read_dir(&root_for_walk).map_err(|e| format!("read_dir: {e}"))?;
        for entry in read.flatten() {
            if looks_like_patient_dir(&entry) {
                let name = entry.file_name().to_string_lossy().to_string();
                out.push((name, entry.path()));
            }
        }
        out.sort_by(|a, b| a.0.cmp(&b.0));
        Ok(out)
    })
    .await
    .map_err(|e| format!("walk task failed: {e}"))??;

    // Cross-reference annotation + seeds saves on the main task (cheap local
    // FS reads). Status rules, in priority order:
    //   complete    — MMD has been run (mmd_method present in saved JSON)
    //   in_progress — any activity: seeds saved, contours saved, or both
    //   not_started — nothing on disk for this patient yet
    let mut patients = Vec::with_capacity(entries.len());
    for (id, path) in entries {
        let path_str = path.to_string_lossy().to_string();
        let (finalized_count, has_mmd) = read_annotation_summary(&app, &id);
        let has_seeds = patient_has_seeds(&app, &path_str);
        let status = if has_mmd {
            "complete"
        } else if has_seeds || finalized_count > 0 {
            "in_progress"
        } else {
            "not_started"
        };
        patients.push(PatientInfo {
            id,
            path: path_str,
            status: status.to_string(),
            finalized_count,
            has_mmd,
        });
    }

    Ok(patients)
}

/// One immediate subdirectory of a patient folder — typically a single DICOM
/// series (e.g. `MonoPlus_70keV`). Used by the patient browser to expand a
/// patient into its series without reading any DICOM headers.
#[derive(serde::Serialize)]
pub struct SeriesDirInfo {
    /// Folder name (e.g. `MonoPlus_70keV`).
    pub name: String,
    /// Absolute path to the series folder.
    pub path: String,
    /// Number of regular files in the folder (≈ DICOM slices).
    pub num_files: usize,
}

/// List immediate subdirectories of a patient folder with file counts.
///
/// Returns folders sorted by name. Skips hidden / underscore-prefixed entries.
/// Fast — does not parse any DICOM headers.
#[tauri::command]
pub async fn list_series_dirs(patient_path: String) -> Result<Vec<SeriesDirInfo>, String> {
    let root = PathBuf::from(&patient_path);
    if !root.is_dir() {
        return Err(format!("not a directory: {patient_path}"));
    }

    tokio::task::spawn_blocking(move || -> Result<Vec<SeriesDirInfo>, String> {
        let mut out = Vec::new();
        let read = std::fs::read_dir(&root).map_err(|e| format!("read_dir: {e}"))?;
        for entry in read.flatten() {
            let Ok(file_type) = entry.file_type() else { continue };
            if !file_type.is_dir() {
                continue;
            }
            let name = entry.file_name().to_string_lossy().to_string();
            if name.starts_with('.') || name.starts_with('_') {
                continue;
            }
            // Count regular files (cheap directory scan, no DICOM parse).
            let num_files = std::fs::read_dir(entry.path())
                .map(|d| {
                    d.flatten()
                        .filter(|e| e.file_type().map(|t| t.is_file()).unwrap_or(false))
                        .count()
                })
                .unwrap_or(0);
            out.push(SeriesDirInfo {
                name,
                path: entry.path().to_string_lossy().to_string(),
                num_files,
            });
        }
        out.sort_by(|a, b| a.name.cmp(&b.name));
        Ok(out)
    })
    .await
    .map_err(|e| format!("list task failed: {e}"))?
}

// ---------------------------------------------------------------------------
// Fast header-only series scan
// ---------------------------------------------------------------------------

#[derive(Clone, serde::Serialize)]
pub struct ProgressEvent {
    pub phase: &'static str,
    pub done: usize,
    pub total: usize,
    /// Optional detail string (e.g. the current series folder name during a
    /// patient-wide load). `None` for phases that don't carry one.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub detail: Option<String>,
}

/// Scan a DICOM folder for series (header-only; no pixel data decoded).
#[tauri::command]
pub async fn scan_series(
    path: String,
    app: AppHandle,
) -> Result<Vec<SeriesDescriptorDto>, String> {
    let _ = app.emit("dicom_load_progress", ProgressEvent { phase: "scanning", done: 0, total: 0, detail: None });
    let dir = PathBuf::from(path);
    let series = dicom_scan::scan_series(&dir)
        .await
        .map_err(|e| e.to_string())?;
    let _ = app.emit("dicom_load_progress", ProgressEvent { phase: "scanned", done: series.len(), total: series.len(), detail: None });
    Ok(series.into_iter().map(SeriesDescriptorDto::from).collect())
}

#[derive(serde::Serialize)]
pub struct SeriesDescriptorDto {
    pub uid: String,
    pub description: String,
    pub image_comments: Option<String>,
    pub rows: u32,
    pub cols: u32,
    pub num_slices: usize,
    pub pixel_spacing: [f64; 2],
    pub slice_spacing: f64,
    pub orientation: [f64; 6],
    pub rescale_slope: f64,
    pub rescale_intercept: f64,
    pub window_center: f64,
    pub window_width: f64,
    pub patient_name: String,
    pub study_description: String,
    /// Absolute file paths in z-sorted order.
    pub file_paths: Vec<String>,
    pub slice_positions_z: Vec<f64>,
    pub image_position_patient: [f64; 3],
}

impl From<SeriesDescriptor> for SeriesDescriptorDto {
    fn from(d: SeriesDescriptor) -> Self {
        Self {
            uid: d.uid,
            description: d.description,
            image_comments: d.image_comments,
            rows: d.rows,
            cols: d.cols,
            num_slices: d.num_slices,
            pixel_spacing: d.pixel_spacing,
            slice_spacing: d.slice_spacing,
            orientation: d.orientation,
            rescale_slope: d.rescale_slope,
            rescale_intercept: d.rescale_intercept,
            window_center: d.window_center,
            window_width: d.window_width,
            patient_name: d.patient_name,
            study_description: d.study_description,
            file_paths: d.file_paths.into_iter().map(|p| p.to_string_lossy().into_owned()).collect(),
            slice_positions_z: d.slice_positions_z,
            image_position_patient: d.image_position_patient,
        }
    }
}

// ---------------------------------------------------------------------------
// Bulk binary load
// ---------------------------------------------------------------------------

/// Load a single series as one framed binary response:
///   [u32 LE: metadata_json_length] [metadata_json] [i16 LE voxel bytes]
/// Frontend receives this as an ArrayBuffer.
///
/// Also emits `dicom_load_progress` Tauri events during decode so the frontend
/// can show a progress bar, and populates `AppState.volume` so legacy commands
/// (CPR, annotation, MMD) continue to work during gradual migration.
#[tauri::command]
pub async fn load_series(
    dir: String,
    uid: String,
    app: AppHandle,
    state: State<'_, Mutex<AppState>>,
) -> Result<Response, String> {
    push_recent(&app, &dir);
    let dir_path = PathBuf::from(dir);
    let cache_key = (dir_path.to_string_lossy().into_owned(), uid.clone());

    // Rust-side LRU cache lookup. Cloning the Arc is a refcount bump; the
    // underlying voxel buffer is not copied.
    let cached_hit = {
        let mut guard = state.lock().map_err(|e| format!("state lock poisoned: {e}"))?;
        guard.volume_cache.get(&cache_key)
    };

    if let Some(cached) = cached_hit {
        // Fast path: reuse the decoded volume. Swap it into `state.volume` so
        // CPR / FAI / MMD operate on this patient's data after the reload.
        {
            let mut guard = state.lock().map_err(|e| format!("state lock poisoned: {e}"))?;
            guard.volume = Some(cached.volume.clone());
            guard.current_volume_key = Some(cache_key.clone());
            guard.last_metadata = Some(cached.metadata.clone());
        }

        let _ = app.emit(
            "dicom_load_progress",
            ProgressEvent {
                phase: "done",
                done: cached.metadata.num_slices,
                total: cached.metadata.num_slices,
                detail: None,
            },
        );

        // bytemuck::cast_slice is a zero-cost reinterpret over the borrow;
        // .to_vec() is still required for Response::new. The expensive thing
        // we skipped is the parallel DICOM pixel decode.
        let voxel_bytes: Vec<u8> = bytemuck::cast_slice(&cached.voxels_i16[..]).to_vec();
        let framed = encode_frame(&cached.metadata, &voxel_bytes)?;
        return Ok(Response::new(framed));
    }

    // Miss path: do the full decode.
    let app_cb = app.clone();
    let progress: Box<dyn Fn(usize, usize) + Send + Sync> = Box::new(move |done, total| {
        let _ = app_cb.emit(
            "dicom_load_progress",
            ProgressEvent { phase: "decoding", done, total, detail: None },
        );
    });

    let vol = dicom_load::load_series(&dir_path, &uid, Some(progress))
        .await
        .map_err(|e| e.to_string())?;

    // Mirror into legacy AppState so CPR / annotation / MMD keep working.
    bridge_into_state(&vol, &state)?;

    // Break `vol` into pieces now so we can share the voxel buffer via Arc
    // between the cache entry and the framed IPC response without a second
    // ~150 MB memcpy. The `Arc::clone` below is a refcount bump.
    let metadata = vol.metadata;
    let voxels_arc: Arc<Vec<i16>> = Arc::new(vol.voxels_i16);

    // Record identity, metadata, and insert into the LRU cache. Reads the
    // freshly-built `state.volume` (set by bridge_into_state) so the cache
    // stores the exact same f32 Array3 consumers use.
    {
        let mut guard = state.lock().map_err(|e| format!("state lock poisoned: {e}"))?;
        guard.current_volume_key = Some(cache_key.clone());
        guard.last_metadata = Some(metadata.clone());

        let cached_volume = guard
            .volume
            .clone()
            .expect("bridge_into_state populated state.volume");
        guard.volume_cache.insert(
            cache_key,
            CachedVolume {
                metadata: metadata.clone(),
                voxels_i16: Arc::clone(&voxels_arc),
                volume: cached_volume,
            },
        );
    }

    // Emit terminal event so frontend can switch to "finalizing" state.
    let _ = app.emit("dicom_load_progress", ProgressEvent {
        phase: "done",
        done: metadata.num_slices,
        total: metadata.num_slices,
        detail: None,
    });

    let voxel_bytes: Vec<u8> = bytemuck::cast_slice(&voxels_arc[..]).to_vec();
    let framed = encode_frame(&metadata, &voxel_bytes)?;
    Ok(Response::new(framed))
}

/// Query whether the Rust-side volume cache can serve the (dir, uid) request
/// without a decode. On hit, swaps the cached `LoadedVolume` into
/// `state.volume` (so downstream CPR / FAI / MMD calls operate on the
/// right patient) and returns its `VolumeMetadata`. On miss returns `None`
/// and the caller should fall back to `load_series`.
#[tauri::command]
pub async fn reuse_loaded_volume(
    dir: String,
    uid: String,
    state: State<'_, Mutex<AppState>>,
) -> Result<Option<PipelineVolumeMetadata>, String> {
    let key = (dir, uid);
    let mut guard = state.lock().map_err(|e| format!("state lock poisoned: {e}"))?;
    let Some(cached) = guard.volume_cache.get(&key) else {
        return Ok(None);
    };
    let metadata = cached.metadata.clone();
    guard.volume = Some(cached.volume.clone());
    guard.current_volume_key = Some(key);
    guard.last_metadata = Some(metadata.clone());
    Ok(Some(metadata))
}

// ---------------------------------------------------------------------------
// Dual-energy fast-path
// ---------------------------------------------------------------------------

/// Load two DICOM series in parallel (one per energy) and populate the
/// dual-energy state slot. `low_dir` / `high_dir` are folder paths whose
/// names must contain a keV label, e.g. `MonoPlus_70keV` — lab-internal
/// data has `ImageComments` stripped and `SeriesDescription` mislabeled,
/// so folder name is the only reliable source.
///
/// Also mirrors the low-energy volume into `state.volume` so CPR /
/// annotation / FAI continue to operate on the same frame of reference.
///
/// Returns the framed binary response of the LOW-energy volume so the
/// frontend can build its primary cornerstone3D volume exactly as with
/// `load_series`; the high-energy voxels stay resident only in
/// `state.dual_energy` for MMD.
#[tauri::command]
pub async fn load_dual_energy(
    low_dir: String,
    high_dir: String,
    app: AppHandle,
    state: State<'_, Mutex<AppState>>,
) -> Result<Response, String> {
    let low_kev = parse_kev_from_folder(&low_dir)
        .ok_or_else(|| format!(
            "cannot extract keV from low-energy folder '{}'. \
             Rename to include 'NNkeV' (e.g. MonoPlus_70keV).",
            folder_name(&low_dir),
        ))?;
    let high_kev = parse_kev_from_folder(&high_dir)
        .ok_or_else(|| format!(
            "cannot extract keV from high-energy folder '{}'. \
             Rename to include 'NNkeV' (e.g. MonoPlus_150keV).",
            folder_name(&high_dir),
        ))?;
    if (low_kev - high_kev).abs() < 1.0 {
        return Err(format!("low and high energies are both {low_kev} keV"));
    }

    push_recent(&app, &low_dir);

    let low_path = PathBuf::from(&low_dir);
    let high_path = PathBuf::from(&high_dir);

    // Progress callback: weight low 0..50%, high 50..100%.
    let app_cb_low = app.clone();
    let low_progress: Box<dyn Fn(usize, usize) + Send + Sync> = Box::new(move |done, total| {
        let scaled = done.saturating_mul(50) / total.max(1);
        let _ = app_cb_low.emit(
            "dicom_load_progress",
            ProgressEvent { phase: "decoding", done: scaled, total: 100, detail: None },
        );
    });
    let app_cb_high = app.clone();
    let high_progress: Box<dyn Fn(usize, usize) + Send + Sync> = Box::new(move |done, total| {
        let scaled = 50 + done.saturating_mul(50) / total.max(1);
        let _ = app_cb_high.emit(
            "dicom_load_progress",
            ProgressEvent { phase: "decoding", done: scaled, total: 100, detail: None },
        );
    });

    // Find the (single) series in each folder, then decode both in parallel.
    let (low_vol, high_vol) = tokio::try_join!(
        load_first_series(low_path.clone(), Some(low_progress)),
        load_first_series(high_path.clone(), Some(high_progress)),
    )?;

    // Volumes must be on the same voxel grid for MMD.
    if (low_vol.metadata.rows, low_vol.metadata.cols, low_vol.metadata.num_slices)
        != (high_vol.metadata.rows, high_vol.metadata.cols, high_vol.metadata.num_slices)
    {
        return Err(format!(
            "low/high volume shape mismatch: \
             low {}×{}×{}, high {}×{}×{}",
            low_vol.metadata.num_slices, low_vol.metadata.rows, low_vol.metadata.cols,
            high_vol.metadata.num_slices, high_vol.metadata.rows, high_vol.metadata.cols,
        ));
    }

    // Mirror low into state.volume via the existing bridge, so CPR / FAI work.
    bridge_into_state(&low_vol, &state)?;

    // Convert both i16 HU → f32 Array3, store as DualEnergyVolume.
    let de = build_dual_energy_volume(&low_vol, &high_vol, low_kev, high_kev)?;
    {
        let mut guard = state.lock().map_err(|e| format!("state lock poisoned: {e}"))?;
        guard.dual_energy = Some(de);
    }

    let _ = app.emit("dicom_load_progress", ProgressEvent {
        phase: "done",
        done: 100,
        total: 100,
        detail: None,
    });

    let voxel_bytes: Vec<u8> = bytemuck::cast_slice(&low_vol.voxels_i16).to_vec();
    let framed = encode_frame(&low_vol.metadata, &voxel_bytes)?;
    Ok(Response::new(framed))
}

/// Scan a folder for DICOM series and load the first one found. Errors if the
/// folder is empty or contains no DICOM images.
async fn load_first_series(
    dir: PathBuf,
    on_progress: Option<Box<dyn Fn(usize, usize) + Send + Sync>>,
) -> Result<PipelineLoadedVolume, String> {
    let series = dicom_scan::scan_series(&dir)
        .await
        .map_err(|e| format!("scan {}: {e}", dir.display()))?;
    let first = series.into_iter().next().ok_or_else(|| {
        format!("no DICOM series found in {}", dir.display())
    })?;
    dicom_load::load_series(&dir, &first.uid, on_progress)
        .await
        .map_err(|e| format!("load {}: {e}", dir.display()))
}

/// Extract a keV number from a folder path's trailing component.
/// Matches patterns like `MonoPlus_70keV`, `70 keV`, `150keV_Soft`, case-insensitive.
fn parse_kev_from_folder(path: &str) -> Option<f64> {
    use regex::Regex;
    static REGEX: std::sync::OnceLock<Regex> = std::sync::OnceLock::new();
    let re = REGEX.get_or_init(|| {
        Regex::new(r"(?i)(\d+(?:\.\d+)?)\s*keV").unwrap()
    });
    let name = folder_name(path);
    re.captures(&name)?.get(1)?.as_str().parse::<f64>().ok()
}

fn folder_name(path: &str) -> String {
    Path::new(path)
        .file_name()
        .map(|f| f.to_string_lossy().into_owned())
        .unwrap_or_else(|| path.to_string())
}

fn build_dual_energy_volume(
    low: &PipelineLoadedVolume,
    high: &PipelineLoadedVolume,
    low_kev: f64,
    high_kev: f64,
) -> Result<pcat_pipeline::dicom_loader::DualEnergyVolume, String> {
    let m = &low.metadata;
    let nz = m.num_slices;
    let ny = m.rows as usize;
    let nx = m.cols as usize;
    let shape = (nz, ny, nx);

    let low_f32: Vec<f32> = low.voxels_i16.iter().map(|&v| v as f32).collect();
    let high_f32: Vec<f32> = high.voxels_i16.iter().map(|&v| v as f32).collect();
    let low_arr = Array3::from_shape_vec(shape, low_f32)
        .map_err(|e| format!("low shape: {e}"))?;
    let high_arr = Array3::from_shape_vec(shape, high_f32)
        .map_err(|e| format!("high shape: {e}"))?;

    let iop = m.orientation;
    let row = [iop[0], iop[1], iop[2]];
    let col = [iop[3], iop[4], iop[5]];
    let normal = [
        row[1] * col[2] - row[2] * col[1],
        row[2] * col[0] - row[0] * col[2],
        row[0] * col[1] - row[1] * col[0],
    ];
    let direction = [
        row[0], row[1], row[2],
        col[0], col[1], col[2],
        normal[0], normal[1], normal[2],
    ];
    let spacing = [m.slice_spacing, m.pixel_spacing[0], m.pixel_spacing[1]];
    // ZYX order; the LPS x/y components were silently dropped before, which
    // miscomputed voxel indices for any acquisition not centered at isocenter.
    let ipp = m.image_position_patient;
    let origin = [
        m.slice_positions_z.first().copied().unwrap_or(ipp[2]),
        ipp[1],
        ipp[0],
    ];

    Ok(pcat_pipeline::dicom_loader::DualEnergyVolume {
        low: Arc::new(low_arr),
        high: Arc::new(high_arr),
        low_energy_kev: low_kev,
        high_energy_kev: high_kev,
        spacing,
        origin,
        direction,
        patient_name: m.patient_name.clone(),
        study_description: m.study_description.clone(),
    })
}

fn bridge_into_state(
    vol: &PipelineLoadedVolume,
    state: &State<'_, Mutex<AppState>>,
) -> Result<(), String> {
    let meta = &vol.metadata;
    let nz = meta.num_slices;
    let ny = meta.rows as usize;
    let nx = meta.cols as usize;

    // Convert i16 HU → f32 for the legacy Array3<f32> consumers.
    let data_f32: Vec<f32> = vol.voxels_i16.iter().map(|&v| v as f32).collect();
    let arr = Array3::from_shape_vec((nz, ny, nx), data_f32)
        .map_err(|e| format!("volume shape mismatch: {e}"))?;

    // Direction: row-major 3x3. Use IOP row × IOP col × normal.
    let iop = meta.orientation;
    let iop_row = [iop[0], iop[1], iop[2]];
    let iop_col = [iop[3], iop[4], iop[5]];
    let normal = [
        iop_row[1] * iop_col[2] - iop_row[2] * iop_col[1],
        iop_row[2] * iop_col[0] - iop_row[0] * iop_col[2],
        iop_row[0] * iop_col[1] - iop_row[1] * iop_col[0],
    ];
    let direction = [
        iop_row[0], iop_row[1], iop_row[2],
        iop_col[0], iop_col[1], iop_col[2],
        normal[0], normal[1], normal[2],
    ];

    let spacing = [meta.slice_spacing, meta.pixel_spacing[0], meta.pixel_spacing[1]];
    // ZYX order; the LPS x/y components were silently dropped before. Sampling
    // code (CPR, ROI, radial-angular) reads `origin[1]` and `origin[2]` and
    // would shift voxel indices by the unrecorded patient x/y offset.
    let ipp = meta.image_position_patient;
    let origin = [
        meta.slice_positions_z.first().copied().unwrap_or(ipp[2]),
        ipp[1],
        ipp[0],
    ];

    let legacy = StateLoadedVolume {
        data: Arc::new(arr),
        spacing,
        origin,
        direction,
        window_center: meta.window_center,
        window_width: meta.window_width,
        patient_name: meta.patient_name.clone(),
        study_description: meta.study_description.clone(),
    };

    let mut guard = state.lock().map_err(|e| format!("state lock poisoned: {e}"))?;
    guard.volume = Some(legacy);
    Ok(())
}

// ---------------------------------------------------------------------------
// Patient-level loader + active-volume switcher
// ---------------------------------------------------------------------------

/// One series in a patient load, surfaced to the frontend switcher.
#[derive(serde::Serialize)]
pub struct LoadedSeriesDescriptor {
    /// Folder name, e.g. `MonoPlus_70keV` or `CCTA_Soft`.
    pub name: String,
    /// Absolute path to the series folder (used as cache key + switch key).
    pub path: String,
    /// DICOM SeriesInstanceUID — second half of the cache key.
    pub uid: String,
    /// SeriesDescription from the DICOM header (may be mislabeled for MonoPlus).
    pub series_description: String,
    /// keV parsed from the folder name, if present. `None` for CaScore/CCTA.
    pub kev: Option<f64>,
    /// Number of slices (metadata-only, no full decode).
    pub num_slices: usize,
    /// Shape [rows, cols] — useful for grouping volumes on the same grid.
    pub rows: usize,
    pub cols: usize,
}

#[derive(serde::Serialize)]
pub struct PatientLoadResult {
    pub series: Vec<LoadedSeriesDescriptor>,
    /// Index into `series` for the volume now in `state.volume`.
    pub active_index: usize,
    /// Errors hit while loading individual series — the command does not
    /// abort the whole patient on a single decode failure. Entries are
    /// `"<folder name>: <error>"`.
    pub failures: Vec<String>,
}

/// Load every DICOM series under `patient_dir` into the volume cache so the
/// user can switch between modalities (CaScore, CCTA, multiple MonoPlus keV)
/// without redecoding. Picks a default "active" volume — the CCTA if
/// recognizable, else the lowest keV MonoPlus, else the first series.
///
/// Also auto-pairs the first two MonoPlus keV series as the MMD dual-energy
/// volume (replaces `state.dual_energy`).
#[tauri::command]
pub async fn load_patient_all(
    patient_dir: String,
    app: AppHandle,
    state: State<'_, Mutex<AppState>>,
) -> Result<PatientLoadResult, String> {
    let root = PathBuf::from(&patient_dir);
    if !root.is_dir() {
        return Err(format!("not a directory: {patient_dir}"));
    }

    // Discover series subfolders (cheap, no DICOM parse).
    let subdirs = tokio::task::spawn_blocking({
        let root = root.clone();
        move || -> Result<Vec<(String, PathBuf)>, String> {
            let mut out = Vec::new();
            let read = std::fs::read_dir(&root).map_err(|e| format!("read_dir: {e}"))?;
            for entry in read.flatten() {
                let Ok(ty) = entry.file_type() else { continue };
                if !ty.is_dir() { continue; }
                let name = entry.file_name().to_string_lossy().to_string();
                if name.starts_with('.') || name.starts_with('_') { continue; }
                out.push((name, entry.path()));
            }
            out.sort_by(|a, b| a.0.cmp(&b.0));
            Ok(out)
        }
    })
    .await
    .map_err(|e| format!("subdir task failed: {e}"))??;

    if subdirs.is_empty() {
        return Err(format!("no series subfolders under {patient_dir}"));
    }

    push_recent(&app, &patient_dir);

    let total = subdirs.len();
    let mut descriptors: Vec<LoadedSeriesDescriptor> = Vec::new();
    let mut failures: Vec<String> = Vec::new();

    for (i, (name, series_dir)) in subdirs.into_iter().enumerate() {
        // Phase "patient_series" signals the start of a new series and carries
        // its folder name as `detail` so the footer can show "Loading X/Y: name".
        let _ = app.emit(
            "dicom_load_progress",
            ProgressEvent {
                phase: "patient_series",
                done: i,
                total,
                detail: Some(name.clone()),
            },
        );

        // Scan to get the series UID (first series in the folder).
        let scan = match dicom_scan::scan_series(&series_dir).await {
            Ok(s) => s,
            Err(e) => {
                failures.push(format!("{name}: scan failed: {e}"));
                continue;
            }
        };
        let Some(first) = scan.into_iter().next() else {
            failures.push(format!("{name}: no DICOM series found"));
            continue;
        };
        let uid = first.uid.clone();
        let path_str = series_dir.to_string_lossy().into_owned();
        let cache_key = (path_str.clone(), uid.clone());

        // Skip decode if already cached.
        let already_cached = {
            let mut guard = state.lock().map_err(|e| format!("state lock poisoned: {e}"))?;
            guard.volume_cache.get(&cache_key).is_some()
        };

        let metadata: PipelineVolumeMetadata = if already_cached {
            // Already resident — just read the cached metadata back.
            let mut guard = state.lock().map_err(|e| format!("state lock poisoned: {e}"))?;
            guard.volume_cache.get(&cache_key).unwrap().metadata.clone()
        } else {
            // Fine-grained per-slice progress so the footer bar animates
            // smoothly within each series. The coarse "patient_series" event
            // emitted just before this series started already carries the
            // series index + name, so we only need `decoding` here.
            let app_cb = app.clone();
            let progress: Box<dyn Fn(usize, usize) + Send + Sync> = Box::new(move |done, total_slices| {
                let _ = app_cb.emit(
                    "dicom_load_progress",
                    ProgressEvent { phase: "decoding", done, total: total_slices, detail: None },
                );
            });

            let vol = match dicom_load::load_series(&series_dir, &uid, Some(progress)).await {
                Ok(v) => v,
                Err(e) => {
                    failures.push(format!("{name}: decode failed: {e}"));
                    continue;
                }
            };

            if let Err(e) = bridge_into_state(&vol, &state) {
                failures.push(format!("{name}: bridge failed: {e}"));
                continue;
            }

            let meta_clone = vol.metadata.clone();
            let voxels_arc: Arc<Vec<i16>> = Arc::new(vol.voxels_i16);
            {
                let mut guard = state.lock().map_err(|e| format!("state lock poisoned: {e}"))?;
                guard.current_volume_key = Some(cache_key.clone());
                guard.last_metadata = Some(meta_clone.clone());
                let cached_volume = guard
                    .volume
                    .clone()
                    .expect("bridge_into_state populated state.volume");
                guard.volume_cache.insert(
                    cache_key.clone(),
                    CachedVolume {
                        metadata: meta_clone.clone(),
                        voxels_i16: Arc::clone(&voxels_arc),
                        volume: cached_volume,
                    },
                );
            }
            meta_clone
        };

        descriptors.push(LoadedSeriesDescriptor {
            name: name.clone(),
            path: path_str,
            uid,
            series_description: metadata.series_description.clone(),
            kev: parse_kev_from_folder(&name),
            num_slices: metadata.num_slices,
            rows: metadata.rows as usize,
            cols: metadata.cols as usize,
        });
    }

    if descriptors.is_empty() {
        return Err(format!(
            "failed to load any series from {patient_dir}: {:?}",
            failures
        ));
    }

    // Pick the default active volume: CCTA-like name first, else lowest keV,
    // else first in sort order.
    let active_index = descriptors
        .iter()
        .position(|d| {
            let n = d.name.to_ascii_lowercase();
            n.contains("ccta")
        })
        .or_else(|| {
            // Lowest-keV MonoPlus.
            let mut best: Option<(usize, f64)> = None;
            for (i, d) in descriptors.iter().enumerate() {
                if let Some(k) = d.kev {
                    if best.map(|(_, bk)| k < bk).unwrap_or(true) {
                        best = Some((i, k));
                    }
                }
            }
            best.map(|(i, _)| i)
        })
        .unwrap_or(0);

    // Bridge the active one into state.volume (may already be there if it
    // was the last-loaded series; the cache-get-then-write is cheap).
    {
        let active = &descriptors[active_index];
        let key = (active.path.clone(), active.uid.clone());
        let mut guard = state.lock().map_err(|e| format!("state lock poisoned: {e}"))?;
        if let Some(cached) = guard.volume_cache.get(&key) {
            guard.volume = Some(cached.volume.clone());
            guard.current_volume_key = Some(key);
            guard.last_metadata = Some(cached.metadata.clone());
        }
    }

    // Auto-pair the two lowest-keV MonoPlus series into state.dual_energy
    // so MMD can run without a separate dual-energy load step. If the user
    // loaded a patient without a keV pair, leave dual_energy alone.
    let mut kev_entries: Vec<(usize, f64)> = descriptors
        .iter()
        .enumerate()
        .filter_map(|(i, d)| d.kev.map(|k| (i, k)))
        .collect();
    kev_entries.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
    if kev_entries.len() >= 2 {
        let (low_idx, low_kev) = kev_entries[0];
        let (high_idx, high_kev) = kev_entries[kev_entries.len() - 1];
        let low_key = (
            descriptors[low_idx].path.clone(),
            descriptors[low_idx].uid.clone(),
        );
        let high_key = (
            descriptors[high_idx].path.clone(),
            descriptors[high_idx].uid.clone(),
        );

        // Pull both voxel buffers out from cache and rebuild the f32
        // Array3 pair for DualEnergyVolume.
        let (low_meta, low_voxels, high_meta, high_voxels) = {
            let mut guard = state.lock().map_err(|e| format!("state lock poisoned: {e}"))?;
            let low = guard.volume_cache.get(&low_key);
            let high = guard.volume_cache.get(&high_key);
            match (low, high) {
                (Some(l), Some(h)) => (
                    l.metadata.clone(),
                    Arc::clone(&l.voxels_i16),
                    h.metadata.clone(),
                    Arc::clone(&h.voxels_i16),
                ),
                _ => {
                    // One missing — skip DE pairing.
                    drop(guard);
                    let _ = app.emit(
                        "dicom_load_progress",
                        ProgressEvent { phase: "done", done: total, total, detail: None },
                    );
                    return Ok(PatientLoadResult { series: descriptors, active_index, failures });
                }
            }
        };

        let low_vol = PipelineLoadedVolume { metadata: low_meta, voxels_i16: (*low_voxels).clone() };
        let high_vol = PipelineLoadedVolume { metadata: high_meta, voxels_i16: (*high_voxels).clone() };

        if low_vol.metadata.rows == high_vol.metadata.rows
            && low_vol.metadata.cols == high_vol.metadata.cols
            && low_vol.metadata.num_slices == high_vol.metadata.num_slices
        {
            match build_dual_energy_volume(&low_vol, &high_vol, low_kev, high_kev) {
                Ok(de) => {
                    let mut guard = state.lock().map_err(|e| format!("state lock poisoned: {e}"))?;
                    guard.dual_energy = Some(de);
                }
                Err(e) => {
                    failures.push(format!("dual-energy pairing failed: {e}"));
                }
            }
        } else {
            failures.push(
                "dual-energy pairing skipped: keV series have different voxel grids".into(),
            );
        }
    }

    let _ = app.emit(
        "dicom_load_progress",
        ProgressEvent { phase: "done", done: total, total, detail: None },
    );

    Ok(PatientLoadResult { series: descriptors, active_index, failures })
}

/// Switch the active volume to a previously-loaded series (must be resident
/// in the volume cache — call `load_patient_all` or `load_series` first).
///
/// Returns the framed low-energy-style binary bundle so the frontend can
/// rebuild the cornerstone3D volume identical to `load_series`.
#[tauri::command]
pub async fn set_active_volume(
    dir: String,
    uid: String,
    state: State<'_, Mutex<AppState>>,
) -> Result<Response, String> {
    let key = (dir.clone(), uid.clone());

    // Fast path: the series is already decoded + cached.
    let cached_hit = {
        let mut guard = state.lock().map_err(|e| format!("state lock poisoned: {e}"))?;
        if let Some(cached) = guard.volume_cache.get(&key) {
            guard.volume = Some(cached.volume.clone());
            guard.current_volume_key = Some(key.clone());
            guard.last_metadata = Some(cached.metadata.clone());
            Some((cached.metadata.clone(), Arc::clone(&cached.voxels_i16)))
        } else {
            None
        }
    };

    let (metadata, voxels) = match cached_hit {
        Some(hit) => hit,
        None => {
            // Lazy decode-on-demand: `load_patient_all` scans every series but
            // only decodes the ones needed up front (active + dual-energy pair).
            // Switching to any other series decodes it here on first access,
            // then it stays cached. This is what keeps the initial load fast.
            let vol = dicom_load::load_series(Path::new(&dir), &uid, None)
                .await
                .map_err(|e| format!("decode failed: {e}"))?;
            bridge_into_state(&vol, &state)?;
            let meta = vol.metadata.clone();
            let voxels_arc: Arc<Vec<i16>> = Arc::new(vol.voxels_i16);
            let mut guard = state.lock().map_err(|e| format!("state lock poisoned: {e}"))?;
            guard.current_volume_key = Some(key.clone());
            guard.last_metadata = Some(meta.clone());
            let cached_volume = guard
                .volume
                .clone()
                .expect("bridge_into_state populated state.volume");
            guard.volume_cache.insert(
                key.clone(),
                CachedVolume {
                    metadata: meta.clone(),
                    voxels_i16: Arc::clone(&voxels_arc),
                    volume: cached_volume,
                },
            );
            (meta, voxels_arc)
        }
    };

    let voxel_bytes: Vec<u8> = bytemuck::cast_slice(&voxels[..]).to_vec();
    let framed = encode_frame(&metadata, &voxel_bytes)?;
    Ok(Response::new(framed))
}
