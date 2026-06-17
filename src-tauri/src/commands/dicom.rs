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
use crate::commands::flag::{flag_file_path, read_flag_at};

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

/// Panic-safe stderr line. Plain `eprintln!` calls `.expect("failed printing
/// to stderr")` internally and ABORTS the whole process (SIGABRT) if the write
/// fails — e.g. a broken output pipe, or a bundled `.app` launched with no
/// connected console. We hit exactly that crash. This ignores the write error.
macro_rules! log_line {
    ($($arg:tt)*) => {{
        use std::io::Write as _;
        let _ = writeln!(std::io::stderr(), $($arg)*);
    }};
}

/// Sanitize a path string into a safe filename component.
fn sanitize_for_filename(s: &str) -> String {
    s.chars()
        .map(|c| if c.is_ascii_alphanumeric() || c == '-' || c == '_' || c == '.' { c } else { '_' })
        .collect::<String>()
        .replace("..", "_")
}

/// Per-patient file key derived from the **full** DICOM folder path, so two
/// patients whose folders share a last component (e.g. every patient's
/// `MonoPlus_70keV`) never collide. The last component is prefixed for human
/// browseability; the sanitized full path after it guarantees uniqueness.
/// SINGLE source of this key for seeds, sessions, AND annotations — route all
/// per-patient persistence through here so the schemes cannot drift apart.
pub(crate) fn patient_file_key(dicom_path: &str) -> String {
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
    let path = dir.join(patient_file_key(&dicom_path));
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
    let path = dir.join(patient_file_key(&dicom_path));
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
    /// Always 0. The annotations-based finalized-contour count is no longer
    /// derived; this field is kept only for API stability.
    pub finalized_count: usize,
    /// Whether MMD has been run, derived from the saved session bundle
    /// (`sessions/{key}.json`): true when `mmd.summary` is non-null.
    pub has_mmd: bool,
    /// Whether the patient's data has been flagged as problematic / unusable.
    pub flagged: bool,
    /// Optional note explaining the flag (None when not flagged).
    pub flag_note: Option<String>,
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

/// Whether a parsed session bundle contains FAI results and an MMD summary.
/// Contract with src/lib/session.ts `saveSession`: the bundle has
/// `fai: object|null` (a non-empty object means FAI ran) and
/// `mmd: { summary: any|null }` (non-null summary means MMD ran). These three
/// field names are the contract between session.ts and this reader; the tests
/// above pin it.
fn session_flags(bundle: &serde_json::Value) -> (bool, bool) {
    let has_fai = bundle
        .get("fai")
        .and_then(|v| v.as_object())
        .map(|o| !o.is_empty())
        .unwrap_or(false);
    let has_mmd = bundle
        .get("mmd")
        .and_then(|m| m.get("summary"))
        .map(|s| !s.is_null())
        .unwrap_or(false);
    (has_fai, has_mmd)
}

/// Derived patient status label. Mirrored (kept trivially simple) in
/// src/lib/patientStatus.ts `derivePatientStatus` for the live footer; this Rust
/// version is the canonical, tested one driving the patient-list badges.
fn status_for(complete: bool, has_seeds: bool) -> &'static str {
    if complete {
        "complete"
    } else if has_seeds {
        "in_progress"
    } else {
        "not_started"
    }
}

/// Scan saved session bundles for this patient. Sessions are keyed per *series*
/// (`patient_file_key` of a series subfolder) and a patient folder holds several
/// series, so we substring-match the sanitized patient path exactly like
/// `patient_has_seeds`. Returns (complete, has_mmd): `complete` iff ANY single
/// session has both FAI and MMD (that series is fully analyzed).
fn patient_session_summary(app: &tauri::AppHandle, patient_path: &str) -> (bool, bool) {
    let dir = app
        .path()
        .app_data_dir()
        .expect("app data dir")
        .join("sessions");
    let needle = sanitize_for_filename(patient_path);
    if needle.is_empty() {
        return (false, false);
    }
    let Ok(read) = std::fs::read_dir(&dir) else {
        return (false, false);
    };
    let mut complete = false;
    let mut has_mmd = false;
    for entry in read.flatten() {
        let name = entry.file_name().to_string_lossy().to_string();
        if !name.contains(&needle) {
            continue;
        }
        let Ok(data) = std::fs::read_to_string(entry.path()) else {
            continue;
        };
        let Ok(json) = serde_json::from_str::<serde_json::Value>(&data) else {
            continue;
        };
        let (f, m) = session_flags(&json);
        complete |= f && m;
        has_mmd |= m;
    }
    (complete, has_mmd)
}

/// Walk `root_dir` and return a sorted list of patient folders with status badges.
///
/// Status is derived from the per-series session bundle (`sessions/{key}.json`):
/// - `complete`     — any session has BOTH FAI results AND an MMD summary in the same bundle.
/// - `in_progress`  — the patient has seeds saved (but no complete session yet).
/// - `not_started`  — neither condition holds.
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

    // Derive status from session bundles + seeds (cheap local FS reads).
    // Session bundles live at sessions/{key}.json, one per series. Status rules:
    //   complete    — any session bundle has BOTH FAI results and MMD summary
    //   in_progress — seeds exist for this patient (no complete session yet)
    //   not_started — nothing on disk for this patient
    let flags_base = app
        .path()
        .app_data_dir()
        .expect("app data dir")
        .join("flags");
    let mut patients = Vec::with_capacity(entries.len());
    for (id, path) in entries {
        let path_str = path.to_string_lossy().to_string();
        let has_seeds = patient_has_seeds(&app, &path_str);
        let (complete, has_mmd) = patient_session_summary(&app, &path_str);
        let status = status_for(complete, has_seeds);
        let flag = read_flag_at(&flag_file_path(&flags_base, &path_str))
            .ok()
            .flatten();
        patients.push(PatientInfo {
            id,
            path: path_str,
            status: status.to_string(),
            // Finalized contour count is not persisted in the session bundle;
            // the old annotations-based count was already always 0. Kept for API
            // stability.
            finalized_count: 0,
            has_mmd,
            flagged: flag.as_ref().map(|f| f.flagged).unwrap_or(false),
            flag_note: flag.and_then(|f| if f.note.is_empty() { None } else { Some(f.note) }),
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

    // A single-series open is NOT a dual-energy pair. Drop any stale
    // `dual_energy` + `wl_calibration` from a previous load so water/lipid can't
    // silently run on a different series' data (it needs two energies).
    {
        let mut guard = state.lock().map_err(|e| format!("state lock poisoned: {e}"))?;
        guard.dual_energy = None;
        guard.wl_calibration = None;
    }

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
        guard.wl_calibration = None; // new pair → discard any prior calibration
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
    let hm = &high.metadata;
    let ny = m.rows as usize;
    let nx = m.cols as usize;
    if hm.rows as usize != ny || hm.cols as usize != nx {
        return Err(format!(
            "in-plane grids differ: {ny}×{nx} ({low_kev} keV) vs {}×{} ({high_kev} keV)",
            hm.rows, hm.cols
        ));
    }
    let slice_len = ny * nx;

    // The two VMI reconstructions come from the same scan and share slice
    // z-positions, but a recon can have a slice or two more/fewer at an end
    // (e.g. 429 vs 428). Pair slices by z-position (nearest within tolerance) so
    // the dual-energy volume is co-registered even when num_slices differs,
    // rather than rejecting the pair outright.
    const Z_TOL_MM: f64 = 0.10;
    let lz = &m.slice_positions_z;
    let hz = &hm.slice_positions_z;
    if lz.len() != m.num_slices || hz.len() != hm.num_slices {
        return Err("slice_positions_z length mismatch — cannot align dual-energy".into());
    }
    let mut pairs: Vec<(usize, usize)> = Vec::with_capacity(lz.len().min(hz.len()));
    let mut used_high = vec![false; hz.len()];
    for (li, &z) in lz.iter().enumerate() {
        let mut best: Option<(usize, f64)> = None;
        for (hj, &hzj) in hz.iter().enumerate() {
            if used_high[hj] {
                continue;
            }
            let d = (hzj - z).abs();
            if best.map(|(_, bd)| d < bd).unwrap_or(true) {
                best = Some((hj, d));
            }
        }
        if let Some((hj, d)) = best {
            if d <= Z_TOL_MM {
                used_high[hj] = true;
                pairs.push((li, hj));
            }
        }
    }
    if pairs.len() < 2 {
        return Err(format!(
            "z-overlap too small: {} common slices ({low_kev} keV has {}, {high_kev} keV has {})",
            pairs.len(),
            lz.len(),
            hz.len()
        ));
    }

    let nz = pairs.len();
    let shape = (nz, ny, nx);
    let mut low_f32 = vec![0f32; nz * slice_len];
    let mut high_f32 = vec![0f32; nz * slice_len];
    for (out_k, &(li, hj)) in pairs.iter().enumerate() {
        let lo = &low.voxels_i16[li * slice_len..(li + 1) * slice_len];
        let hi = &high.voxels_i16[hj * slice_len..(hj + 1) * slice_len];
        let dst = out_k * slice_len;
        for t in 0..slice_len {
            low_f32[dst + t] = lo[t] as f32;
            high_f32[dst + t] = hi[t] as f32;
        }
    }
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
    // Z is the first MATCHED slice (the built volume starts there, not at the
    // low series' first slice which may have been trimmed by the z-alignment).
    let ipp = m.image_position_patient;
    let origin = [
        lz.get(pairs[0].0).copied().unwrap_or(ipp[2]),
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

    // Fresh patient — drop the previous patient's dual-energy pair + calibration.
    // If this patient has a keV pair it's re-established below; if not,
    // `dual_energy` correctly stays None (water/lipid will refuse to run).
    {
        let mut guard = state.lock().map_err(|e| format!("state lock poisoned: {e}"))?;
        guard.dual_energy = None;
        guard.wl_calibration = None;
    }

    let total = subdirs.len();
    let mut descriptors: Vec<LoadedSeriesDescriptor> = Vec::new();
    // Parallel to `descriptors` — the folder each series lives in, for lazy decode.
    let mut series_dirs: Vec<PathBuf> = Vec::new();
    let mut failures: Vec<String> = Vec::new();
    let scan_started = std::time::Instant::now();

    // Phase 1 — SCAN every series (header-only, fast, NO pixel decode) to build
    // the volume-switcher list. Decoding all series up front was the slow path
    // (a NAEOTOM patient has many keV series + CCTA + CaScore, each decoded
    // sequentially). We decode only what's needed in Phase 2 and lazy-decode the
    // rest on demand in `set_active_volume`.
    for (i, (name, series_dir)) in subdirs.into_iter().enumerate() {
        let _ = app.emit(
            "dicom_load_progress",
            ProgressEvent {
                phase: "scanning",
                done: i,
                total,
                detail: Some(name.clone()),
            },
        );

        // Quick-scan: one header + a file count, NOT a header read of every
        // slice. The full per-slice scan happens lazily inside the Phase-2
        // decode (only for the series actually decoded), so reading all headers
        // here was redundant — and dominated the load time on SMB.
        let qs = match dicom_scan::quick_scan_series(&series_dir).await {
            Ok(Some(q)) => q,
            Ok(None) => {
                failures.push(format!("{name}: no DICOM series found"));
                continue;
            }
            Err(e) => {
                failures.push(format!("{name}: scan failed: {e}"));
                continue;
            }
        };

        descriptors.push(LoadedSeriesDescriptor {
            name: name.clone(),
            path: series_dir.to_string_lossy().into_owned(),
            uid: qs.uid,
            series_description: qs.description,
            kev: parse_kev_from_folder(&name),
            num_slices: qs.num_slices,
            rows: qs.rows as usize,
            cols: qs.cols as usize,
        });
        series_dirs.push(series_dir);
    }
    log_line!(
        "[load-timing] load_patient_all: scanned {} series (header-only) in {:.2?}",
        descriptors.len(),
        scan_started.elapsed()
    );

    if descriptors.is_empty() {
        return Err(format!(
            "failed to load any series from {patient_dir}: {:?}",
            failures
        ));
    }

    // CCTA-like series (high-res contrast scan), if present — used for the
    // centerline, so it's decoded up front even though it isn't the default view.
    let ccta_index = descriptors.iter().position(|d| {
        let n = d.name.to_ascii_lowercase();
        n.contains("ccta")
    });
    // Lowest-keV MonoPlus index, if any.
    let lowest_kev_index = {
        let mut best: Option<(usize, f64)> = None;
        for (i, d) in descriptors.iter().enumerate() {
            if let Some(k) = d.kev {
                if best.map(|(_, bk)| k < bk).unwrap_or(true) {
                    best = Some((i, k));
                }
            }
        }
        best.map(|(i, _)| i)
    };

    // Default the displayed volume to the lowest-keV series (e.g. 70 keV): the
    // low-keV image has the highest soft-tissue/fat contrast, which is where
    // pericoronary-fat annotation is done. Fall back to CCTA, then first.
    let active_index = lowest_kev_index.or(ccta_index).unwrap_or(0);

    // Phase 2 — decode ONLY the series the clinician needs right away: the
    // displayed one (lowest keV), the dual-energy keV pair (lowest + highest, for
    // MMD / water-lipid), and the CCTA (for the centerline). Everything else
    // (e.g. CA_SCORING) stays scanned-but-not-decoded and is decoded on first
    // switch in `set_active_volume`. This is what keeps the load fast.
    let mut need: Vec<usize> = vec![active_index];
    if let Some(c) = ccta_index {
        need.push(c);
    }
    {
        let mut kevs: Vec<(usize, f64)> = descriptors
            .iter()
            .enumerate()
            .filter_map(|(i, d)| d.kev.map(|k| (i, k)))
            .collect();
        kevs.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
        if kevs.len() >= 2 {
            need.push(kevs[0].0); // lowest keV
            need.push(kevs[kevs.len() - 1].0); // highest keV
        }
    }
    need.sort_unstable();
    need.dedup();

    let decode_started = std::time::Instant::now();
    let n_need = need.len();
    for (k, &idx) in need.iter().enumerate() {
        // Coarse per-series progress for the decode phase (the long leg).
        let _ = app.emit(
            "dicom_load_progress",
            ProgressEvent {
                phase: "patient_series",
                done: k,
                total: n_need,
                detail: Some(descriptors[idx].name.clone()),
            },
        );

        let cache_key = (descriptors[idx].path.clone(), descriptors[idx].uid.clone());
        let already_cached = {
            let mut guard = state.lock().map_err(|e| format!("state lock poisoned: {e}"))?;
            guard.volume_cache.get(&cache_key).is_some()
        };
        if already_cached {
            continue;
        }

        let app_cb = app.clone();
        let progress: Box<dyn Fn(usize, usize) + Send + Sync> = Box::new(move |done, total_slices| {
            let _ = app_cb.emit(
                "dicom_load_progress",
                ProgressEvent { phase: "decoding", done, total: total_slices, detail: None },
            );
        });

        let vol = match dicom_load::load_series(&series_dirs[idx], &descriptors[idx].uid, Some(progress)).await {
            Ok(v) => v,
            Err(e) => {
                failures.push(format!("{}: decode failed: {e}", descriptors[idx].name));
                continue;
            }
        };
        if let Err(e) = bridge_into_state(&vol, &state) {
            failures.push(format!("{}: bridge failed: {e}", descriptors[idx].name));
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
    }
    log_line!(
        "[load-timing] load_patient_all: decoded {} needed series (active + dual-energy pair) in {:.2?}",
        need.len(),
        decode_started.elapsed()
    );

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

        // build_dual_energy_volume checks the in-plane grid and aligns the two
        // series by slice z-position, so a num_slices difference (e.g. 429 vs
        // 428 from different VMI recons) is handled rather than rejected.
        match build_dual_energy_volume(&low_vol, &high_vol, low_kev, high_kev) {
            Ok(de) => {
                let n = de.low.shape()[0];
                let mut guard = state.lock().map_err(|e| format!("state lock poisoned: {e}"))?;
                guard.dual_energy = Some(de);
                if low_vol.metadata.num_slices != high_vol.metadata.num_slices {
                    log_line!(
                        "[dual-energy] paired {} keV ({} slices) + {} keV ({} slices) on {} common z-positions",
                        low_kev, low_vol.metadata.num_slices, high_kev, high_vol.metadata.num_slices, n
                    );
                }
            }
            Err(e) => {
                failures.push(format!("dual-energy pairing skipped: {e}"));
            }
        }
    }

    let _ = app.emit(
        "dicom_load_progress",
        ProgressEvent { phase: "done", done: total, total, detail: None },
    );

    Ok(PatientLoadResult { series: descriptors, active_index, failures })
}

/// Make `(dir, uid)` the active volume in Rust state (so Water/Lipid, MMD and
/// the FAI pipeline operate on it), decoding it on first access. Returns the
/// metadata plus a shared handle to the decoded voxels — callers decide whether
/// to ship the voxels to the frontend.
///
/// This is the single source of the activation logic: both `set_active_volume`
/// (returns voxels) and `set_active_volume_meta` (metadata only) go through it.
async fn activate_volume(
    dir: &str,
    uid: &str,
    state: &State<'_, Mutex<AppState>>,
) -> Result<(PipelineVolumeMetadata, Arc<Vec<i16>>), String> {
    let key = (dir.to_string(), uid.to_string());

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

    if let Some(hit) = cached_hit {
        return Ok(hit);
    }

    // Lazy decode-on-demand: `load_patient_all` scans every series but only
    // decodes the ones needed up front (active + dual-energy pair). Switching
    // to any other series decodes it here on first access, then it stays
    // cached. This is what keeps the initial load fast.
    let vol = dicom_load::load_series(Path::new(dir), uid, None)
        .await
        .map_err(|e| format!("decode failed: {e}"))?;
    bridge_into_state(&vol, state)?;
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
        key,
        CachedVolume {
            metadata: meta.clone(),
            voxels_i16: Arc::clone(&voxels_arc),
            volume: cached_volume,
        },
    );
    Ok((meta, voxels_arc))
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
    let (metadata, voxels) = activate_volume(&dir, &uid, &state).await?;
    let voxel_bytes: Vec<u8> = bytemuck::cast_slice(&voxels[..]).to_vec();
    let framed = encode_frame(&metadata, &voxel_bytes)?;
    Ok(Response::new(framed))
}

/// Activate a volume in Rust state and return ONLY its metadata (JSON), never
/// the voxels. The frontend calls this when cornerstone3D already holds the
/// volume in its GPU/scalar cache, so re-shipping 150–220 MB of voxels over IPC
/// — the dominant cost of switching between already-viewed series — is skipped
/// entirely. On a cornerstone cache miss the frontend falls back to
/// `set_active_volume` to obtain the voxels.
#[tauri::command]
pub async fn set_active_volume_meta(
    dir: String,
    uid: String,
    state: State<'_, Mutex<AppState>>,
) -> Result<PipelineVolumeMetadata, String> {
    let (metadata, _voxels) = activate_volume(&dir, &uid, &state).await?;
    Ok(metadata)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Regression: every patient folder contains a `MonoPlus_70keV` series, so a
    /// last-component-only key collided and loaded one patient's saved session
    /// onto another. The key must be patient-unique (full path embedded).
    #[test]
    fn patient_file_key_unique_across_patients_with_same_series_name() {
        let a = patient_file_key("/data/512143294/MonoPlus_70keV");
        let b = patient_file_key("/data/510829769/MonoPlus_70keV");
        assert_ne!(a, b, "different patients sharing a series name must not collide");
        assert!(a.ends_with(".json") && b.ends_with(".json"));
    }

    #[test]
    fn session_flags_complete_bundle() {
        // Mirrors src/lib/session.ts saveSession bundle shape.
        let bundle: serde_json::Value = serde_json::from_str(
            r#"{ "version":1, "savedAt":"t",
                 "seeds":{"activeVessel":"RCA","vessels":{}},
                 "fai":{"RCA":{"fai_mean_hu":-75.0}},
                 "mmd":{"summary":{"n":3},"surfaces":[]},
                 "wl":null }"#,
        )
        .unwrap();
        assert_eq!(session_flags(&bundle), (true, true));
    }

    #[test]
    fn session_flags_fai_only_is_not_complete() {
        let bundle: serde_json::Value =
            serde_json::from_str(r#"{"fai":{"RCA":{}},"mmd":{"summary":null}}"#).unwrap();
        assert_eq!(session_flags(&bundle), (true, false));
    }

    #[test]
    fn session_flags_empty_bundle() {
        let bundle: serde_json::Value =
            serde_json::from_str(r#"{"fai":null,"mmd":{"summary":null}}"#).unwrap();
        assert_eq!(session_flags(&bundle), (false, false));
    }

    #[test]
    fn status_for_rules() {
        // (complete, has_seeds) -> label. This is the bug regression: a bundle
        // with FAI+MMD must yield "complete".
        assert_eq!(status_for(true, true), "complete");
        assert_eq!(status_for(true, false), "complete");
        assert_eq!(status_for(false, true), "in_progress");
        assert_eq!(status_for(false, false), "not_started");
    }

    /// Minimal single-column volume: `z_positions.len()` slices, in-plane
    /// `rows×cols`, every voxel = `fill`.
    fn make_vol(z_positions: &[f64], rows: u32, cols: u32, fill: i16) -> PipelineLoadedVolume {
        let n = z_positions.len();
        let slice_len = (rows as usize) * (cols as usize);
        PipelineLoadedVolume {
            metadata: PipelineVolumeMetadata {
                series_uid: "uid".into(),
                series_description: "desc".into(),
                image_comments: None,
                rows,
                cols,
                num_slices: n,
                pixel_spacing: [1.0, 1.0],
                slice_spacing: 1.0,
                orientation: [1.0, 0.0, 0.0, 0.0, 1.0, 0.0],
                window_center: 40.0,
                window_width: 400.0,
                patient_name: "p".into(),
                study_description: "s".into(),
                slice_positions_z: z_positions.to_vec(),
                image_position_patient: [0.0, 0.0, z_positions.first().copied().unwrap_or(0.0)],
            },
            voxels_i16: vec![fill; n * slice_len],
        }
    }

    #[test]
    fn dual_energy_aligns_mismatched_slice_counts() {
        // 70 keV has 4 slices (z=0,1,2,3); 150 keV has 3 (z=0,1,2) — like the
        // real 429-vs-428 case. Pairing on z keeps the 3 common slices.
        let low = make_vol(&[0.0, 1.0, 2.0, 3.0], 1, 1, 70);
        let high = make_vol(&[0.0, 1.0, 2.0], 1, 1, 150);
        let de = build_dual_energy_volume(&low, &high, 70.0, 150.0)
            .expect("should pair on common z-positions");
        assert_eq!(de.low.shape(), &[3, 1, 1], "keeps the 3 common slices");
        assert_eq!(de.high.shape(), &[3, 1, 1]);
        assert!(de.low.iter().all(|&v| v == 70.0), "low co-registered");
        assert!(de.high.iter().all(|&v| v == 150.0), "high co-registered");
        assert!((de.origin[0] - 0.0).abs() < 1e-9, "origin z = first matched slice");
    }

    #[test]
    fn dual_energy_rejects_different_inplane_grid() {
        let low = make_vol(&[0.0, 1.0], 1, 2, 70); // 1×2 in-plane
        let high = make_vol(&[0.0, 1.0], 1, 1, 150); // 1×1 in-plane
        assert!(
            build_dual_energy_volume(&low, &high, 70.0, 150.0).is_err(),
            "different in-plane grids must be rejected"
        );
    }

    #[test]
    fn dual_energy_errors_when_no_z_overlap() {
        let low = make_vol(&[0.0, 1.0], 1, 1, 70);
        let high = make_vol(&[100.0, 101.0], 1, 1, 150); // no shared z
        assert!(build_dual_energy_volume(&low, &high, 70.0, 150.0).is_err());
    }
}
