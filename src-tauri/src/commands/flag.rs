//! Patient-level data-quality flag ("Flagged" — data has a problem / unusable).
//! The only human-authored fact in the patient-status feature, so it is the only
//! thing persisted: one tiny sidecar `flags/{patient_file_key}.json` per patient.
//!
//! `dicom_path` everywhere in the app is a *series* subfolder; flags are
//! patient-level, so we key off the parent (patient) folder — every series of a
//! patient shares one flag file. `list_patients` reads the same key.

use std::path::{Path, PathBuf};
use tauri::Manager;

use crate::commands::dicom::patient_file_key;

#[derive(serde::Serialize, serde::Deserialize, Clone)]
pub struct PatientFlag {
    pub flagged: bool,
    pub note: String,
    pub flagged_at: String,
}

/// The patient folder that owns a loaded *series* path (its parent dir).
pub(crate) fn patient_dir_of(series_path: &str) -> String {
    Path::new(series_path)
        .parent()
        .map(|p| p.to_string_lossy().to_string())
        .filter(|s| !s.is_empty())
        .unwrap_or_else(|| series_path.to_string())
}

/// `<base>/<patient_file_key(patient_dir)>` — the same key `list_patients` uses,
/// so a flag set while viewing any series is found when listing the patient.
pub(crate) fn flag_file_path(base_dir: &Path, patient_dir: &str) -> PathBuf {
    base_dir.join(patient_file_key(patient_dir))
}

/// Read a flag file, returning None if absent. Errors only on real IO/parse
/// failures (never swallows them into None).
pub(crate) fn read_flag_at(path: &Path) -> Result<Option<PatientFlag>, String> {
    match std::fs::read_to_string(path) {
        Ok(s) => {
            let f: PatientFlag =
                serde_json::from_str(&s).map_err(|e| format!("parse failed: {e}"))?;
            Ok(Some(f))
        }
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => Ok(None),
        Err(e) => Err(format!("read failed: {e}")),
    }
}

fn flags_dir(app: &tauri::AppHandle) -> PathBuf {
    let dir = app
        .path()
        .app_data_dir()
        .expect("app data dir")
        .join("flags");
    let _ = std::fs::create_dir_all(&dir);
    dir
}

#[tauri::command]
pub async fn set_patient_flag(
    app: tauri::AppHandle,
    dicom_path: String,
    flagged: bool,
    note: String,
    flagged_at: String,
) -> Result<(), String> {
    let patient_dir = patient_dir_of(&dicom_path);
    let path = flag_file_path(&flags_dir(&app), &patient_dir);
    if flagged {
        let flag = PatientFlag {
            flagged: true,
            note,
            flagged_at,
        };
        let json =
            serde_json::to_string_pretty(&flag).map_err(|e| format!("serialize failed: {e}"))?;
        std::fs::write(&path, json).map_err(|e| format!("write failed: {e}"))?;
    } else {
        match std::fs::remove_file(&path) {
            Ok(()) => {}
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => {}
            Err(e) => return Err(format!("remove failed: {e}")),
        }
    }
    Ok(())
}

#[tauri::command]
pub async fn get_patient_flag(
    app: tauri::AppHandle,
    dicom_path: String,
) -> Result<Option<PatientFlag>, String> {
    let patient_dir = patient_dir_of(&dicom_path);
    let path = flag_file_path(&flags_dir(&app), &patient_dir);
    read_flag_at(&path)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn flag_is_patient_level_two_series_one_file() {
        let base = Path::new("/base");
        let a = flag_file_path(base, &patient_dir_of("/data/P001/MonoPlus_70keV"));
        let b = flag_file_path(base, &patient_dir_of("/data/P001/CCTA"));
        assert_eq!(a, b, "two series of one patient must map to one flag file");
        assert!(a.starts_with("/base"));
    }

    #[test]
    fn flag_distinct_patients_distinct_files() {
        let base = Path::new("/base");
        let a = flag_file_path(base, &patient_dir_of("/data/P001/MonoPlus_70keV"));
        let b = flag_file_path(base, &patient_dir_of("/data/P002/MonoPlus_70keV"));
        assert_ne!(a, b);
    }

    #[test]
    fn flag_round_trip_write_read_delete() {
        let dir = std::env::temp_dir().join("pcat_flag_test_round_trip");
        let _ = std::fs::create_dir_all(&dir);
        let path = dir.join("p.json");
        let _ = std::fs::remove_file(&path);
        assert!(read_flag_at(&path).unwrap().is_none());

        let flag = PatientFlag {
            flagged: true,
            note: "vessel jump ~slice 40".into(),
            flagged_at: "2026-06-17T00:00:00Z".into(),
        };
        std::fs::write(&path, serde_json::to_string(&flag).unwrap()).unwrap();
        let got = read_flag_at(&path).unwrap().unwrap();
        assert!(got.flagged);
        assert_eq!(got.note, "vessel jump ~slice 40");

        std::fs::remove_file(&path).unwrap();
        assert!(read_flag_at(&path).unwrap().is_none());
    }
}
