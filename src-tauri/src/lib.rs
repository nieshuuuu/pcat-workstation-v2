pub mod commands;
mod state;
pub mod volume_cache;

use state::AppState;
use std::sync::Mutex;
use tauri::Manager;

#[cfg_attr(mobile, tauri::mobile_entry_point)]
pub fn run() {
    tauri::Builder::default()
        .plugin(tauri_plugin_dialog::init())
        .setup(|app| {
            app.manage(Mutex::new(AppState::new()));
            {
                use tauri::Manager;
                let home = std::env::var("HOME").unwrap_or_default();
                let add = app.path().app_data_dir();
                let mut log = format!("APP_DATA_DIR={add:?}\n");
                let root = "/Volumes/Molloilab/Shu Nie/UCI NAEOTOM CCTA Data";
                if let Ok(rd) = std::fs::read_dir(root) {
                    for e in rd.flatten() {
                        if e.file_name().to_string_lossy() != "512339528" { continue; }
                        let pp = e.path().to_string_lossy().to_string();
                        let needle: String = pp.chars()
                            .map(|c| if c.is_ascii_alphanumeric() || c=='-' || c=='_' || c=='.' { c } else { '_' })
                            .collect::<String>().replace("..", "_");
                        log.push_str(&format!("patient_path={pp}\nneedle={needle}\n"));
                        if let Ok(ad) = &add {
                            let sess = ad.join("sessions");
                            log.push_str(&format!("sessions_dir={sess:?} exists={}\n", sess.is_dir()));
                            let (mut nm, mut complete) = (0, false);
                            if let Ok(srd) = std::fs::read_dir(&sess) {
                                for se in srd.flatten() {
                                    let n = se.file_name().to_string_lossy().to_string();
                                    if !n.contains(&needle) { continue; }
                                    nm += 1;
                                    match std::fs::read_to_string(se.path()).ok()
                                        .and_then(|d| serde_json::from_str::<serde_json::Value>(&d).ok()) {
                                        Some(j) => {
                                            let f = j.get("fai").and_then(|v| v.as_object()).map(|o| !o.is_empty()).unwrap_or(false);
                                            let m = j.get("mmd").and_then(|x| x.get("summary")).map(|s| !s.is_null()).unwrap_or(false);
                                            log.push_str(&format!("MATCH {n} fai={f} mmd={m}\n"));
                                            complete |= f && m;
                                        }
                                        None => log.push_str(&format!("READ/PARSE FAIL {n}\n")),
                                    }
                                }
                            } else { log.push_str("sessions read_dir FAILED\n"); }
                            log.push_str(&format!("nm={nm} complete={complete}\n"));
                        }
                    }
                } else { log.push_str("root read_dir FAILED\n"); }
                let _ = std::fs::write(format!("{home}/pcat_appctx.txt"), log);
            }
            Ok(())
        })
        .invoke_handler(tauri::generate_handler![
            commands::dicom::open_dicom_dialog,
            commands::dicom::get_recent_dicoms,
            commands::dicom::save_seeds,
            commands::dicom::load_seeds,
            commands::dicom::scan_series,
            commands::dicom::load_series,
            commands::dicom::reuse_loaded_volume,
            commands::dicom::load_dual_energy,
            commands::dicom::load_patient_all,
            commands::dicom::set_active_volume,
            commands::dicom::set_active_volume_meta,
            commands::dicom::list_patients,
            commands::dicom::list_series_dirs,
            commands::cpr::build_cpr_frame,
            commands::cpr::render_cpr_image,
            commands::cpr::render_stretched_cpr_image,
            commands::cpr::render_cross_sections,
            commands::cpr::compute_cpr_image,
            commands::cpr::compute_cross_section_image,
            commands::cpr::compute_cross_sections_batch,
            commands::cpr::get_cpr_projection_info,
            commands::pipeline::run_pipeline,
            commands::annotation::generate_annotation_targets,
            commands::annotation::init_snake,
            commands::annotation::evolve_snake,
            commands::annotation::finalize_contour,
            commands::annotation::use_vessel_wall_as_contour,
            commands::annotation::run_mmd_on_roi,
            commands::annotation::sample_surfaces,
            commands::annotation::get_mmd_overlay,
            commands::annotation::save_annotations,
            commands::annotation::load_annotations,
            commands::annotation::save_session,
            commands::annotation::load_session,
            commands::annotation::export_mmd_csv,
            commands::water_lipid::run_water_lipid,
            commands::water_lipid::get_wl_slice,
            commands::water_lipid::restore_wl_calibration,
            commands::flag::set_patient_flag,
            commands::flag::get_patient_flag,
        ])
        .run(tauri::generate_context!())
        .expect("error while running tauri application");
}
