//! End-to-end test: exercise `load_dicom_directory` against the real local
//! patient `57955439` for every series subfolder. Verifies the loader works
//! before plumbing it through Tauri commands.
//!
//! Run with: `cargo test --test load_patient_57955439 -- --nocapture`

use std::path::Path;

use pcat_pipeline::dicom_loader::{load_dicom_directory, scan_dicom_series};

const PATIENT_ROOT: &str = "/Users/shunie/Developer/PCAT/UCI NAEOTOM CCTA Data/57955439";

#[test]
fn each_series_loads_cleanly() {
    let series_dirs = [
        "MonoPlus_70keV",
        "MonoPlus_150keV",
        "CCTA_Bv44__CORONARY  Bv44 Q2 75%",
        "CA_SCORING",
    ];

    for name in &series_dirs {
        let path = format!("{}/{}", PATIENT_ROOT, name);
        let p = Path::new(&path);
        if !p.exists() {
            println!("SKIP  {} — path missing", name);
            continue;
        }

        match load_dicom_directory(p) {
            Ok(vol) => {
                let shape = vol.data.shape();
                println!(
                    "OK    {:30} shape=[{},{},{}]  spacing=[{:.3},{:.3},{:.3}]  origin=[{:.2},{:.2},{:.2}]",
                    name,
                    shape[0],
                    shape[1],
                    shape[2],
                    vol.spacing[0],
                    vol.spacing[1],
                    vol.spacing[2],
                    vol.origin[0],
                    vol.origin[1],
                    vol.origin[2],
                );
                let mid = vol.data[[shape[0] / 2, shape[1] / 2, shape[2] / 2]];
                println!("      center voxel HU = {:.1}", mid);
                assert!(shape[0] > 1 && shape[1] > 1 && shape[2] > 1, "degenerate volume");
            }
            Err(e) => {
                println!("FAIL  {} — {}", name, e);
                panic!("load_dicom_directory failed for {}: {}", name, e);
            }
        }
    }
}

#[test]
fn scan_patient_root_lists_all_series() {
    let p = Path::new(PATIENT_ROOT);
    if !p.exists() {
        println!("SKIP — patient root missing");
        return;
    }

    let series = scan_dicom_series(p).expect("scan_dicom_series must succeed");
    println!("Found {} series under patient root:", series.len());
    for s in &series {
        let kev = s.kev_label.map(|k| format!("{:.0} keV", k)).unwrap_or_else(|| "—".into());
        println!(
            "  uid={}  desc={:?}  slices={}  kev={}",
            &s.series_uid[..s.series_uid.len().min(40)],
            s.description,
            s.num_slices,
            kev,
        );
    }
    assert!(series.len() >= 2, "expected at least 2 series under patient root");
}
