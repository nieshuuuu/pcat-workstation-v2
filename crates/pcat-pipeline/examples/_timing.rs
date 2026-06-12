//! THROWAWAY timing harness — measures scan + decode time on a patient folder
//! to locate the volume-loading bottleneck. Deleted after use.

use std::path::Path;
use std::time::Instant;

use pcat_pipeline::{dicom_load, dicom_scan};

fn main() {
    let rt = tokio::runtime::Runtime::new().unwrap();
    rt.block_on(async {
        let pat = "/Users/shunie/Developer/PCAT/UCI NAEOTOM CCTA Data/57955439";
        let subdirs = [
            "CA_SCORING",
            "CCTA_Bv44__CORONARY  Bv44 Q2 75%",
            "MonoPlus_150keV",
            "MonoPlus_70keV",
        ];

        println!("\n=== SCAN (header-only) per series ===");
        let t_all = Instant::now();
        for sd in subdirs {
            let dir = format!("{pat}/{sd}");
            let t = Instant::now();
            match dicom_scan::scan_series(Path::new(&dir)).await {
                Ok(s) => println!(
                    "  scan {:<32} {:>8.2?}  ({} series, {} slices)",
                    sd,
                    t.elapsed(),
                    s.len(),
                    s.first().map(|d| d.num_slices).unwrap_or(0)
                ),
                Err(e) => println!("  scan {sd:<32} ERR {e}"),
            }
        }
        println!("  scan ALL 4 series: {:.2?}", t_all.elapsed());

        println!("\n=== DECODE (load_series) ===");
        for sd in ["MonoPlus_70keV", "CCTA_Bv44__CORONARY  Bv44 Q2 75%"] {
            let dir = format!("{pat}/{sd}");
            let scan = match dicom_scan::scan_series(Path::new(&dir)).await {
                Ok(s) => s,
                Err(e) => {
                    println!("  decode {sd}: scan ERR {e}");
                    continue;
                }
            };
            let uid = scan[0].uid.clone();
            let n = scan[0].num_slices.max(1);
            let t = Instant::now();
            match dicom_load::load_series(Path::new(&dir), &uid, None).await {
                Ok(vol) => {
                    let dt = t.elapsed();
                    println!(
                        "  decode {:<32} {:>8.2?}  ({} slices, {:.1} ms/slice, {} MB)",
                        sd,
                        dt,
                        vol.metadata.num_slices,
                        dt.as_secs_f64() * 1000.0 / n as f64,
                        vol.voxels_i16.len() * 2 / (1024 * 1024)
                    );
                }
                Err(e) => println!("  decode {sd}: ERR {e}"),
            }
        }

        println!("\n=== what lazy load_patient_all decodes (active CCTA + 70/150 pair) ===");
        let t = Instant::now();
        for sd in ["CCTA_Bv44__CORONARY  Bv44 Q2 75%", "MonoPlus_70keV", "MonoPlus_150keV"] {
            let dir = format!("{pat}/{sd}");
            let scan = dicom_scan::scan_series(Path::new(&dir)).await.unwrap();
            let _ = dicom_load::load_series(Path::new(&dir), &scan[0].uid, None).await.unwrap();
        }
        println!("  active + DE pair (3 series) total: {:.2?}", t.elapsed());
    });
}
