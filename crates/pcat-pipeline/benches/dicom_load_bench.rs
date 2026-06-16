//! Wall-clock benchmark for scan_series + load_series.
//!
//! By default, runs against the synthetic fixture (30 slices × 64²).
//! To measure real SMB performance, set BENCH_DICOM_DIR:
//!
//!     BENCH_DICOM_DIR=/Volumes/labshare/patient_xyz cargo bench -p pcat-pipeline

use std::path::PathBuf;
use std::time::Instant;

use pcat_pipeline::dicom_load::load_series;
use pcat_pipeline::dicom_scan::scan_series;

#[path = "../tests/common/mod.rs"]
mod fixture;

fn main() {
    let runtime = tokio::runtime::Runtime::new().unwrap();
    runtime.block_on(async {
        let (dir_owner, dir) = resolve_dir();
        eprintln!("Benchmarking dir: {}", dir.display());

        let t0 = Instant::now();
        let series = scan_series(&dir).await.expect("scan failed");
        let t_scan = t0.elapsed();

        println!("scan_series: {:.3} s ({} series)", t_scan.as_secs_f64(), series.len());
        for (i, s) in series.iter().enumerate() {
            println!("  [{i}] {} ({} slices, {}×{})", s.description, s.num_slices, s.rows, s.cols);
        }

        let first = series.first().expect("at least one series");
        let t1 = Instant::now();
        let vol = load_series(&dir, &first.uid, None).await.expect("load failed");
        let t_load = t1.elapsed();
        let total = t_scan + t_load;
        let mb = (vol.voxels_i16.len() * 2) / (1024 * 1024);
        println!(
            "load_series: {:.3} s ({} slices, {} MB)",
            t_load.as_secs_f64(), vol.metadata.num_slices, mb,
        );
        println!("TOTAL cold load: {:.3} s", total.as_secs_f64());

        drop(dir_owner);
    });
}

fn resolve_dir() -> (Option<tempfile::TempDir>, PathBuf) {
    if let Ok(d) = std::env::var("BENCH_DICOM_DIR") {
        return (None, PathBuf::from(d));
    }
    let td = tempfile::tempdir().expect("tempdir");
    fixture::write_mini_ct(td.path());
    let p = td.path().to_path_buf();
    (Some(td), p)
}
