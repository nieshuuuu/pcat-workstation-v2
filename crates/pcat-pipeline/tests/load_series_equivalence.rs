//! Byte-equivalence guard for `load_series`. Captures the decoded volume's
//! shape, z-ordering, and a voxel checksum on a real local series so a change
//! to the decode path (e.g. single-pass read) can be proven not to alter the
//! reconstructed volume.

use std::path::Path;

use pcat_pipeline::dicom_load::load_series;
use pcat_pipeline::dicom_scan::scan_series;

const SERIES: &str = "/Users/shunie/Developer/PCAT/UCI NAEOTOM CCTA Data/57955439/MonoPlus_70keV";

#[tokio::test]
async fn load_series_reference_checksum() {
    let dir = Path::new(SERIES);
    if !dir.exists() {
        eprintln!("SKIP — series missing");
        return;
    }
    let descs = scan_series(dir).await.expect("scan");
    let uid = descs[0].uid.clone();
    let vol = load_series(dir, &uid, None).await.expect("load");
    let m = &vol.metadata;
    let len = vol.voxels_i16.len();
    let sum: i64 = vol.voxels_i16.iter().map(|&v| v as i64).sum();
    let first = vol.voxels_i16[0];
    let mid = vol.voxels_i16[len / 2];
    let last = *vol.voxels_i16.last().unwrap();
    eprintln!(
        "VALUES num_slices={} rows={} cols={} len={} sum={} first={} mid={} last={} z0={} zlast={} spacing={}",
        m.num_slices, m.rows, m.cols, len, sum, first, mid, last,
        m.slice_positions_z[0], m.slice_positions_z.last().unwrap(), m.slice_spacing
    );

    // Reference captured from the scan-then-decode implementation on this series.
    // The single-pass loader MUST reproduce these exactly (byte-identical volume,
    // identical z-ordering) — the integer checksums are the strong guarantee.
    assert_eq!(m.num_slices, 339, "num_slices");
    assert_eq!(m.rows, 512, "rows");
    assert_eq!(m.cols, 513, "cols");
    assert_eq!(len, 89_040_384, "voxel count");
    assert_eq!(sum, -22_095_243_454, "voxel checksum (sum)");
    assert_eq!(first, -913, "first voxel");
    assert_eq!(mid, -488, "mid voxel");
    assert_eq!(last, -92, "last voxel");
    assert!((m.slice_positions_z[0] - 1487.9196).abs() < 1e-3, "z0");
    assert!((m.slice_positions_z.last().unwrap() - 1607.0804).abs() < 1e-3, "zlast");
    assert!((m.slice_spacing - 0.3525).abs() < 1e-3, "spacing");
}
