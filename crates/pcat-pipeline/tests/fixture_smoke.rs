//! Smoke test: fixture generator produces files that read_header can parse.

mod common;

use pcat_pipeline::dicom_scan::read_header;

#[test]
fn fixture_headers_parse() {
    let dir = tempfile::tempdir().unwrap();
    let paths = common::write_mini_ct(dir.path());
    assert_eq!(paths.len(), common::FIXTURE_SLICES);

    for (i, p) in paths.iter().enumerate() {
        let h = read_header(p).unwrap().expect("fixture should parse");
        assert_eq!(h.rows, common::FIXTURE_ROWS);
        assert_eq!(h.cols, common::FIXTURE_COLS);
        assert_eq!(h.series_uid, "1.2.826.0.1.3680043.99.1");
        assert_eq!(h.rescale_slope, 1.0);
        assert_eq!(h.rescale_intercept, -1024.0);
        let expected_z = i as f64 * common::FIXTURE_SPACING_Z;
        assert!((h.image_position_z.unwrap() - expected_z).abs() < 1e-9);
    }
}

#[test]
fn fixture_mixed_series() {
    let dir = tempfile::tempdir().unwrap();
    let a = common::write_mini_ct(dir.path());
    let b = common::write_second_series(dir.path());
    assert_eq!(a.len() + b.len(), 60);

    let h0 = read_header(&a[0]).unwrap().unwrap();
    let h1 = read_header(&b[0]).unwrap().unwrap();
    assert_ne!(h0.series_uid, h1.series_uid);
    assert_eq!(h1.image_comments.as_deref(), Some("E = 70 keV"));
}
