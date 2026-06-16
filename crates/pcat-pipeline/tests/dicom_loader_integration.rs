mod common;

use pcat_pipeline::dicom_scan::scan_series;

#[tokio::test]
async fn scan_returns_one_series_for_mini_ct() {
    let dir = tempfile::tempdir().unwrap();
    common::write_mini_ct(dir.path());

    let series = scan_series(dir.path()).await.unwrap();
    assert_eq!(series.len(), 1);

    let s = &series[0];
    assert_eq!(s.uid, "1.2.826.0.1.3680043.99.1");
    assert_eq!(s.description, "MiniCT");
    assert_eq!(s.num_slices, common::FIXTURE_SLICES);
    assert_eq!(s.rows, common::FIXTURE_ROWS);
    assert_eq!(s.cols, common::FIXTURE_COLS);
    assert_eq!(s.file_paths.len(), common::FIXTURE_SLICES);

    // Files are z-sorted
    let zs: Vec<f64> = (0..common::FIXTURE_SLICES)
        .map(|i| i as f64 * common::FIXTURE_SPACING_Z)
        .collect();
    for (i, p) in s.file_paths.iter().enumerate() {
        assert!(
            p.to_string_lossy().contains(&format!("slice_{:03}", i)),
            "slice at index {i} is {p:?}, expected contains slice_{:03}",
            i
        );
        assert!((s.slice_positions_z[i] - zs[i]).abs() < 1e-9);
    }
}

#[tokio::test]
async fn scan_returns_two_series_for_dual_folder() {
    let dir = tempfile::tempdir().unwrap();
    common::write_mini_ct(dir.path());
    common::write_second_series(dir.path());

    let mut series = scan_series(dir.path()).await.unwrap();
    series.sort_by(|a, b| a.uid.cmp(&b.uid));
    assert_eq!(series.len(), 2);
    assert_eq!(series[1].description, "MonoPlus 70 keV");
    assert_eq!(series[1].image_comments.as_deref(), Some("E = 70 keV"));
}

#[tokio::test]
async fn scan_silently_skips_non_dicom_files() {
    let dir = tempfile::tempdir().unwrap();
    common::write_mini_ct(dir.path());
    std::fs::write(dir.path().join(".DS_Store"), b"junk").unwrap();
    std::fs::write(dir.path().join("README.txt"), b"hello").unwrap();

    let series = scan_series(dir.path()).await.unwrap();
    assert_eq!(series.len(), 1);
    assert_eq!(series[0].num_slices, common::FIXTURE_SLICES);
}

#[tokio::test]
async fn scan_errors_when_folder_missing() {
    let err = scan_series(std::path::Path::new("/nonexistent/path/xyz"))
        .await
        .unwrap_err();
    let s = format!("{err}");
    assert!(s.contains("folder not readable") || s.contains("No such file"));
}

#[tokio::test]
async fn decode_slice_returns_expected_hu() {
    use pcat_pipeline::dicom_decode::decode_slice_i16;
    use pcat_pipeline::dicom_scan::scan_series;

    let dir = tempfile::tempdir().unwrap();
    common::write_mini_ct(dir.path());

    let series = scan_series(dir.path()).await.unwrap();
    let desc = &series[0];
    let px = decode_slice_i16(
        &desc.file_paths[0],
        desc.rescale_slope,
        desc.rescale_intercept,
        desc.rows,
        desc.cols,
    )
    .unwrap();

    assert_eq!(px.len(), (desc.rows * desc.cols) as usize);
    // Fixture center: (1024 - 0*30 + 0*5) * 1 + (-1024) = 0 HU exactly
    let center_idx = ((desc.rows / 2) * desc.cols + desc.cols / 2) as usize;
    assert!(
        (px[center_idx] - 0).abs() < 2,
        "center HU should be ~0, got {}",
        px[center_idx]
    );
    // Corner: (1024 - sqrt(2*32^2)*30 + 0) * 1 + (-1024) = very negative
    assert!(px[0] < -200, "corner should be negative HU, got {}", px[0]);
}

#[tokio::test]
async fn load_series_returns_correct_volume() {
    use pcat_pipeline::dicom_load::{load_series, LoadedVolume};
    use pcat_pipeline::dicom_scan::scan_series;

    let dir = tempfile::tempdir().unwrap();
    common::write_mini_ct(dir.path());

    let series = scan_series(dir.path()).await.unwrap();
    let vol: LoadedVolume = load_series(dir.path(), &series[0].uid, None).await.unwrap();

    assert_eq!(vol.metadata.num_slices, common::FIXTURE_SLICES);
    assert_eq!(vol.metadata.rows, common::FIXTURE_ROWS);
    assert_eq!(vol.metadata.cols, common::FIXTURE_COLS);
    let slice_len = (common::FIXTURE_ROWS * common::FIXTURE_COLS) as usize;
    assert_eq!(vol.voxels_i16.len(), slice_len * common::FIXTURE_SLICES);

    // Center of first slice = HU ~0 (see Task 7 test)
    let center_in_first = (common::FIXTURE_ROWS / 2) as usize
        * common::FIXTURE_COLS as usize
        + (common::FIXTURE_COLS / 2) as usize;
    assert!((vol.voxels_i16[center_in_first] - 0).abs() < 2);

    // Center of last slice gets z-offset: (1024 + 29*5) - 1024 = 145 HU
    let last_slice_offset = slice_len * (common::FIXTURE_SLICES - 1);
    assert!(
        (vol.voxels_i16[last_slice_offset + center_in_first] - 145).abs() < 2,
        "last slice center HU ~145, got {}",
        vol.voxels_i16[last_slice_offset + center_in_first]
    );
}

#[tokio::test]
async fn load_series_errors_on_unknown_uid() {
    let dir = tempfile::tempdir().unwrap();
    common::write_mini_ct(dir.path());

    let err = pcat_pipeline::dicom_load::load_series(dir.path(), "not-a-real-uid", None)
        .await
        .unwrap_err();
    assert!(format!("{err}").contains("not-a-real-uid"));
}

#[tokio::test]
async fn load_series_rejects_oversized_volume() {
    use pcat_pipeline::dicom_load::check_volume_size_mb;
    let mb = 1024 * 1024 * 1024usize; // 1 GB hypothetical
    assert!(check_volume_size_mb(mb).is_err());
}
