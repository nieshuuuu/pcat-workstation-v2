//! Smoke test for the bulk load data shape. Runs the pipeline directly
//! (not via Tauri) and asserts the frame can be parsed back.

use pcat_pipeline::dicom_load::load_series;
use pcat_pipeline::dicom_scan::scan_series;

#[path = "../../crates/pcat-pipeline/tests/common/mod.rs"]
mod fixture;

#[tokio::test]
async fn bulk_load_frame_round_trip() {
    let dir = tempfile::tempdir().unwrap();
    fixture::write_mini_ct(dir.path());

    let series = scan_series(dir.path()).await.unwrap();
    let vol = load_series(dir.path(), &series[0].uid, None).await.unwrap();

    let payload: Vec<u8> = bytemuck::cast_slice(&vol.voxels_i16).to_vec();
    let frame = pcat_workstation_v2_lib::commands::framed::encode_frame(&vol.metadata, &payload)
        .expect("frame");

    let meta_len = u32::from_le_bytes(frame[..4].try_into().unwrap()) as usize;
    let meta_json = &frame[4..4 + meta_len];
    let body = &frame[4 + meta_len..];

    let meta: pcat_pipeline::dicom_load::VolumeMetadata =
        serde_json::from_slice(meta_json).unwrap();
    assert_eq!(meta.num_slices, fixture::FIXTURE_SLICES);
    assert_eq!(body.len(), vol.voxels_i16.len() * 2);
}
