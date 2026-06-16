//! Test fixtures — generates valid synthetic DICOM files for integration tests.
//!
//! The fixture is a 30-slice mini-CT at 64x64 resolution with a simple radial
//! gradient pattern. Each slice is a valid Explicit-VR Little Endian DICOM file
//! with uncompressed 16-bit pixel data, rescale slope/intercept set, and
//! ImagePositionPatient encoding the Z axis.

use std::path::{Path, PathBuf};

use dicom::core::{DataElement, PrimitiveValue, VR};
use dicom::dictionary_std::tags;
use dicom::object::{FileMetaTableBuilder, InMemDicomObject};

/// Dimensions of the synthetic fixture.
pub const FIXTURE_ROWS: u32 = 64;
pub const FIXTURE_COLS: u32 = 64;
pub const FIXTURE_SLICES: usize = 30;
pub const FIXTURE_SPACING_XY: f64 = 0.5;
pub const FIXTURE_SPACING_Z: f64 = 1.0;

/// Write 30 synthetic DICOM slices into `dir`. Returns the list of created
/// file paths in z-ascending order. Uses series UID `1.2.826.0.1.3680043.99.1`.
pub fn write_mini_ct(dir: &Path) -> Vec<PathBuf> {
    write_series(dir, "1.2.826.0.1.3680043.99.1", "MiniCT", None, 0)
}

/// Write a second series to the same folder (for multi-series tests).
pub fn write_second_series(dir: &Path) -> Vec<PathBuf> {
    write_series(
        dir,
        "1.2.826.0.1.3680043.99.2",
        "MonoPlus 70 keV",
        Some("E = 70 keV"),
        1000,
    )
}

fn write_series(
    dir: &Path,
    series_uid: &str,
    description: &str,
    image_comments: Option<&str>,
    filename_offset: usize,
) -> Vec<PathBuf> {
    let study_uid = "1.2.826.0.1.3680043.99.0";

    (0..FIXTURE_SLICES)
        .map(|z| {
            let sop_uid = format!("1.2.826.0.1.3680043.99.9.{}.{}", series_uid, z);
            let mut obj = InMemDicomObject::new_empty();

            // Patient/Study/Series
            obj.put(elem(tags::PATIENT_NAME, VR::PN, "SYNTH^PHANTOM"));
            obj.put(elem(tags::PATIENT_ID, VR::LO, "FIXTURE-1"));
            obj.put(elem(tags::STUDY_INSTANCE_UID, VR::UI, study_uid));
            obj.put(elem(tags::STUDY_DESCRIPTION, VR::LO, "Synthetic study"));
            obj.put(elem(tags::SERIES_INSTANCE_UID, VR::UI, series_uid));
            obj.put(elem(tags::SERIES_DESCRIPTION, VR::LO, description));
            obj.put(elem(tags::SOP_INSTANCE_UID, VR::UI, &sop_uid));
            obj.put(elem(tags::SOP_CLASS_UID, VR::UI, "1.2.840.10008.5.1.4.1.1.2"));
            obj.put(elem(tags::MODALITY, VR::CS, "CT"));

            // Image pixel module
            obj.put(elem_int(tags::ROWS, VR::US, FIXTURE_ROWS));
            obj.put(elem_int(tags::COLUMNS, VR::US, FIXTURE_COLS));
            obj.put(elem_int(tags::BITS_ALLOCATED, VR::US, 16u32));
            obj.put(elem_int(tags::BITS_STORED, VR::US, 16u32));
            obj.put(elem_int(tags::HIGH_BIT, VR::US, 15u32));
            obj.put(elem_int(tags::PIXEL_REPRESENTATION, VR::US, 1u32));
            obj.put(elem_int(tags::SAMPLES_PER_PIXEL, VR::US, 1u32));
            obj.put(elem(
                tags::PHOTOMETRIC_INTERPRETATION,
                VR::CS,
                "MONOCHROME2",
            ));

            // Geometry
            obj.put(elem_multi_f64(
                tags::PIXEL_SPACING,
                &[FIXTURE_SPACING_XY, FIXTURE_SPACING_XY],
            ));
            obj.put(elem_multi_f64(
                tags::IMAGE_ORIENTATION_PATIENT,
                &[1.0, 0.0, 0.0, 0.0, 1.0, 0.0],
            ));
            let z_mm = (z as f64) * FIXTURE_SPACING_Z;
            obj.put(elem_multi_f64(
                tags::IMAGE_POSITION_PATIENT,
                &[0.0, 0.0, z_mm],
            ));
            obj.put(elem(tags::INSTANCE_NUMBER, VR::IS, &format!("{}", z + 1)));
            obj.put(elem(
                tags::SLICE_THICKNESS,
                VR::DS,
                &format!("{}", FIXTURE_SPACING_Z),
            ));

            // HU rescale (CT: slope 1, intercept -1024)
            obj.put(elem(tags::RESCALE_SLOPE, VR::DS, "1"));
            obj.put(elem(tags::RESCALE_INTERCEPT, VR::DS, "-1024"));

            // Window defaults
            obj.put(elem(tags::WINDOW_CENTER, VR::DS, "40"));
            obj.put(elem(tags::WINDOW_WIDTH, VR::DS, "400"));

            // Optional ImageComments (for MonoPlus keV fixture)
            if let Some(c) = image_comments {
                obj.put(elem(IMAGE_COMMENTS_TAG, VR::LT, c));
            }

            // Pixel data: radial gradient centered on slice midpoint
            let mut pixels = vec![0i16; (FIXTURE_ROWS * FIXTURE_COLS) as usize];
            for r in 0..FIXTURE_ROWS {
                for c in 0..FIXTURE_COLS {
                    let dx = (c as f64) - (FIXTURE_COLS as f64) / 2.0;
                    let dy = (r as f64) - (FIXTURE_ROWS as f64) / 2.0;
                    let d = (dx * dx + dy * dy).sqrt();
                    // HU-like: +1024 at center fades to 0 at edge, plus z offset
                    let raw = (1024.0 - d * 30.0 + (z as f64) * 5.0)
                        .max(i16::MIN as f64)
                        .min(i16::MAX as f64) as i16;
                    pixels[(r * FIXTURE_COLS + c) as usize] = raw;
                }
            }
            let pixel_bytes: Vec<u8> = pixels
                .iter()
                .flat_map(|v| v.to_le_bytes())
                .collect();
            obj.put(DataElement::new(
                tags::PIXEL_DATA,
                VR::OW,
                PrimitiveValue::from(pixel_bytes),
            ));

            // Build file meta (Explicit VR Little Endian)
            let meta = FileMetaTableBuilder::new()
                .transfer_syntax("1.2.840.10008.1.2.1") // Explicit VR Little Endian
                .media_storage_sop_class_uid("1.2.840.10008.5.1.4.1.1.2")
                .media_storage_sop_instance_uid(&sop_uid);

            let file_obj = obj.with_meta(meta).unwrap();
            let path = dir.join(format!("slice_{:03}.dcm", z + filename_offset));
            file_obj.write_to_file(&path).unwrap();
            path
        })
        .collect()
}

/// IMAGE_COMMENTS tag re-export for the fixture helpers to avoid crate-private imports.
pub const IMAGE_COMMENTS_TAG: dicom::core::Tag = dicom::core::Tag(0x0020, 0x4000);

// Helpers to build DataElements succinctly.
fn elem(tag: dicom::core::Tag, vr: VR, s: &str) -> DataElement<InMemDicomObject> {
    DataElement::new(tag, vr, PrimitiveValue::from(s.to_string()))
}
fn elem_int<T: Into<u32>>(tag: dicom::core::Tag, vr: VR, v: T) -> DataElement<InMemDicomObject> {
    DataElement::new(tag, vr, PrimitiveValue::from(v.into() as u32))
}
fn elem_multi_f64(tag: dicom::core::Tag, values: &[f64]) -> DataElement<InMemDicomObject> {
    DataElement::new(
        tag,
        VR::DS,
        PrimitiveValue::from(
            values
                .iter()
                .map(|v| v.to_string())
                .collect::<Vec<_>>()
                .join("\\"),
        ),
    )
}
