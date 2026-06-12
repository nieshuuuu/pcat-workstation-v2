//! Header-only DICOM scanning primitives.
//!
//! `read_header` opens a DICOM file, stops parsing at the PixelData tag, and
//! returns a `SliceHeader` with only the tags we care about for indexing and
//! grouping. On an SMB share this transfers ~4 KB per file instead of ~512 KB.

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use dicom::core::Tag;
use dicom::dictionary_std::tags;
use dicom::object::{FileDicomObject, InMemDicomObject, OpenFileOptions};

use crate::dicom_errors::DicomLoadError;

/// DICOM tag for ImageComments (0020,4000) — MonoPlus keV truth per lab finding.
pub const IMAGE_COMMENTS: Tag = Tag(0x0020, 0x4000);

/// Per-file header fields sufficient to (a) group files by series and (b) later
/// load pixel data without re-parsing the header.
#[derive(Debug, Clone)]
pub struct SliceHeader {
    pub path: PathBuf,
    pub series_uid: String,
    pub series_description: String,
    pub image_comments: Option<String>,
    pub instance_number: Option<i32>,
    pub image_position_z: Option<f64>,
    /// Full ImagePositionPatient (0020,0032) in patient LPS mm — `[x, y, z]`.
    /// `None` if the tag is absent or malformed. The legacy `image_position_z`
    /// field above mirrors `image_position_patient[2]` and is kept because slice
    /// sorting is in terms of z only.
    pub image_position_patient: Option<[f64; 3]>,
    pub rows: u32,
    pub cols: u32,
    pub rescale_slope: f64,
    pub rescale_intercept: f64,
    pub pixel_spacing: [f64; 2],
    pub orientation: [f64; 6],
    pub patient_name: String,
    pub study_description: String,
    pub window_center: f64,
    pub window_width: f64,
}

/// Read only the file header up to PixelData. Returns Ok(Some(header)) for valid
/// image-DICOM files, Ok(None) for files that are not DICOM or have no pixel data,
/// and Err for hard I/O errors.
pub fn read_header(path: &Path) -> Result<Option<SliceHeader>, DicomLoadError> {
    let obj = match OpenFileOptions::new()
        .read_until(tags::PIXEL_DATA)
        .open_file(path)
    {
        Ok(o) => o,
        // Not a valid DICOM file — skip silently.
        Err(_) => return Ok(None),
    };

    // If there is no SeriesInstanceUID, this is not an image series instance.
    let series_uid = match read_string(&obj, tags::SERIES_INSTANCE_UID) {
        Some(s) if !s.is_empty() => s,
        _ => return Ok(None),
    };

    let rows = read_u32(&obj, tags::ROWS).unwrap_or(0);
    let cols = read_u32(&obj, tags::COLUMNS).unwrap_or(0);
    if rows == 0 || cols == 0 {
        // Non-image DICOM (SR, KO, PR, etc.) — skip.
        return Ok(None);
    }

    let pixel_spacing = read_multi_f64(&obj, tags::PIXEL_SPACING);
    let orient = read_multi_f64(&obj, tags::IMAGE_ORIENTATION_PATIENT);
    let ipp = read_multi_f64(&obj, tags::IMAGE_POSITION_PATIENT);

    let image_position_patient = if ipp.len() >= 3 {
        Some([ipp[0], ipp[1], ipp[2]])
    } else {
        None
    };

    Ok(Some(SliceHeader {
        path: path.to_path_buf(),
        series_uid,
        series_description: read_string(&obj, tags::SERIES_DESCRIPTION).unwrap_or_default(),
        image_comments: read_string(&obj, IMAGE_COMMENTS),
        instance_number: read_i32(&obj, tags::INSTANCE_NUMBER),
        image_position_z: image_position_patient.map(|p| p[2]),
        image_position_patient,
        rows,
        cols,
        rescale_slope: read_f64(&obj, tags::RESCALE_SLOPE).unwrap_or(1.0),
        rescale_intercept: read_f64(&obj, tags::RESCALE_INTERCEPT).unwrap_or(0.0),
        pixel_spacing: if pixel_spacing.len() >= 2 {
            [pixel_spacing[0], pixel_spacing[1]]
        } else {
            [1.0, 1.0]
        },
        orientation: if orient.len() >= 6 {
            [orient[0], orient[1], orient[2], orient[3], orient[4], orient[5]]
        } else {
            [1.0, 0.0, 0.0, 0.0, 1.0, 0.0]
        },
        patient_name: read_string(&obj, tags::PATIENT_NAME).unwrap_or_default(),
        study_description: read_string(&obj, tags::STUDY_DESCRIPTION).unwrap_or_default(),
        window_center: read_f64(&obj, tags::WINDOW_CENTER).unwrap_or(40.0),
        window_width: read_f64(&obj, tags::WINDOW_WIDTH).unwrap_or(400.0),
    }))
}

fn read_string(obj: &FileDicomObject<InMemDicomObject>, tag: Tag) -> Option<String> {
    obj.element(tag)
        .ok()
        .and_then(|e| e.to_str().ok())
        .map(|s| s.trim().to_string())
}

fn read_f64(obj: &FileDicomObject<InMemDicomObject>, tag: Tag) -> Option<f64> {
    obj.element(tag).ok().and_then(|e| e.to_float64().ok())
}

fn read_multi_f64(obj: &FileDicomObject<InMemDicomObject>, tag: Tag) -> Vec<f64> {
    obj.element(tag)
        .ok()
        .and_then(|e| e.to_multi_float64().ok())
        .unwrap_or_default()
}

fn read_u32(obj: &FileDicomObject<InMemDicomObject>, tag: Tag) -> Option<u32> {
    obj.element(tag).ok().and_then(|e| e.to_int::<u32>().ok())
}

fn read_i32(obj: &FileDicomObject<InMemDicomObject>, tag: Tag) -> Option<i32> {
    obj.element(tag).ok().and_then(|e| e.to_int::<i32>().ok())
}

/// Partition a flat list of slice headers into a map keyed by SeriesInstanceUID,
/// with each group sorted by `image_position_z` (falling back to
/// `instance_number`, then preserving input order for fully-unordered series).
pub fn group_by_series(headers: Vec<SliceHeader>) -> HashMap<String, Vec<SliceHeader>> {
    let mut groups: HashMap<String, Vec<SliceHeader>> = HashMap::new();
    for h in headers {
        groups.entry(h.series_uid.clone()).or_default().push(h);
    }
    for slices in groups.values_mut() {
        // Stable sort preserves original order for tied keys; pairs with the
        // "no key" fallback branch below so unordered series stay in input order.
        slices.sort_by(|a, b| {
            match (a.image_position_z, b.image_position_z) {
                (Some(za), Some(zb)) => za.partial_cmp(&zb).unwrap_or(std::cmp::Ordering::Equal),
                _ => match (a.instance_number, b.instance_number) {
                    (Some(ia), Some(ib)) => ia.cmp(&ib),
                    _ => std::cmp::Ordering::Equal,
                },
            }
        });
    }
    groups
}

use futures::stream::{self, StreamExt};
use tokio::sync::Semaphore;

/// Concurrent header opens. Empirically 32–64 is the SMB sweet spot; 48 is a
/// conservative middle value.
// Concurrent header reads. The scan is I/O-latency bound on network mounts
// (SMB) — each open is a slow round-trip — so we read far more than CPU cores in
// flight to overlap the latency.
const SCAN_CONCURRENCY: usize = 64;

/// Public descriptor for a single series, used by Tauri commands and frontend.
#[derive(Debug, Clone, serde::Serialize)]
pub struct SeriesDescriptor {
    pub uid: String,
    pub description: String,
    pub image_comments: Option<String>,
    pub rows: u32,
    pub cols: u32,
    pub num_slices: usize,
    pub pixel_spacing: [f64; 2],
    pub slice_spacing: f64,
    pub orientation: [f64; 6],
    pub rescale_slope: f64,
    pub rescale_intercept: f64,
    pub window_center: f64,
    pub window_width: f64,
    pub patient_name: String,
    pub study_description: String,
    pub file_paths: Vec<PathBuf>,
    /// Per-slice z position (parallel to `file_paths`).
    pub slice_positions_z: Vec<f64>,
    /// Full `ImagePositionPatient` of the first slice in patient LPS mm,
    /// `[x, y, z]`. Defaults to `[0, 0, 0]` when the tag is absent on every
    /// slice. Used by callers to anchor the volume in patient space — the
    /// previous practice of dropping x and y would silently miscompute voxel
    /// indices for any acquisition that wasn't centered at isocenter.
    pub image_position_patient: [f64; 3],
}

/// Walk a folder (non-recursively) and return one SeriesDescriptor per
/// SeriesInstanceUID found. Header-only; does not touch pixel data.
pub async fn scan_series(dir: &Path) -> Result<Vec<SeriesDescriptor>, DicomLoadError> {
    let mut paths: Vec<PathBuf> = Vec::new();
    let mut entries = tokio::fs::read_dir(dir).await?;
    while let Some(entry) = entries.next_entry().await? {
        let file_type = entry.file_type().await?;
        if file_type.is_file() {
            paths.push(entry.path());
        }
    }

    let scanned = paths.len();
    if scanned == 0 {
        return Err(DicomLoadError::NoDicoms { scanned: 0, skipped: 0 });
    }

    let sem = std::sync::Arc::new(Semaphore::new(SCAN_CONCURRENCY));
    let headers: Vec<Option<SliceHeader>> = stream::iter(paths.into_iter().map(|p| {
        let sem = sem.clone();
        async move {
            let _permit = sem.acquire_owned().await.expect("semaphore closed");
            tokio::task::spawn_blocking(move || read_header(&p).ok().flatten())
                .await
                .ok()
                .flatten()
        }
    }))
    .buffer_unordered(SCAN_CONCURRENCY)
    .collect()
    .await;

    let skipped = headers.iter().filter(|h| h.is_none()).count();
    let valid: Vec<SliceHeader> = headers.into_iter().flatten().collect();
    if valid.is_empty() {
        return Err(DicomLoadError::NoDicoms { scanned, skipped });
    }

    let groups = group_by_series(valid);
    let mut descriptors: Vec<SeriesDescriptor> = groups
        .into_iter()
        .map(|(uid, slices)| descriptor_from_slices(uid, slices))
        .collect();
    descriptors.sort_by(|a, b| a.uid.cmp(&b.uid));
    Ok(descriptors)
}

fn descriptor_from_slices(uid: String, slices: Vec<SliceHeader>) -> SeriesDescriptor {
    let first = &slices[0];
    let rows = first.rows;
    let cols = first.cols;
    let pixel_spacing = first.pixel_spacing;
    let orientation = first.orientation;
    let rescale_slope = first.rescale_slope;
    let rescale_intercept = first.rescale_intercept;
    let window_center = first.window_center;
    let window_width = first.window_width;
    let patient_name = first.patient_name.clone();
    let study_description = first.study_description.clone();
    let description = first.series_description.clone();
    let image_comments = first.image_comments.clone();

    let file_paths: Vec<PathBuf> = slices.iter().map(|h| h.path.clone()).collect();
    let slice_positions_z: Vec<f64> = slices
        .iter()
        .enumerate()
        .map(|(i, h)| h.image_position_z.unwrap_or(i as f64))
        .collect();

    // Infer slice spacing from the first two positions (or default 1.0).
    let slice_spacing = if slice_positions_z.len() >= 2 {
        (slice_positions_z[1] - slice_positions_z[0]).abs().max(1e-6)
    } else {
        1.0
    };

    let image_position_patient = first
        .image_position_patient
        .unwrap_or([0.0, 0.0, slice_positions_z.first().copied().unwrap_or(0.0)]);

    SeriesDescriptor {
        uid,
        description,
        image_comments,
        rows,
        cols,
        num_slices: slices.len(),
        pixel_spacing,
        slice_spacing,
        orientation,
        rescale_slope,
        rescale_intercept,
        window_center,
        window_width,
        patient_name,
        study_description,
        file_paths,
        slice_positions_z,
        image_position_patient,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn non_dicom_file_returns_none() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("notdicom.txt");
        let mut f = std::fs::File::create(&path).unwrap();
        f.write_all(b"This is not a DICOM file.").unwrap();
        drop(f);

        let result = read_header(&path).unwrap();
        assert!(result.is_none());
    }

    #[test]
    fn empty_file_returns_none() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("empty.dcm");
        std::fs::File::create(&path).unwrap();
        let result = read_header(&path).unwrap();
        assert!(result.is_none());
    }

    // Note: reading an actual DICOM file is covered by integration tests
    // in tests/dicom_loader_integration.rs, after the fixture generator is
    // landed in Task 4.

    fn h(uid: &str, z: Option<f64>, inst: Option<i32>, path: &str) -> SliceHeader {
        SliceHeader {
            path: std::path::PathBuf::from(path),
            series_uid: uid.to_string(),
            series_description: String::new(),
            image_comments: None,
            instance_number: inst,
            image_position_z: z,
            image_position_patient: z.map(|zz| [0.0, 0.0, zz]),
            rows: 64,
            cols: 64,
            rescale_slope: 1.0,
            rescale_intercept: -1024.0,
            pixel_spacing: [1.0, 1.0],
            orientation: [1.0, 0.0, 0.0, 0.0, 1.0, 0.0],
            patient_name: String::new(),
            study_description: String::new(),
            window_center: 40.0,
            window_width: 400.0,
        }
    }

    #[test]
    fn groups_by_series_and_sorts_by_ipp() {
        let headers = vec![
            h("A", Some(2.0), None, "a2"),
            h("B", Some(0.0), None, "b0"),
            h("A", Some(0.0), None, "a0"),
            h("A", Some(1.0), None, "a1"),
            h("B", Some(1.0), None, "b1"),
        ];
        let groups = group_by_series(headers);
        assert_eq!(groups.len(), 2);

        let a = groups.get("A").unwrap();
        assert_eq!(
            a.iter().map(|h| h.path.to_string_lossy().into_owned()).collect::<Vec<_>>(),
            vec!["a0", "a1", "a2"],
        );

        let b = groups.get("B").unwrap();
        assert_eq!(
            b.iter().map(|h| h.path.to_string_lossy().into_owned()).collect::<Vec<_>>(),
            vec!["b0", "b1"],
        );
    }

    #[test]
    fn falls_back_to_instance_number_when_ipp_missing() {
        let headers = vec![
            h("A", None, Some(3), "c"),
            h("A", None, Some(1), "a"),
            h("A", None, Some(2), "b"),
        ];
        let groups = group_by_series(headers);
        let a = groups.get("A").unwrap();
        assert_eq!(
            a.iter().map(|h| h.path.to_string_lossy().into_owned()).collect::<Vec<_>>(),
            vec!["a", "b", "c"],
        );
    }

    #[test]
    fn headers_with_no_sort_key_preserve_input_order() {
        let headers = vec![
            h("A", None, None, "first"),
            h("A", None, None, "second"),
        ];
        let groups = group_by_series(headers);
        let a = groups.get("A").unwrap();
        assert_eq!(
            a.iter().map(|h| h.path.to_string_lossy().into_owned()).collect::<Vec<_>>(),
            vec!["first", "second"],
        );
    }
}
