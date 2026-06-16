use std::collections::HashMap;
use std::path::Path;
use std::sync::Arc;

use dicom::core::Tag;
use dicom::dictionary_std::tags;
use dicom::object::{open_file, FileDicomObject, InMemDicomObject};
use dicom_pixeldata::PixelDecoder;
use ndarray::Array3;
use walkdir::WalkDir;

use crate::error::AppError;
use crate::types::LoadedVolume;

// ---------------------------------------------------------------------------
// ImageComments tag — not in the standard dictionary
// ---------------------------------------------------------------------------

/// DICOM tag for ImageComments (0020,4000).
const IMAGE_COMMENTS: Tag = Tag(0x0020, 0x4000);

// ---------------------------------------------------------------------------
// Internal types
// ---------------------------------------------------------------------------

/// Per-slice metadata extracted from a DICOM file.
struct SliceInfo {
    obj: FileDicomObject<InMemDicomObject>,
    z_position: f64,
    rescale_slope: f64,
    rescale_intercept: f64,
}

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

/// Metadata about a single DICOM series discovered during a directory scan.
#[derive(Debug, Clone, serde::Serialize)]
pub struct SeriesInfo {
    pub series_uid: String,
    pub description: String,
    pub num_slices: usize,
    pub kev_label: Option<f64>,
}

/// A pair of co-registered CT volumes acquired at two different X-ray energies.
pub struct DualEnergyVolume {
    pub low: Arc<Array3<f32>>,
    pub high: Arc<Array3<f32>>,
    pub low_energy_kev: f64,
    pub high_energy_kev: f64,
    pub spacing: [f64; 3],
    pub origin: [f64; 3],
    pub direction: [f64; 9],
    pub patient_name: String,
    pub study_description: String,
}

// ---------------------------------------------------------------------------
// Tag helper functions
// ---------------------------------------------------------------------------

/// Read a multi-value f64 tag, returning an empty vec on failure.
fn read_multi_f64(obj: &FileDicomObject<InMemDicomObject>, tag: Tag) -> Vec<f64> {
    obj.element(tag)
        .ok()
        .and_then(|e| e.to_multi_float64().ok())
        .unwrap_or_default()
}

/// Read a single f64 tag with a fallback default.
fn read_f64(obj: &FileDicomObject<InMemDicomObject>, tag: Tag, default: f64) -> f64 {
    obj.element(tag)
        .ok()
        .and_then(|e| e.to_float64().ok())
        .unwrap_or(default)
}

/// Read a string tag with a fallback default.
fn read_string(obj: &FileDicomObject<InMemDicomObject>, tag: Tag, default: &str) -> String {
    obj.element(tag)
        .ok()
        .and_then(|e| e.to_str().ok())
        .map(|s| s.trim().to_string())
        .unwrap_or_else(|| default.to_string())
}

/// Read a u16 tag with a fallback default.
fn read_u16(obj: &FileDicomObject<InMemDicomObject>, tag: Tag, default: u16) -> u16 {
    obj.element(tag)
        .ok()
        .and_then(|e| e.to_int::<u16>().ok())
        .unwrap_or(default)
}

// ---------------------------------------------------------------------------
// keV parsing
// ---------------------------------------------------------------------------

/// Try to extract a keV energy value from DICOM metadata strings.
///
/// Strategy:
///   1. Check ImageComments for pattern `E = <number> keV` (case-insensitive).
///   2. Check SeriesDescription for a number immediately followed by `keV`
///      (e.g. "MonoPlus 70 keV", "70keV").
fn parse_kev(image_comments: &str, series_description: &str) -> Option<f64> {
    // Strategy 1: "E = 150 keV" in ImageComments
    if let Some(val) = parse_kev_from_image_comments(image_comments) {
        return Some(val);
    }
    // Strategy 2: "<number> keV" or "<number>keV" in SeriesDescription
    parse_kev_from_description(series_description)
}

/// Parse `E = <number> keV` from ImageComments.
fn parse_kev_from_image_comments(s: &str) -> Option<f64> {
    let lower = s.to_ascii_lowercase();
    // Look for "e" followed by optional whitespace, "=", optional whitespace, a number,
    // optional whitespace, and "kev".
    let re_like = |s: &str| -> Option<f64> {
        let bytes = s.as_bytes();
        let len = bytes.len();
        let mut i = 0;
        while i < len {
            // Find 'e'
            if bytes[i] == b'e' {
                let mut j = i + 1;
                // skip whitespace
                while j < len && bytes[j] == b' ' {
                    j += 1;
                }
                // expect '='
                if j < len && bytes[j] == b'=' {
                    j += 1;
                    // skip whitespace
                    while j < len && bytes[j] == b' ' {
                        j += 1;
                    }
                    // parse number
                    let num_start = j;
                    while j < len && (bytes[j].is_ascii_digit() || bytes[j] == b'.') {
                        j += 1;
                    }
                    if j > num_start {
                        let num_str = &s[num_start..j];
                        if let Ok(val) = num_str.parse::<f64>() {
                            // skip whitespace
                            while j < len && bytes[j] == b' ' {
                                j += 1;
                            }
                            // check for "kev"
                            if j + 3 <= len && &s[j..j + 3] == "kev" {
                                return Some(val);
                            }
                        }
                    }
                }
            }
            i += 1;
        }
        None
    };
    re_like(&lower)
}

/// Parse `<number> keV` or `<number>keV` from SeriesDescription.
fn parse_kev_from_description(s: &str) -> Option<f64> {
    let lower = s.to_ascii_lowercase();
    // Find "kev" and look backwards for a number
    let mut search_from = 0;
    while let Some(pos) = lower[search_from..].find("kev") {
        let abs_pos = search_from + pos;
        // Walk backwards from "kev", skipping whitespace, then parse digits/dot
        let before = &lower[..abs_pos];
        let trimmed = before.trim_end();
        if !trimmed.is_empty() {
            // Find the start of the number at the end of trimmed
            let num_end = trimmed.len();
            let mut num_start = num_end;
            let tbytes = trimmed.as_bytes();
            while num_start > 0
                && (tbytes[num_start - 1].is_ascii_digit() || tbytes[num_start - 1] == b'.')
            {
                num_start -= 1;
            }
            if num_start < num_end {
                if let Ok(val) = trimmed[num_start..num_end].parse::<f64>() {
                    return Some(val);
                }
            }
        }
        search_from = abs_pos + 3;
    }
    None
}

// ---------------------------------------------------------------------------
// Pixel data extraction
// ---------------------------------------------------------------------------

/// Decode pixel data from a DICOM object, apply rescale slope/intercept to
/// convert to Hounsfield Units, and clamp outliers.
///
/// Returns a flat `Vec<f32>` of length `rows * cols`.
fn decode_slice_hu(
    obj: &FileDicomObject<InMemDicomObject>,
    rescale_slope: f64,
    rescale_intercept: f64,
    expected_len: usize,
) -> Result<Vec<f32>, AppError> {
    let pixel_data = obj
        .decode_pixel_data()
        .map_err(|e| AppError::Dicom(format!("pixel decode error: {e}")))?;

    let pixel_repr = read_u16(obj, tags::PIXEL_REPRESENTATION, 0);

    let hu_values: Vec<f32> = if pixel_repr == 1 {
        // Signed 16-bit
        let raw: Vec<i16> = pixel_data
            .data()
            .chunks_exact(2)
            .map(|c| i16::from_le_bytes([c[0], c[1]]))
            .collect();
        raw.iter()
            .map(|&v| {
                let hu = v as f64 * rescale_slope + rescale_intercept;
                clamp_hu(hu) as f32
            })
            .collect()
    } else {
        // Unsigned 16-bit
        let raw: Vec<u16> = pixel_data
            .data()
            .chunks_exact(2)
            .map(|c| u16::from_le_bytes([c[0], c[1]]))
            .collect();
        raw.iter()
            .map(|&v| {
                let hu = v as f64 * rescale_slope + rescale_intercept;
                clamp_hu(hu) as f32
            })
            .collect()
    };

    if hu_values.len() != expected_len {
        return Err(AppError::Dicom(format!(
            "pixel count mismatch: got {} expected {}",
            hu_values.len(),
            expected_len
        )));
    }

    Ok(hu_values)
}

/// Clamp HU values: anything <= -8192 becomes -1024, anything > 3095 becomes 3095.
#[inline]
fn clamp_hu(hu: f64) -> f64 {
    if hu <= -8192.0 {
        -1024.0
    } else if hu > 3095.0 {
        3095.0
    } else {
        hu
    }
}

// ---------------------------------------------------------------------------
// Direction matrix from ImageOrientationPatient
// ---------------------------------------------------------------------------

/// Build a 3x3 direction cosine matrix (row-major, 9 elements) from
/// ImageOrientationPatient (6 elements: row_x, row_y, row_z, col_x, col_y, col_z).
/// The third row is the cross product of the row and column direction vectors.
fn build_direction_matrix(iop: &[f64]) -> [f64; 9] {
    if iop.len() < 6 {
        // Fallback: identity matrix
        return [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
    }
    let (rx, ry, rz) = (iop[0], iop[1], iop[2]);
    let (cx, cy, cz) = (iop[3], iop[4], iop[5]);

    // Cross product: normal = row x col
    let nx = ry * cz - rz * cy;
    let ny = rz * cx - rx * cz;
    let nz = rx * cy - ry * cx;

    [rx, ry, rz, cx, cy, cz, nx, ny, nz]
}

// ---------------------------------------------------------------------------
// Core volume-building helper (shared by load_dicom_directory and load_series_by_uid)
// ---------------------------------------------------------------------------

/// Build a `LoadedVolume` from a pre-collected vector of `SliceInfo`.
///
/// The caller is responsible for filtering which DICOM files to include;
/// this function handles sorting, metadata extraction, pixel decoding, and
/// array construction.
fn build_volume_from_slices(mut slices: Vec<SliceInfo>) -> Result<LoadedVolume, AppError> {
    if slices.is_empty() {
        return Err(AppError::Dicom(
            "no valid DICOM image slices found".to_string(),
        ));
    }

    // Sort by Z position ascending.
    slices.sort_by(|a, b| {
        a.z_position
            .partial_cmp(&b.z_position)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    // Extract metadata from the first slice (representative).
    let first = &slices[0].obj;

    let rows = read_u16(first, tags::ROWS, 0) as usize;
    let cols = read_u16(first, tags::COLUMNS, 0) as usize;
    if rows == 0 || cols == 0 {
        return Err(AppError::Dicom("invalid image dimensions".to_string()));
    }

    let pixel_spacing = read_multi_f64(first, tags::PIXEL_SPACING);
    let (sy, sx) = if pixel_spacing.len() >= 2 {
        (pixel_spacing[0], pixel_spacing[1])
    } else {
        (1.0, 1.0)
    };

    let sz = if slices.len() > 1 {
        (slices[1].z_position - slices[0].z_position).abs()
    } else {
        1.0
    };

    let ipp_first = read_multi_f64(first, tags::IMAGE_POSITION_PATIENT);
    let origin = if ipp_first.len() >= 3 {
        [ipp_first[2], ipp_first[1], ipp_first[0]] // [oz, oy, ox]
    } else {
        [0.0, 0.0, 0.0]
    };

    let iop = read_multi_f64(first, tags::IMAGE_ORIENTATION_PATIENT);
    let direction = build_direction_matrix(&iop);

    let window_center = read_f64(first, tags::WINDOW_CENTER, 40.0);
    let window_width = read_f64(first, tags::WINDOW_WIDTH, 400.0);
    let patient_name = read_string(first, tags::PATIENT_NAME, "");
    let study_description = read_string(first, tags::STUDY_DESCRIPTION, "");

    // Decode pixel data for every slice and stack into Array3.
    let n_slices = slices.len();
    let slice_len = rows * cols;
    let mut volume_data: Vec<f32> = Vec::with_capacity(n_slices * slice_len);

    for (i, si) in slices.iter().enumerate() {
        let hu = decode_slice_hu(&si.obj, si.rescale_slope, si.rescale_intercept, slice_len)
            .map_err(|e| AppError::Dicom(format!("slice {i}: {e}")))?;
        volume_data.extend_from_slice(&hu);
    }

    let data = Array3::from_shape_vec((n_slices, rows, cols), volume_data)
        .map_err(|e| AppError::Dicom(format!("failed to build 3D array: {e}")))?;

    Ok(LoadedVolume {
        data: Arc::new(data),
        spacing: [sz, sy, sx],
        origin,
        direction,
        window_center,
        window_width,
        patient_name,
        study_description,
    })
}

/// Scan a directory for DICOM files and collect `SliceInfo` entries,
/// optionally filtered by a SeriesInstanceUID.
fn collect_slices(
    dir: &Path,
    series_uid_filter: Option<&str>,
) -> Vec<SliceInfo> {
    let mut slices: Vec<SliceInfo> = Vec::new();

    for entry in WalkDir::new(dir).into_iter().filter_map(|e| e.ok()) {
        if !entry.file_type().is_file() {
            continue;
        }
        let obj = match open_file(entry.path()) {
            Ok(o) => o,
            Err(_) => continue, // skip non-DICOM files silently
        };

        // If a UID filter is provided, skip files that don't match.
        if let Some(uid) = series_uid_filter {
            let file_uid = read_string(&obj, tags::SERIES_INSTANCE_UID, "");
            if file_uid != uid {
                continue;
            }
        }

        let ipp = read_multi_f64(&obj, tags::IMAGE_POSITION_PATIENT);
        if ipp.len() < 3 {
            continue; // not a valid image slice
        }

        let rescale_slope = read_f64(&obj, tags::RESCALE_SLOPE, 1.0);
        let rescale_intercept = read_f64(&obj, tags::RESCALE_INTERCEPT, 0.0);

        slices.push(SliceInfo {
            obj,
            z_position: ipp[2],
            rescale_slope,
            rescale_intercept,
        });
    }

    slices
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Load all DICOM slices from `dir`, sort by Z position, apply HU rescaling,
/// and return a `LoadedVolume`.
pub fn load_dicom_directory(dir: &Path) -> Result<LoadedVolume, AppError> {
    let slices = collect_slices(dir, None);
    build_volume_from_slices(slices)
}

/// Load only the DICOM slices matching a specific SeriesInstanceUID from `dir`.
fn load_series_by_uid(dir: &Path, series_uid: &str) -> Result<LoadedVolume, AppError> {
    let slices = collect_slices(dir, Some(series_uid));
    if slices.is_empty() {
        return Err(AppError::Dicom(format!(
            "no DICOM slices found for SeriesInstanceUID '{series_uid}'"
        )));
    }
    build_volume_from_slices(slices)
}

/// Scan a DICOM directory, group files by SeriesInstanceUID, and return
/// metadata per series.
///
/// This reads only the header of each file (no pixel data) and is therefore
/// fast even for directories with many slices.
pub fn scan_dicom_series(dir: &Path) -> Result<Vec<SeriesInfo>, AppError> {
    // Map: series_uid -> (description, image_comments, count)
    let mut map: HashMap<String, (String, String, usize)> = HashMap::new();

    for entry in WalkDir::new(dir).into_iter().filter_map(|e| e.ok()) {
        if !entry.file_type().is_file() {
            continue;
        }
        let obj = match open_file(entry.path()) {
            Ok(o) => o,
            Err(_) => continue,
        };

        let series_uid = read_string(&obj, tags::SERIES_INSTANCE_UID, "");
        if series_uid.is_empty() {
            continue;
        }

        let entry_val = map.entry(series_uid).or_insert_with(|| {
            let desc = read_string(&obj, tags::SERIES_DESCRIPTION, "");
            let comments = read_string(&obj, IMAGE_COMMENTS, "");
            (desc, comments, 0)
        });
        entry_val.2 += 1;
    }

    let mut result: Vec<SeriesInfo> = map
        .into_iter()
        .map(|(uid, (desc, comments, count))| {
            let kev_label = parse_kev(&comments, &desc);
            SeriesInfo {
                series_uid: uid,
                description: desc,
                num_slices: count,
                kev_label,
            }
        })
        .collect();

    result.sort_by(|a, b| a.description.cmp(&b.description));

    Ok(result)
}

/// Load two DICOM series (low and high energy) from the same directory.
///
/// The two series must have matching geometry (same rows, cols, num_slices,
/// spacing within 0.01 mm tolerance).
pub fn load_dual_energy(
    dir: &Path,
    low_series_uid: &str,
    high_series_uid: &str,
    low_kev: f64,
    high_kev: f64,
) -> Result<DualEnergyVolume, AppError> {
    let low_vol = load_series_by_uid(dir, low_series_uid)?;
    let high_vol = load_series_by_uid(dir, high_series_uid)?;

    // Verify matching geometry: dimensions
    let low_shape = low_vol.data.shape();
    let high_shape = high_vol.data.shape();
    if low_shape != high_shape {
        return Err(AppError::Dicom(format!(
            "geometry mismatch: low volume shape {:?} != high volume shape {:?}",
            low_shape, high_shape
        )));
    }

    // Verify matching geometry: spacing (within 0.01 mm tolerance)
    const TOL: f64 = 0.01;
    for i in 0..3 {
        if (low_vol.spacing[i] - high_vol.spacing[i]).abs() > TOL {
            return Err(AppError::Dicom(format!(
                "spacing mismatch on axis {i}: low={:.4} mm vs high={:.4} mm (tolerance {TOL} mm)",
                low_vol.spacing[i], high_vol.spacing[i]
            )));
        }
    }

    Ok(DualEnergyVolume {
        low: low_vol.data,
        high: high_vol.data,
        low_energy_kev: low_kev,
        high_energy_kev: high_kev,
        spacing: low_vol.spacing,
        origin: low_vol.origin,
        direction: low_vol.direction,
        patient_name: low_vol.patient_name,
        study_description: low_vol.study_description,
    })
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_kev_image_comments() {
        assert_eq!(parse_kev("E = 150 keV", ""), Some(150.0));
        assert_eq!(parse_kev("E=70keV", ""), Some(70.0));
        assert_eq!(parse_kev("e = 40 kev", ""), Some(40.0));
        assert_eq!(parse_kev("E =  100  keV", ""), Some(100.0));
    }

    #[test]
    fn test_parse_kev_series_description() {
        assert_eq!(parse_kev("", "MonoPlus 70 keV"), Some(70.0));
        assert_eq!(parse_kev("", "70keV"), Some(70.0));
        assert_eq!(parse_kev("", "MonoPlus 150 keV Soft Tissue"), Some(150.0));
    }

    #[test]
    fn test_parse_kev_image_comments_takes_precedence() {
        // ImageComments should take precedence over SeriesDescription
        assert_eq!(
            parse_kev("E = 150 keV", "MonoPlus 70 keV"),
            Some(150.0)
        );
    }

    #[test]
    fn test_parse_kev_none() {
        assert_eq!(parse_kev("", ""), None);
        assert_eq!(parse_kev("no energy here", "plain description"), None);
    }

    #[test]
    fn test_parse_kev_decimal() {
        assert_eq!(parse_kev("E = 70.5 keV", ""), Some(70.5));
        assert_eq!(parse_kev("", "MonoPlus 70.5 keV"), Some(70.5));
    }
}
