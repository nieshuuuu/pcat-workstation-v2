//! Pixel decode: open a DICOM file, decode pixel data with the Modality LUT
//! applied by `dicom-pixeldata`, and clamp to i16 (HU).

use std::path::Path;

use dicom::object::open_file;
use dicom_pixeldata::PixelDecoder;

use crate::dicom_errors::DicomLoadError;

/// Decode one slice's pixel data as HU (Modality LUT applied by the library)
/// and return a flat row-major `Vec<i16>` of length `rows * cols`.
///
/// `_rescale_slope` and `_rescale_intercept` are kept for API parity with the
/// series descriptor but are applied internally by `dicom-pixeldata`'s default
/// Modality LUT; passing them here is a no-op.
pub fn decode_slice_i16(
    path: &Path,
    _rescale_slope: f64,
    _rescale_intercept: f64,
    expected_rows: u32,
    expected_cols: u32,
) -> Result<Vec<i16>, DicomLoadError> {
    let obj = open_file(path).map_err(|e| DicomLoadError::ParseFailed {
        path: path.to_path_buf(),
        reason: e.to_string(),
    })?;
    let decoded = obj
        .decode_pixel_data()
        .map_err(|e| DicomLoadError::DecodeFailed {
            path: path.to_path_buf(),
            reason: e.to_string(),
        })?;

    let rows = decoded.rows();
    let cols = decoded.columns();
    if rows != expected_rows || cols != expected_cols {
        return Err(DicomLoadError::InconsistentDims {
            path: path.to_path_buf(),
            rows_got: rows,
            cols_got: cols,
            rows_want: expected_rows,
            cols_want: expected_cols,
        });
    }

    let ndarr = decoded
        .to_ndarray::<i32>()
        .map_err(|e| DicomLoadError::DecodeFailed {
            path: path.to_path_buf(),
            reason: e.to_string(),
        })?;
    let mut out = Vec::with_capacity((rows * cols) as usize);
    for &v in ndarr.iter() {
        out.push(v.clamp(i16::MIN as i32, i16::MAX as i32) as i16);
    }
    Ok(out)
}
