use std::path::PathBuf;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum DicomLoadError {
    #[error("folder not readable: {0}")]
    IoError(#[from] std::io::Error),

    #[error("no DICOM files found ({scanned} scanned, {skipped} failed header parse)")]
    NoDicoms { scanned: usize, skipped: usize },

    #[error("series {uid} not found in folder")]
    SeriesNotFound { uid: String },

    #[error("slice {path:?} has dims ({rows_got}x{cols_got}), expected ({rows_want}x{cols_want})")]
    InconsistentDims {
        path: PathBuf,
        rows_got: u32,
        cols_got: u32,
        rows_want: u32,
        cols_want: u32,
    },

    #[error("unsupported transfer syntax: {0}")]
    UnsupportedTransferSyntax(String),

    #[error("volume too large: {requested_mb} MB exceeds {limit_mb} MB limit")]
    VolumeTooLarge { requested_mb: usize, limit_mb: usize },

    #[error("pixel decode failed at {path:?}: {reason}")]
    DecodeFailed { path: PathBuf, reason: String },

    #[error("dicom parse error at {path:?}: {reason}")]
    ParseFailed { path: PathBuf, reason: String },
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn display_series_not_found() {
        let e = DicomLoadError::SeriesNotFound { uid: "1.2.3".into() };
        assert_eq!(format!("{e}"), "series 1.2.3 not found in folder");
    }

    #[test]
    fn display_no_dicoms() {
        let e = DicomLoadError::NoDicoms { scanned: 42, skipped: 5 };
        assert!(format!("{e}").contains("42 scanned"));
        assert!(format!("{e}").contains("5 failed"));
    }

    #[test]
    fn io_error_converts() {
        let io = std::io::Error::new(std::io::ErrorKind::NotFound, "x");
        let err: DicomLoadError = io.into();
        assert!(matches!(err, DicomLoadError::IoError(_)));
    }
}
