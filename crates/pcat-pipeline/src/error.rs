use thiserror::Error;

#[allow(dead_code)]
#[derive(Debug, Error)]
pub enum AppError {
    #[error("DICOM error: {0}")]
    Dicom(String),
    #[error("No volume loaded")]
    NoVolume,
    #[error("Invalid argument: {0}")]
    InvalidArg(String),
    #[error("Pipeline error: {0}")]
    Pipeline(String),
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),
}

impl serde::Serialize for AppError {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        serializer.serialize_str(&self.to_string())
    }
}
