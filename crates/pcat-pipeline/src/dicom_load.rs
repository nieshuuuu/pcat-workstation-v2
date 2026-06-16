//! Parallel pixel decode + volume assembly.
//!
//! Given a folder and a series UID, walks headers (to locate files), then
//! decodes pixel data in parallel (rayon) and returns a densely packed i16
//! volume in z-major order.

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};

use rayon::prelude::*;

use crate::dicom_errors::DicomLoadError;
use crate::dicom_scan::{descriptors_from_headers, read_slice_with_pixels, SeriesDescriptor, SliceHeader};

/// 4 GB soft limit (conservative — covers 1000-slice 512² i16 at 1.5 GB).
const VOLUME_SIZE_LIMIT_MB: usize = 4096;

/// Concurrent file reads during pixel decode. The decode is I/O-LATENCY bound on
/// network mounts (SMB), where each file open is a slow round-trip; the default
/// rayon pool (~num_cores) leaves the link mostly idle. Far more threads overlap
/// the per-file latency. The CPU decode is ~0.5 ms/slice, so on local disk the
/// extra threads finish near-instantly (negligible overhead).
const READ_CONCURRENCY: usize = 64;

/// Metadata subset that travels with a loaded volume's pixel bytes.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct VolumeMetadata {
    pub series_uid: String,
    pub series_description: String,
    pub image_comments: Option<String>,
    pub rows: u32,
    pub cols: u32,
    pub num_slices: usize,
    pub pixel_spacing: [f64; 2],
    pub slice_spacing: f64,
    pub orientation: [f64; 6],
    pub window_center: f64,
    pub window_width: f64,
    pub patient_name: String,
    pub study_description: String,
    pub slice_positions_z: Vec<f64>,
    /// `ImagePositionPatient` of the first slice in patient LPS mm, `[x, y, z]`.
    /// Required for correct voxel-index conversion in any sampler that consumes
    /// patient-space coordinates (CPR, ROI, radial-angular).
    pub image_position_patient: [f64; 3],
}

impl From<&SeriesDescriptor> for VolumeMetadata {
    fn from(d: &SeriesDescriptor) -> Self {
        Self {
            series_uid: d.uid.clone(),
            series_description: d.description.clone(),
            image_comments: d.image_comments.clone(),
            rows: d.rows,
            cols: d.cols,
            num_slices: d.num_slices,
            pixel_spacing: d.pixel_spacing,
            slice_spacing: d.slice_spacing,
            orientation: d.orientation,
            window_center: d.window_center,
            window_width: d.window_width,
            patient_name: d.patient_name.clone(),
            study_description: d.study_description.clone(),
            slice_positions_z: d.slice_positions_z.clone(),
            image_position_patient: d.image_position_patient,
        }
    }
}

#[derive(Debug)]
pub struct LoadedVolume {
    pub metadata: VolumeMetadata,
    /// Tightly packed i16, z-major order: `voxels_i16[z * rows * cols + r * cols + c]`.
    pub voxels_i16: Vec<i16>,
}

/// Return Err if a planned volume would exceed the size limit.
pub fn check_volume_size_mb(requested_mb: usize) -> Result<(), DicomLoadError> {
    if requested_mb > VOLUME_SIZE_LIMIT_MB {
        Err(DicomLoadError::VolumeTooLarge {
            requested_mb,
            limit_mb: VOLUME_SIZE_LIMIT_MB,
        })
    } else {
        Ok(())
    }
}

/// Load a single series by UID. Runs `scan_series` to find the descriptor, then
/// rayon-parallel decodes all slices.
///
/// `on_progress`, if provided, is called as `(done, total)` periodically during
/// pixel decode. It will be called once with `(0, total)` before decode begins,
/// then at most every `total / 50` slices, and once more at completion.
pub async fn load_series(
    dir: &Path,
    uid: &str,
    on_progress: Option<Box<dyn Fn(usize, usize) + Send + Sync>>,
) -> Result<LoadedVolume, DicomLoadError> {
    // Single directory round-trip to enumerate candidate files. Header AND
    // pixels then come from ONE open per file below (see `read_slice_with_pixels`),
    // instead of the old scan-then-decode which opened every file twice — the
    // second open was the dominant cost on SMB.
    let mut file_paths: Vec<PathBuf> = Vec::new();
    let mut entries = tokio::fs::read_dir(dir).await?;
    while let Some(entry) = entries.next_entry().await? {
        if entry.file_type().await?.is_file() {
            file_paths.push(entry.path());
        }
    }
    if file_paths.is_empty() {
        return Err(DicomLoadError::NoDicoms { scanned: 0, skipped: 0 });
    }
    file_paths.sort();
    let total_files = file_paths.len();

    // Emit initial progress (`total_files` ≈ slice count for a single-series
    // folder; the exact count is known after grouping below).
    if let Some(ref cb) = on_progress {
        cb(0, total_files);
    }
    let progress = on_progress.map(Arc::new);
    let uid_owned = uid.to_string();
    let dir_owned = dir.to_path_buf();

    let (metadata, voxels) = tokio::task::spawn_blocking(
        move || -> Result<(VolumeMetadata, Vec<i16>), DicomLoadError> {
            let counter = Arc::new(AtomicUsize::new(0));
            let step = (total_files / 50).max(1);

            // One open per file → (header, pixels). High-concurrency pool so the
            // SMB round-trips overlap. Propagates the first hard decode error.
            let run = || -> Result<Vec<(SliceHeader, Vec<i16>)>, DicomLoadError> {
                file_paths
                    .par_iter()
                    .map(|p| {
                        let res = read_slice_with_pixels(p);
                        let done = counter.fetch_add(1, Ordering::Relaxed) + 1;
                        if let Some(ref cb) = progress {
                            if done % step == 0 || done == total_files {
                                cb(done, total_files);
                            }
                        }
                        res
                    })
                    .collect::<Result<Vec<Option<(SliceHeader, Vec<i16>)>>, DicomLoadError>>()
                    .map(|v| v.into_iter().flatten().collect())
            };
            let slices: Vec<(SliceHeader, Vec<i16>)> = match rayon::ThreadPoolBuilder::new()
                .num_threads(READ_CONCURRENCY)
                .build()
            {
                Ok(pool) => pool.install(run),
                Err(_) => run(),
            }?;

            if slices.is_empty() {
                return Err(DicomLoadError::NoDicoms {
                    scanned: total_files,
                    skipped: total_files,
                });
            }

            // Reuse scan's grouping + z-ordering on the headers we just read, so
            // the reconstructed volume is byte-identical to the scan-then-decode
            // path. Pixels are looked up by path in the descriptor's z-order.
            let mut pixels_by_path: HashMap<PathBuf, Vec<i16>> =
                HashMap::with_capacity(slices.len());
            let mut headers: Vec<SliceHeader> = Vec::with_capacity(slices.len());
            for (h, px) in slices {
                pixels_by_path.insert(h.path.clone(), px);
                headers.push(h);
            }

            let desc = descriptors_from_headers(headers)
                .into_iter()
                .find(|d| d.uid == uid_owned)
                .ok_or(DicomLoadError::SeriesNotFound { uid: uid_owned })?;

            let slice_len = (desc.rows as usize) * (desc.cols as usize);
            let total_voxels = slice_len * desc.num_slices;
            let total_bytes_mb = (total_voxels * 2) / (1024 * 1024);
            check_volume_size_mb(total_bytes_mb)?;

            let mut out = vec![0i16; total_voxels];
            for (z, p) in desc.file_paths.iter().enumerate() {
                let px = pixels_by_path
                    .remove(p)
                    .ok_or_else(|| DicomLoadError::ParseFailed {
                        path: p.clone(),
                        reason: "slice pixels missing after single-pass decode".to_string(),
                    })?;
                if px.len() != slice_len {
                    return Err(DicomLoadError::InconsistentDims {
                        path: p.clone(),
                        rows_got: (px.len() / (desc.cols.max(1) as usize)) as u32,
                        cols_got: desc.cols,
                        rows_want: desc.rows,
                        cols_want: desc.cols,
                    });
                }
                out[z * slice_len..(z + 1) * slice_len].copy_from_slice(&px);
            }

            Ok((VolumeMetadata::from(&desc), out))
        },
    )
    .await
    .map_err(|e| DicomLoadError::ParseFailed {
        path: dir_owned,
        reason: format!("decode task panicked: {e}"),
    })??;

    Ok(LoadedVolume {
        metadata,
        voxels_i16: voxels,
    })
}
