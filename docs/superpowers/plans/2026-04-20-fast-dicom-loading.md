# Fast DICOM Loading Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the current serial-IO, double-walk, 384-IPC-roundtrip DICOM loader with a stateless three-stage pipeline (parallel async header scan → parallel pixel decode → single bulk `ArrayBuffer` transfer), reducing cold SMB load time from 40–85 s to ≤ 10 s.

**Architecture:** Three stages in `pcat-pipeline`:
1. `scan_series(dir) -> Vec<SeriesDescriptor>` — `tokio::fs::read_dir` + `buffer_unordered(48)` + `dicom-rs OpenFileOptions::read_until(PIXEL_DATA)`.
2. `load_series(dir, uid) -> LoadedVolume` — `rayon::par_iter` over sorted file paths + `decode_pixel_data` + rescale to `i16`.
3. Tauri command returns a framed `ArrayBuffer` via `tauri::ipc::Response` (no JSON encoding).

Frontend uses cornerstone3D's boring `createLocalVolume` — no streaming, no custom loader, no per-slice IPC.

**Tech Stack:** Rust (tokio, rayon, futures, dicom-rs 0.9, thiserror 2, ndarray 0.16), Tauri v2, Svelte 5, cornerstone3D.

**Spec:** `docs/specs/2026-04-20-fast-dicom-loading-design.md`.

---

## File Structure

### New files
- `crates/pcat-pipeline/tests/common/mod.rs` — test helpers (fixture generator).
- `crates/pcat-pipeline/tests/dicom_loader_integration.rs` — integration tests against synthetic fixture.
- `crates/pcat-pipeline/benches/dicom_load_bench.rs` — benchmark harness.
- `tests/e2e/dicom-load.spec.ts` — Playwright end-to-end test.

### Modified files
- `crates/pcat-pipeline/Cargo.toml` — add `tokio`, `futures`, `tempfile` (dev).
- `crates/pcat-pipeline/src/dicom_loader.rs` — full rewrite (retains existing public symbols that are consumed elsewhere; removes internal helpers that are no longer used).
- `src-tauri/src/commands/dicom.rs` — add `scan_series`, `load_series`; keep `open_dicom_dialog`, `get_recent_dicoms`, `list_patients`, `list_series_dirs`; remove old `load_dicom`, `load_dual_energy`, `scan_series` (old signature).
- `src-tauri/src/commands/mod.rs` — remove `volume` module.
- `src-tauri/src/lib.rs` (or wherever commands are registered) — update command registration list.
- `src-tauri/src/state.rs` — simplify `AppState.volume` if it carried ndarray-specific data.
- `src/lib/api.ts` — add typed `scanSeries` (new shape), `loadSeries`; remove `getSlice`, `loadDicom`, `loadDualEnergy`.
- `src/lib/cornerstone/volumeLoader.ts` — rewrite to accept a pre-loaded `ArrayBuffer` + metadata.
- `src/lib/stores/volumeStore.svelte.ts` — update `VolumeMetadata` shape if needed.
- `src/App.svelte` — swap the old `loadDicom` → `getSlice` × N flow for `scanSeries` → `loadSeries` → `loadVolume`. Replace `loadDualEnergy` with two parallel `loadSeries` calls.

### Deleted files
- `src-tauri/src/commands/volume.rs` — entire file (`get_slice` is dead).

### Responsibility boundaries
- `dicom_loader.rs` owns directory enumeration, header parsing, pixel decode, sorting, and grouping. Pure Rust, no Tauri imports.
- `commands/dicom.rs` owns Tauri command signatures, framing, error-to-String conversion. No DICOM logic of its own.
- `api.ts` owns the framed-response parser and TypeScript types.
- `volumeLoader.ts` owns cornerstone3D integration only.
- `App.svelte` owns the user flow (picker → series list → load → display).

---

## Cross-cutting conventions

1. **Each task ends with a commit** using a semantic-commit style: `feat(scope): …`, `test(scope): …`, `refactor(scope): …`, `chore(scope): …`.
2. **Test-first:** every task writes a failing test before any production code, runs it to confirm red, writes the minimal implementation, runs to confirm green, then commits.
3. **Run commands from the repo root** (`/Users/shunie/Developer/PCAT/pcat-workstation-v2`) unless a task specifies otherwise.
4. **Never skip hooks.** If a pre-commit hook fails, fix the root cause — don't `--no-verify`.
5. **Use `cargo test -p pcat-pipeline`** to scope tests to the pipeline crate. Use `cargo test` at workspace root for full workspace.

---

## Task 1: Add async + fixture dependencies to pcat-pipeline

**Files:**
- Modify: `crates/pcat-pipeline/Cargo.toml`

- [ ] **Step 1: Add deps**

Edit `crates/pcat-pipeline/Cargo.toml`. Replace the `[dependencies]` section with:

```toml
[dependencies]
ndarray = { version = "0.16", features = ["serde"] }
nalgebra = "0.33"
dicom = "0.9"
dicom-pixeldata = { version = "0.9", features = ["ndarray"] }
bytemuck = { version = "1", features = ["derive"] }
walkdir = "2"
serde = { version = "1", features = ["derive"] }
serde_json = "1"
rayon = "1"
thiserror = "2"
tokio = { version = "1", features = ["rt-multi-thread", "fs", "sync", "macros"] }
futures = "0.3"

[dev-dependencies]
tempfile = "3"
```

- [ ] **Step 2: Verify the workspace still builds**

Run: `cargo build -p pcat-pipeline`
Expected: builds clean (no new code yet, just new deps resolve).

- [ ] **Step 3: Commit**

```bash
git add crates/pcat-pipeline/Cargo.toml Cargo.lock
git commit -m "chore(pipeline): add tokio, futures, tempfile for async DICOM loader"
```

---

## Task 2: DicomLoadError enum

**Files:**
- Create: `crates/pcat-pipeline/src/dicom_errors.rs`
- Modify: `crates/pcat-pipeline/src/lib.rs`

- [ ] **Step 1: Write the failing test**

Create `crates/pcat-pipeline/src/dicom_errors.rs`:

```rust
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
```

Add to `crates/pcat-pipeline/src/lib.rs` — insert after `pub mod dicom_loader;`:

```rust
pub mod dicom_errors;
```

- [ ] **Step 2: Run tests to verify they pass (code is complete)**

Run: `cargo test -p pcat-pipeline dicom_errors::tests`
Expected: 3 passes.

- [ ] **Step 3: Commit**

```bash
git add crates/pcat-pipeline/src/dicom_errors.rs crates/pcat-pipeline/src/lib.rs
git commit -m "feat(pipeline): add DicomLoadError enum"
```

---

## Task 3: SliceHeader + read_header (header-only DICOM read)

**Files:**
- Create: `crates/pcat-pipeline/src/dicom_scan.rs`
- Modify: `crates/pcat-pipeline/src/lib.rs`

This introduces a new module `dicom_scan` that holds the header-only scan primitives, keeping the rewrite isolated from the legacy `dicom_loader.rs` until we swap over.

- [ ] **Step 1: Write the failing test**

Create `crates/pcat-pipeline/src/dicom_scan.rs`:

```rust
//! Header-only DICOM scanning primitives.
//!
//! `read_header` opens a DICOM file, stops parsing at the PixelData tag, and
//! returns a `SliceHeader` with only the tags we care about for indexing and
//! grouping. On an SMB share this transfers ~4 KB per file instead of ~512 KB.

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

    Ok(Some(SliceHeader {
        path: path.to_path_buf(),
        series_uid,
        series_description: read_string(&obj, tags::SERIES_DESCRIPTION).unwrap_or_default(),
        image_comments: read_string(&obj, IMAGE_COMMENTS),
        instance_number: read_i32(&obj, tags::INSTANCE_NUMBER),
        image_position_z: ipp.get(2).copied(),
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
}
```

Add to `crates/pcat-pipeline/src/lib.rs` after `pub mod dicom_errors;`:

```rust
pub mod dicom_scan;
```

- [ ] **Step 2: Run tests**

Run: `cargo test -p pcat-pipeline dicom_scan::tests`
Expected: 2 passes.

- [ ] **Step 3: Commit**

```bash
git add crates/pcat-pipeline/src/dicom_scan.rs crates/pcat-pipeline/src/lib.rs
git commit -m "feat(pipeline): add read_header for DICOM header-only scan"
```

---

## Task 4: Integration test fixture (synthetic mini-CT)

**Files:**
- Create: `crates/pcat-pipeline/tests/common/mod.rs`
- Create: `crates/pcat-pipeline/tests/fixture_smoke.rs`

This writes valid DICOM files using `dicom-rs` so the integration tests run against real parseable DICOMs without committing binaries.

- [ ] **Step 1: Write the fixture helper**

Create `crates/pcat-pipeline/tests/common/mod.rs`:

```rust
//! Test fixtures — generates valid synthetic DICOM files for integration tests.
//!
//! The fixture is a 30-slice mini-CT at 64x64 resolution with a simple radial
//! gradient pattern. Each slice is a valid Explicit-VR Little Endian DICOM file
//! with uncompressed 16-bit pixel data, rescale slope/intercept set, and
//! ImagePositionPatient encoding the Z axis.

use std::path::{Path, PathBuf};

use dicom::core::{DataElement, Length, PrimitiveValue, VR};
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
                .media_storage_sop_instance_uid(&sop_uid)
                .build()
                .unwrap();

            let file_obj = obj.with_meta(meta);
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
```

Create `crates/pcat-pipeline/tests/fixture_smoke.rs`:

```rust
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
```

- [ ] **Step 2: Run the smoke test**

Run: `cargo test -p pcat-pipeline --test fixture_smoke`
Expected: 2 passes.

- [ ] **Step 3: Commit**

```bash
git add crates/pcat-pipeline/tests/common/ crates/pcat-pipeline/tests/fixture_smoke.rs
git commit -m "test(pipeline): add synthetic mini-CT fixture for DICOM loader tests"
```

---

## Task 5: group_by_series with IPP-based sort

**Files:**
- Modify: `crates/pcat-pipeline/src/dicom_scan.rs`

- [ ] **Step 1: Write the failing test**

Append to `crates/pcat-pipeline/src/dicom_scan.rs` — inside the existing `mod tests { ... }`:

```rust
    fn h(uid: &str, z: Option<f64>, inst: Option<i32>, path: &str) -> SliceHeader {
        SliceHeader {
            path: std::path::PathBuf::from(path),
            series_uid: uid.to_string(),
            series_description: String::new(),
            image_comments: None,
            instance_number: inst,
            image_position_z: z,
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
```

Also at the top of `dicom_scan.rs` add the import:

```rust
use std::collections::HashMap;
```

- [ ] **Step 2: Run to confirm failure (function not defined)**

Run: `cargo test -p pcat-pipeline dicom_scan::tests`
Expected: FAIL — `group_by_series` not found.

- [ ] **Step 3: Implement group_by_series**

Append to `crates/pcat-pipeline/src/dicom_scan.rs`:

```rust
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
```

- [ ] **Step 4: Run to confirm passing**

Run: `cargo test -p pcat-pipeline dicom_scan::tests`
Expected: 5 passes (2 from Task 3 + 3 new).

- [ ] **Step 5: Commit**

```bash
git add crates/pcat-pipeline/src/dicom_scan.rs
git commit -m "feat(pipeline): group_by_series with IPP-then-InstanceNumber sort"
```

---

## Task 6: scan_series public API (tokio + bounded concurrency)

**Files:**
- Modify: `crates/pcat-pipeline/src/dicom_scan.rs`
- Create: `crates/pcat-pipeline/tests/dicom_loader_integration.rs`

- [ ] **Step 1: Write the failing integration test**

Create `crates/pcat-pipeline/tests/dicom_loader_integration.rs`:

```rust
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
```

- [ ] **Step 2: Run to confirm failure**

Run: `cargo test -p pcat-pipeline --test dicom_loader_integration`
Expected: FAIL — `scan_series` not found, also `SeriesDescriptor` missing several fields.

- [ ] **Step 3: Add SeriesDescriptor + scan_series**

Append to `crates/pcat-pipeline/src/dicom_scan.rs`:

```rust
use std::path::PathBuf;

use futures::stream::{self, StreamExt};
use tokio::sync::Semaphore;

/// Concurrent header opens. Empirically 32–64 is the SMB sweet spot; 48 is a
/// conservative middle value.
const SCAN_CONCURRENCY: usize = 48;

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
    }
}
```

- [ ] **Step 4: Run to confirm passing**

Run: `cargo test -p pcat-pipeline --test dicom_loader_integration`
Expected: 4 passes.

- [ ] **Step 5: Commit**

```bash
git add crates/pcat-pipeline/src/dicom_scan.rs crates/pcat-pipeline/tests/dicom_loader_integration.rs
git commit -m "feat(pipeline): scan_series with tokio + Semaphore(48) parallel header scan"
```

---

## Task 7: decode_slice_i16 function

**Files:**
- Create: `crates/pcat-pipeline/src/dicom_decode.rs`
- Modify: `crates/pcat-pipeline/src/lib.rs`

- [ ] **Step 1: Write the failing test**

Create `crates/pcat-pipeline/src/dicom_decode.rs`:

```rust
//! Pixel decode: open a DICOM file, read raw pixel data, apply rescale, clamp to i16.

use std::path::Path;

use dicom::object::open_file;
use dicom_pixeldata::PixelDecoder;

use crate::dicom_errors::DicomLoadError;

/// Decode one slice's pixel data and apply the rescale transform,
/// returning a flat row-major `Vec<i16>` of length `rows * cols`.
pub fn decode_slice_i16(
    path: &Path,
    rescale_slope: f64,
    rescale_intercept: f64,
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

    // Convert to ndarray of i16/u16 and apply rescale. We go through an f64
    // intermediate because slope/intercept can be non-integer (rare but legal).
    let ndarr = decoded
        .to_ndarray::<i32>()
        .map_err(|e| DicomLoadError::DecodeFailed {
            path: path.to_path_buf(),
            reason: e.to_string(),
        })?;
    let mut out = Vec::with_capacity((rows * cols) as usize);
    for &raw in ndarr.iter() {
        let hu = (raw as f64) * rescale_slope + rescale_intercept;
        let clamped = hu.clamp(i16::MIN as f64, i16::MAX as f64);
        out.push(clamped.round() as i16);
    }
    Ok(out)
}
```

Add to `crates/pcat-pipeline/src/lib.rs` after `pub mod dicom_scan;`:

```rust
pub mod dicom_decode;
```

Add a decode integration test to `crates/pcat-pipeline/tests/dicom_loader_integration.rs`:

```rust
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
```

- [ ] **Step 2: Run tests**

Run: `cargo test -p pcat-pipeline --test dicom_loader_integration`
Expected: 5 passes (previous 4 + new decode test).

- [ ] **Step 3: Commit**

```bash
git add crates/pcat-pipeline/src/dicom_decode.rs crates/pcat-pipeline/src/lib.rs crates/pcat-pipeline/tests/dicom_loader_integration.rs
git commit -m "feat(pipeline): decode_slice_i16 with rescale + dim check"
```

---

## Task 8: load_series public API (rayon parallel decode)

**Files:**
- Create: `crates/pcat-pipeline/src/dicom_load.rs`
- Modify: `crates/pcat-pipeline/src/lib.rs`
- Modify: `crates/pcat-pipeline/tests/dicom_loader_integration.rs`

- [ ] **Step 1: Write the failing integration test**

Append to `crates/pcat-pipeline/tests/dicom_loader_integration.rs`:

```rust
#[tokio::test]
async fn load_series_returns_correct_volume() {
    use pcat_pipeline::dicom_load::{load_series, LoadedVolume};
    use pcat_pipeline::dicom_scan::scan_series;

    let dir = tempfile::tempdir().unwrap();
    common::write_mini_ct(dir.path());

    let series = scan_series(dir.path()).await.unwrap();
    let vol: LoadedVolume = load_series(dir.path(), &series[0].uid).await.unwrap();

    assert_eq!(
        vol.metadata.num_slices,
        common::FIXTURE_SLICES,
    );
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

    let err = pcat_pipeline::dicom_load::load_series(dir.path(), "not-a-real-uid")
        .await
        .unwrap_err();
    assert!(format!("{err}").contains("not-a-real-uid"));
}

#[tokio::test]
async fn load_series_rejects_oversized_volume() {
    // Synthesize a descriptor with huge dims and pass it through the size guard.
    use pcat_pipeline::dicom_load::check_volume_size_mb;
    let mb = 1024 * 1024 * 1024usize; // 1 GB hypothetical
    assert!(check_volume_size_mb(mb).is_err());
}
```

- [ ] **Step 2: Run to confirm failure**

Run: `cargo test -p pcat-pipeline --test dicom_loader_integration`
Expected: FAIL — `dicom_load` module missing.

- [ ] **Step 3: Implement load_series**

Create `crates/pcat-pipeline/src/dicom_load.rs`:

```rust
//! Parallel pixel decode + volume assembly.
//!
//! Given a folder and a series UID, walks headers (to locate files), then
//! decodes pixel data in parallel (rayon) and returns a densely packed i16
//! volume in z-major order.

use std::path::Path;

use rayon::prelude::*;

use crate::dicom_decode::decode_slice_i16;
use crate::dicom_errors::DicomLoadError;
use crate::dicom_scan::{scan_series, SeriesDescriptor};

/// 4 GB soft limit (conservative — covers 1000-slice 512² i16 at 1.5 GB).
const VOLUME_SIZE_LIMIT_MB: usize = 4096;

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
        }
    }
}

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
pub async fn load_series(
    dir: &Path,
    uid: &str,
) -> Result<LoadedVolume, DicomLoadError> {
    let descriptors = scan_series(dir).await?;
    let desc = descriptors
        .into_iter()
        .find(|d| d.uid == uid)
        .ok_or_else(|| DicomLoadError::SeriesNotFound { uid: uid.to_string() })?;

    let slice_len = (desc.rows as usize) * (desc.cols as usize);
    let total_voxels = slice_len * desc.num_slices;
    let total_bytes_mb = (total_voxels * 2) / (1024 * 1024);
    check_volume_size_mb(total_bytes_mb)?;

    let file_paths = desc.file_paths.clone();
    let rescale_slope = desc.rescale_slope;
    let rescale_intercept = desc.rescale_intercept;
    let rows = desc.rows;
    let cols = desc.cols;

    // Rayon is CPU + I/O mixed here; its work-stealing pool saturates I/O just
    // fine for our file counts. Running inside spawn_blocking keeps the tokio
    // reactor responsive.
    let voxels = tokio::task::spawn_blocking(move || {
        let mut out = vec![0i16; total_voxels];
        // Decode slices in parallel, write each into its z-offset slot.
        let results: Vec<Result<(usize, Vec<i16>), DicomLoadError>> = file_paths
            .par_iter()
            .enumerate()
            .map(|(z, p)| {
                decode_slice_i16(p, rescale_slope, rescale_intercept, rows, cols)
                    .map(|px| (z, px))
            })
            .collect();
        for r in results {
            let (z, px) = r?;
            out[z * slice_len..(z + 1) * slice_len].copy_from_slice(&px);
        }
        Ok::<_, DicomLoadError>(out)
    })
    .await
    .map_err(|e| DicomLoadError::ParseFailed {
        path: dir.to_path_buf(),
        reason: format!("decode task panicked: {e}"),
    })??;

    Ok(LoadedVolume {
        metadata: VolumeMetadata::from(&desc),
        voxels_i16: voxels,
    })
}
```

Add to `crates/pcat-pipeline/src/lib.rs` after `pub mod dicom_decode;`:

```rust
pub mod dicom_load;
```

- [ ] **Step 4: Run tests**

Run: `cargo test -p pcat-pipeline --test dicom_loader_integration`
Expected: 8 passes.

- [ ] **Step 5: Commit**

```bash
git add crates/pcat-pipeline/src/dicom_load.rs crates/pcat-pipeline/src/lib.rs crates/pcat-pipeline/tests/dicom_loader_integration.rs
git commit -m "feat(pipeline): load_series with rayon parallel pixel decode"
```

---

## Task 9: Benchmark harness

**Files:**
- Create: `crates/pcat-pipeline/benches/dicom_load_bench.rs`
- Modify: `crates/pcat-pipeline/Cargo.toml`

- [ ] **Step 1: Add benchmark stanza**

Append to `crates/pcat-pipeline/Cargo.toml`:

```toml
[[bench]]
name = "dicom_load_bench"
harness = false
```

- [ ] **Step 2: Write the bench**

Create `crates/pcat-pipeline/benches/dicom_load_bench.rs`:

```rust
//! Wall-clock benchmark for scan_series + load_series.
//!
//! By default, runs against the synthetic fixture (30 slices × 64²).
//! To measure real SMB performance, set BENCH_DICOM_DIR to a folder of your choice:
//!
//!     BENCH_DICOM_DIR=/Volumes/labshare/patient_xyz cargo bench -p pcat-pipeline

use std::path::PathBuf;
use std::time::Instant;

use pcat_pipeline::dicom_load::load_series;
use pcat_pipeline::dicom_scan::scan_series;

#[path = "../tests/common/mod.rs"]
mod fixture;

fn main() {
    let runtime = tokio::runtime::Runtime::new().unwrap();
    runtime.block_on(async {
        let (dir_owner, dir) = resolve_dir();
        eprintln!("Benchmarking dir: {}", dir.display());

        let t0 = Instant::now();
        let series = scan_series(&dir).await.expect("scan failed");
        let t_scan = t0.elapsed();

        println!("scan_series: {:.3} s ({} series)", t_scan.as_secs_f64(), series.len());
        for (i, s) in series.iter().enumerate() {
            println!("  [{i}] {} ({} slices, {}×{})", s.description, s.num_slices, s.rows, s.cols);
        }

        let first = series.first().expect("at least one series");
        let t1 = Instant::now();
        let vol = load_series(&dir, &first.uid).await.expect("load failed");
        let t_load = t1.elapsed();
        let total = t_scan + t_load;
        let mb = (vol.voxels_i16.len() * 2) / (1024 * 1024);
        println!(
            "load_series: {:.3} s ({} slices, {} MB)",
            t_load.as_secs_f64(), vol.metadata.num_slices, mb,
        );
        println!("TOTAL cold load: {:.3} s", total.as_secs_f64());

        drop(dir_owner);
    });
}

// Returns (TempDir guard, path). TempDir is dropped at end of main to clean up.
fn resolve_dir() -> (Option<tempfile::TempDir>, PathBuf) {
    if let Ok(d) = std::env::var("BENCH_DICOM_DIR") {
        return (None, PathBuf::from(d));
    }
    let td = tempfile::tempdir().expect("tempdir");
    fixture::write_mini_ct(td.path());
    let p = td.path().to_path_buf();
    (Some(td), p)
}
```

- [ ] **Step 3: Run the bench against the fixture**

Run: `cargo bench -p pcat-pipeline --bench dicom_load_bench`
Expected: prints timings. Fixture should finish in well under 1 s total.

- [ ] **Step 4: Commit**

```bash
git add crates/pcat-pipeline/Cargo.toml crates/pcat-pipeline/benches/dicom_load_bench.rs
git commit -m "test(pipeline): benchmark harness for scan+load cold path"
```

---

## Task 10: Tauri framed-response helper

**Files:**
- Create: `src-tauri/src/commands/framed.rs`
- Modify: `src-tauri/src/commands/mod.rs`

- [ ] **Step 1: Write the failing test**

Create `src-tauri/src/commands/framed.rs`:

```rust
//! Framed binary response for Tauri commands.
//!
//! Layout (single buffer):
//!   [u32 LE: metadata_json_length] [metadata_json_bytes] [payload_bytes]
//!
//! Frontend parses the first 4 bytes, reads JSON, treats the rest as the
//! payload ArrayBuffer.

use serde::Serialize;

/// Encode (metadata, payload_bytes) into a single framed Vec<u8>.
pub fn encode_frame<M: Serialize>(metadata: &M, payload: &[u8]) -> Result<Vec<u8>, String> {
    let json = serde_json::to_vec(metadata).map_err(|e| format!("serialize metadata: {e}"))?;
    let meta_len: u32 = json
        .len()
        .try_into()
        .map_err(|_| "metadata json exceeds u32 range".to_string())?;

    let mut out = Vec::with_capacity(4 + json.len() + payload.len());
    out.extend_from_slice(&meta_len.to_le_bytes());
    out.extend_from_slice(&json);
    out.extend_from_slice(payload);
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[derive(Serialize)]
    struct M { a: u32, b: String }

    #[test]
    fn round_trips_payload() {
        let m = M { a: 42, b: "hi".into() };
        let payload = [0x01u8, 0x02, 0x03, 0x04];
        let buf = encode_frame(&m, &payload).unwrap();

        // Read back
        let meta_len = u32::from_le_bytes(buf[..4].try_into().unwrap()) as usize;
        let json = &buf[4..4 + meta_len];
        let body = &buf[4 + meta_len..];
        let parsed: serde_json::Value = serde_json::from_slice(json).unwrap();
        assert_eq!(parsed["a"], 42);
        assert_eq!(parsed["b"], "hi");
        assert_eq!(body, &payload);
    }

    #[test]
    fn encodes_empty_payload() {
        let m = M { a: 1, b: "".into() };
        let buf = encode_frame(&m, &[]).unwrap();
        assert!(buf.len() >= 5);
        let meta_len = u32::from_le_bytes(buf[..4].try_into().unwrap()) as usize;
        assert_eq!(buf.len(), 4 + meta_len);
    }
}
```

Edit `src-tauri/src/commands/mod.rs`. Replace contents with:

```rust
pub mod annotation;
pub mod cpr;
pub mod dicom;
pub mod framed;
pub mod pipeline;
```

(Note: `volume` is removed. If any other file re-exports the `volume` module, they will need a corresponding update in Task 13.)

- [ ] **Step 2: Run tests**

Run: `cargo test -p pcat-workstation-v2 framed`
Expected: 2 passes.

- [ ] **Step 3: Commit**

```bash
git add src-tauri/src/commands/framed.rs src-tauri/src/commands/mod.rs
git commit -m "feat(tauri): framed binary response helper for bulk payloads"
```

---

## Task 11: Tauri `scan_series` command (new shape)

**Files:**
- Modify: `src-tauri/src/commands/dicom.rs`
- Modify: `src-tauri/src/lib.rs` (or wherever commands are registered)

- [ ] **Step 1: Add the new scan_series command**

In `src-tauri/src/commands/dicom.rs`:

Add at the top of the file (after existing `use` lines):

```rust
use pcat_pipeline::dicom_scan::{self, SeriesDescriptor};
use std::path::PathBuf;
```

Append this function at the end of the file:

```rust
/// Scan a DICOM folder for series (new shape — replaces legacy scan).
/// Returns header-only metadata per series; pixel data is not decoded.
#[tauri::command]
pub async fn scan_series_v2(path: String) -> Result<Vec<SeriesDescriptorDto>, String> {
    let dir = PathBuf::from(path);
    let series = tokio::task::spawn_blocking(move || {
        tokio::runtime::Handle::current().block_on(dicom_scan::scan_series(&dir))
    })
    .await
    .map_err(|e| format!("scan task failed: {e}"))?
    .map_err(|e| e.to_string())?;

    Ok(series.into_iter().map(SeriesDescriptorDto::from).collect())
}

#[derive(serde::Serialize)]
pub struct SeriesDescriptorDto {
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
    /// Absolute file paths in z-sorted order.
    pub file_paths: Vec<String>,
    pub slice_positions_z: Vec<f64>,
}

impl From<SeriesDescriptor> for SeriesDescriptorDto {
    fn from(d: SeriesDescriptor) -> Self {
        Self {
            uid: d.uid,
            description: d.description,
            image_comments: d.image_comments,
            rows: d.rows,
            cols: d.cols,
            num_slices: d.num_slices,
            pixel_spacing: d.pixel_spacing,
            slice_spacing: d.slice_spacing,
            orientation: d.orientation,
            rescale_slope: d.rescale_slope,
            rescale_intercept: d.rescale_intercept,
            window_center: d.window_center,
            window_width: d.window_width,
            patient_name: d.patient_name,
            study_description: d.study_description,
            file_paths: d.file_paths.into_iter().map(|p| p.to_string_lossy().into_owned()).collect(),
            slice_positions_z: d.slice_positions_z,
        }
    }
}
```

- [ ] **Step 2: Register the command**

Find the `tauri::Builder::default().invoke_handler(tauri::generate_handler![...])` block (likely in `src-tauri/src/lib.rs` or `main.rs`). Add `commands::dicom::scan_series_v2` to the list.

- [ ] **Step 3: Build + ensure old and new coexist**

Run: `cargo check -p pcat-workstation-v2`
Expected: compiles. (Old `scan_series` is still present; we rename/delete in Task 13.)

- [ ] **Step 4: Commit**

```bash
git add src-tauri/src/commands/dicom.rs src-tauri/src/lib.rs
git commit -m "feat(tauri): scan_series_v2 command returning rich SeriesDescriptor"
```

---

## Task 12: Tauri `load_series` command (bulk binary response)

**Files:**
- Modify: `src-tauri/src/commands/dicom.rs`
- Modify: `src-tauri/src/lib.rs`

- [ ] **Step 1: Add the bulk-load command**

Append to `src-tauri/src/commands/dicom.rs`:

```rust
use pcat_pipeline::dicom_load;
use tauri::ipc::Response;
use crate::commands::framed::encode_frame;

/// Load a single series as one framed binary response:
///   [u32 LE: metadata_json_length] [metadata_json] [i16 LE voxel bytes]
/// Frontend receives this as an ArrayBuffer.
#[tauri::command]
pub async fn load_series_v2(dir: String, uid: String) -> Result<Response, String> {
    let dir_path = PathBuf::from(dir);
    let uid_clone = uid.clone();
    let vol = dicom_load::load_series(&dir_path, &uid_clone)
        .await
        .map_err(|e| e.to_string())?;

    let voxel_bytes: Vec<u8> = bytemuck::cast_slice(&vol.voxels_i16).to_vec();
    let framed = encode_frame(&vol.metadata, &voxel_bytes)?;
    Ok(Response::new(framed))
}
```

- [ ] **Step 2: Register the command**

Add `commands::dicom::load_series_v2` to the `invoke_handler!` macro list.

- [ ] **Step 3: Build**

Run: `cargo check -p pcat-workstation-v2`
Expected: compiles clean.

- [ ] **Step 4: Quick Rust-side integration smoke test**

Create `src-tauri/tests/bulk_load_command.rs`:

```rust
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
    let vol = load_series(dir.path(), &series[0].uid).await.unwrap();

    // Frame it like the Tauri command does.
    let payload: Vec<u8> = bytemuck::cast_slice(&vol.voxels_i16).to_vec();
    let frame = pcat_workstation_v2_lib::commands::framed::encode_frame(&vol.metadata, &payload)
        .expect("frame");

    // Parse it back the way the frontend will.
    let meta_len = u32::from_le_bytes(frame[..4].try_into().unwrap()) as usize;
    let meta_json = &frame[4..4 + meta_len];
    let body = &frame[4 + meta_len..];

    let meta: pcat_pipeline::dicom_load::VolumeMetadata =
        serde_json::from_slice(meta_json).unwrap();
    assert_eq!(meta.num_slices, fixture::FIXTURE_SLICES);
    assert_eq!(body.len(), vol.voxels_i16.len() * 2);
}
```

Ensure `pcat-workstation-v2` re-exports `commands` from `lib.rs`. If not already, add:

```rust
pub mod commands;
```

to `src-tauri/src/lib.rs` (or confirm an existing export).

Add `tempfile = "3"` to `src-tauri/Cargo.toml` under `[dev-dependencies]`.

Run: `cargo test -p pcat-workstation-v2 --test bulk_load_command`
Expected: 1 pass.

- [ ] **Step 5: Commit**

```bash
git add src-tauri/src/commands/dicom.rs src-tauri/src/lib.rs src-tauri/tests/bulk_load_command.rs src-tauri/Cargo.toml
git commit -m "feat(tauri): load_series_v2 returning framed ArrayBuffer response"
```

---

## Task 13: Delete legacy DICOM commands and volume.rs

**Files:**
- Delete: `src-tauri/src/commands/volume.rs`
- Modify: `src-tauri/src/commands/mod.rs`
- Modify: `src-tauri/src/commands/dicom.rs`
- Modify: `src-tauri/src/lib.rs`
- Modify: `src-tauri/src/state.rs`
- Modify: `crates/pcat-pipeline/src/dicom_loader.rs` (strip legacy helpers that are now unused)

- [ ] **Step 1: Find all call sites of legacy commands**

Run: `rg "load_dicom|get_slice|load_dual_energy|scan_series\b" src-tauri/src/ src/ 2>&1 | head -60`
Record the results. This shows every place the old commands are wired.

- [ ] **Step 2: Delete volume.rs**

```bash
git rm src-tauri/src/commands/volume.rs
```

- [ ] **Step 3: Remove legacy commands from dicom.rs**

Edit `src-tauri/src/commands/dicom.rs`:

- Delete `load_dicom` (the entire `#[tauri::command]` block currently around line 68).
- Delete `load_dual_energy` and any associated helpers.
- Delete the old `scan_series` (the one that returns the short `SeriesInfo`).
- **Rename** `scan_series_v2` → `scan_series` and `load_series_v2` → `load_series` (now that there's no name conflict).
- **Delete** `VolumeInfo` struct (no longer returned).

Corresponding renames inside the registered command macro list.

- [ ] **Step 4: Update state.rs**

Edit `src-tauri/src/state.rs`:

- Remove the `volume: Option<LoadedVolume>` field and any fields depending on the old `LoadedVolume` type (the old `ndarray::Array3<f32>` container).
- If `AppState` becomes empty, keep it as `pub struct AppState {}` for now — other code likely still locks it.
- Remove `use` of `ndarray::Array3` if it was only for the volume.

If some downstream code (e.g., `pipeline.rs`) still reads `state.volume`, it must be updated to accept a volume passed in from the frontend (via a new command signature). For now, add `#[allow(dead_code)]` stubs or remove those consumers as needed. Note any flagged callers in a comment for Task 17.

- [ ] **Step 5: Strip legacy helpers from dicom_loader.rs**

Edit `crates/pcat-pipeline/src/dicom_loader.rs`:

- Delete `load_dicom_directory`, `DualEnergyVolume`, and any helper private functions (`read_multi_f64`, `read_f64`, `read_string`, `read_u16` — these are replicated in the new `dicom_scan.rs`; grep to confirm no other file imports them).
- **Keep** `parse_kev`, `parse_kev_from_image_comments`, `parse_kev_from_description` if used elsewhere — grep first.
- If any caller outside this file uses `SeriesInfo` (the legacy type), update them to use the new `SeriesDescriptor`.

Run: `cargo check -p pcat-pipeline -p pcat-workstation-v2`
Fix any remaining compile errors by rewiring callers to the new APIs.

- [ ] **Step 6: Commit**

```bash
git add -A src-tauri/src crates/pcat-pipeline/src
git commit -m "refactor: remove legacy DICOM commands and volume.rs after v2 cutover"
```

---

## Task 14: Update TypeScript types + api.ts

**Files:**
- Modify: `src/lib/api.ts`

- [ ] **Step 1: Redefine types**

Edit `src/lib/api.ts`:

Replace the old `VolumeInfo`, `SeriesInfo`, `DualEnergyInfo`, and the functions `loadDicom`, `getSlice`, `loadDualEnergy`, existing `scanSeries` with new definitions.

Paste this block in place of those declarations (near the top, after `import { invoke }`):

```ts
import { invoke } from '@tauri-apps/api/core';

export interface SeriesDescriptor {
  uid: string;
  description: string;
  image_comments: string | null;
  rows: number;
  cols: number;
  num_slices: number;
  pixel_spacing: [number, number];
  slice_spacing: number;
  orientation: [number, number, number, number, number, number];
  rescale_slope: number;
  rescale_intercept: number;
  window_center: number;
  window_width: number;
  patient_name: string;
  study_description: string;
  file_paths: string[];
  slice_positions_z: number[];
}

export interface VolumeMetadata {
  series_uid: string;
  series_description: string;
  image_comments: string | null;
  rows: number;
  cols: number;
  num_slices: number;
  pixel_spacing: [number, number];
  slice_spacing: number;
  orientation: [number, number, number, number, number, number];
  window_center: number;
  window_width: number;
  patient_name: string;
  study_description: string;
  slice_positions_z: number[];
}

/** Scan a DICOM folder (header-only) for available series. */
export async function scanSeries(dir: string): Promise<SeriesDescriptor[]> {
  return invoke<SeriesDescriptor[]>('scan_series', { path: dir });
}

/**
 * Load one series as a single bulk ArrayBuffer containing framed metadata + i16 voxels.
 * Returns parsed metadata and a Uint8Array over the voxel bytes (zero-copy).
 */
export async function loadSeries(
  dir: string,
  uid: string,
): Promise<{ metadata: VolumeMetadata; voxels: Int16Array }> {
  const buf = await invoke<ArrayBuffer>('load_series', { dir, uid });

  const view = new DataView(buf);
  const metaLen = view.getUint32(0, true);
  const metaBytes = new Uint8Array(buf, 4, metaLen);
  const metadata = JSON.parse(new TextDecoder().decode(metaBytes)) as VolumeMetadata;

  // Voxel data starts after [4 + metaLen] bytes. Int16Array requires alignment
  // on a 2-byte boundary — `buf` is a fresh ArrayBuffer from Tauri, so the
  // byteOffset (4 + metaLen) is typically aligned; if metaLen is odd we copy.
  const voxelOffset = 4 + metaLen;
  let voxels: Int16Array;
  if (voxelOffset % 2 === 0) {
    voxels = new Int16Array(buf, voxelOffset);
  } else {
    const copy = new Uint8Array(buf.byteLength - voxelOffset);
    copy.set(new Uint8Array(buf, voxelOffset));
    voxels = new Int16Array(copy.buffer);
  }
  return { metadata, voxels };
}
```

Delete these functions (no longer exist backend-side):
- `loadDicom`
- `getSlice`
- `loadDualEnergy`
- the previous `scanSeries` and `SeriesInfo` types

Keep everything else (`openDicomDialog`, `getRecentDicoms`, `saveSeeds`, `loadSeeds`, annotation commands, etc.) unchanged.

- [ ] **Step 2: Type-check**

Run: `npm run check` (or the project's equivalent — check `package.json` `scripts`; if absent, run `npx svelte-check`).
Expected: errors come from existing consumers (App.svelte still calls `loadDicom` etc.) — note the failing file/line so Task 16 can address them all in one go.

- [ ] **Step 3: Commit**

```bash
git add src/lib/api.ts
git commit -m "feat(api): scanSeries + loadSeries with framed ArrayBuffer parse"
```

---

## Task 15: Rewrite cornerstone volumeLoader.ts

**Files:**
- Modify: `src/lib/cornerstone/volumeLoader.ts`
- Check: `src/lib/stores/volumeStore.svelte.ts` (for `VolumeMetadata` type consumers)

- [ ] **Step 1: Replace volumeLoader.ts**

Overwrite `src/lib/cornerstone/volumeLoader.ts` with:

```ts
/**
 * Volume loader — takes a pre-loaded bulk ArrayBuffer + metadata from the
 * backend and constructs a cornerstone3D local volume in one call.
 */
import { volumeLoader, cache, type Types } from '@cornerstonejs/core';

import type { VolumeMetadata } from '$lib/api';

type Point3 = Types.Point3;
type Mat3 = Types.Mat3;

/**
 * Construct a cornerstone3D local volume from a loaded series.
 * Returns the cornerstone volume ID.
 */
export function buildVolume(
  volumeKey: string,
  metadata: VolumeMetadata,
  voxels: Int16Array,
): string {
  const csVolumeId = `pcat:${volumeKey}`;

  const existing = cache.getVolume(csVolumeId);
  if (existing) {
    return csVolumeId;
  }

  // Rust metadata uses per-DICOM conventions (rows, cols, num_slices).
  // cornerstone3D uses [X, Y, Z] = [cols, rows, slices].
  const dimensions: Point3 = [metadata.cols, metadata.rows, metadata.num_slices];
  const spacing: Point3 = [
    metadata.pixel_spacing[1],   // sx (column spacing)
    metadata.pixel_spacing[0],   // sy (row spacing)
    metadata.slice_spacing,
  ];
  const origin: Point3 = [0, 0, metadata.slice_positions_z[0] ?? 0];

  // 3x3 direction matrix: rows = (iop_x, iop_y, normal).
  const iopRow: [number, number, number] = [
    metadata.orientation[0], metadata.orientation[1], metadata.orientation[2],
  ];
  const iopCol: [number, number, number] = [
    metadata.orientation[3], metadata.orientation[4], metadata.orientation[5],
  ];
  const normal: [number, number, number] = [
    iopRow[1] * iopCol[2] - iopRow[2] * iopCol[1],
    iopRow[2] * iopCol[0] - iopRow[0] * iopCol[2],
    iopRow[0] * iopCol[1] - iopRow[1] * iopCol[0],
  ];
  const direction = [
    iopRow[0], iopRow[1], iopRow[2],
    iopCol[0], iopCol[1], iopCol[2],
    normal[0], normal[1], normal[2],
  ] as Mat3;

  volumeLoader.createLocalVolume(csVolumeId, {
    scalarData: voxels,
    metadata: {
      BitsAllocated: 16,
      BitsStored: 16,
      SamplesPerPixel: 1,
      HighBit: 15,
      PhotometricInterpretation: 'MONOCHROME2',
      PixelRepresentation: 1,
      Modality: 'CT',
      ImageOrientationPatient: Array.from(direction.slice(0, 6)),
      PixelSpacing: [metadata.pixel_spacing[0], metadata.pixel_spacing[1]],
      FrameOfReferenceUID: `1.2.826.0.1.3680043.8.498.pcat`,
      Columns: metadata.cols,
      Rows: metadata.rows,
      voiLut: [{ windowCenter: metadata.window_center, windowWidth: metadata.window_width }],
      VOILUTFunction: 'LINEAR',
    },
    dimensions,
    spacing,
    origin,
    direction,
  });

  return csVolumeId;
}
```

- [ ] **Step 2: Update consumers of the old signature**

Run: `rg "loadVolume\(|getSlice\(|FETCH_CONCURRENCY|runWithConcurrency" src/`
Record the list. These are all App.svelte / store call sites that Task 16 will update.

Check `src/lib/stores/volumeStore.svelte.ts`:
- If it exports a `VolumeMetadata` type that shadows the new one in `api.ts`, update the store to re-export from `api.ts`:

```ts
export type { VolumeMetadata } from '$lib/api';
```

- [ ] **Step 3: Type-check**

Run: `npx svelte-check --tsconfig tsconfig.json`
Expected: errors only in App.svelte (old flow); no errors in volumeLoader.ts itself.

- [ ] **Step 4: Commit**

```bash
git add src/lib/cornerstone/volumeLoader.ts src/lib/stores/volumeStore.svelte.ts
git commit -m "refactor(cornerstone): buildVolume takes pre-loaded Int16Array"
```

---

## Task 16: Update App.svelte single-series DICOM load flow

**Files:**
- Modify: `src/App.svelte`

- [ ] **Step 1: Find the current load flow**

Run: `rg -n "openDicomDialog|loadDicom\(|loadDualEnergy\(|getSlice\(|scanSeries\(|loadVolume\(" src/App.svelte`
Record line numbers of everywhere the old flow is called.

- [ ] **Step 2: Replace the single-series flow**

Locate in `src/App.svelte` the block (roughly) that looks like:

```ts
const info = await loadDicom(dir);
const volumeKey = ...;
const csId = await loadVolume(info, onProgress);
```

Replace with:

```ts
// Scan first — shows which series are in the folder.
const series = await scanSeries(dir);
if (series.length === 0) {
  throw new Error('No DICOM series found in folder.');
}

// Default behavior: if only one series, auto-select it. Otherwise surface
// a picker (wired in the series-picker UI; if the old UI unconditionally
// loaded the first series, preserve that and leave the picker work for
// a follow-up).
const chosen = series[0];
const { metadata, voxels } = await loadSeries(dir, chosen.uid);
const volumeKey = `${chosen.patient_name}-${chosen.uid}`.replace(/[^A-Za-z0-9-_]/g, '_');
const csId = buildVolume(volumeKey, metadata, voxels);
```

Update the imports at the top of the `<script>` section:

```ts
import { scanSeries, loadSeries, openDicomDialog /* + existing imports */ } from '$lib/api';
import { buildVolume } from '$lib/cornerstone/volumeLoader';
```

Remove the `loadDicom`, `loadVolume`, and any `getSlice` imports.

If App.svelte holds volume metadata in a store (e.g. `volumeStore.setMetadata(info)`), adapt to set it from `metadata` — the field names differ:
- `info.shape` → `[metadata.num_slices, metadata.rows, metadata.cols]`
- `info.spacing` → `[metadata.slice_spacing, metadata.pixel_spacing[0], metadata.pixel_spacing[1]]`
- `info.window_center`/`window_width` → `metadata.window_center`/`window_width`
- `info.patient_name`/`study_description` → same

Update the store and any derived UI strings accordingly.

- [ ] **Step 3: Type-check**

Run: `npx svelte-check --tsconfig tsconfig.json`
Expected: passes (or only dual-energy-related errors remaining for Task 17).

- [ ] **Step 4: Commit**

```bash
git add src/App.svelte src/lib/stores/volumeStore.svelte.ts
git commit -m "refactor(ui): single-series load flow uses scanSeries + loadSeries"
```

---

## Task 17: Update App.svelte dual-energy / MMD flow

**Files:**
- Modify: `src/App.svelte`
- Modify: any MMD-specific components that called `loadDualEnergy`

- [ ] **Step 1: Find dual-energy call sites**

Run: `rg -n "loadDualEnergy\(|DualEnergyInfo\b" src/`
Record results.

- [ ] **Step 2: Replace dual-energy load with two parallel loadSeries calls**

In each call site, replace:

```ts
const de = await loadDualEnergy(path, lowUid, highUid, lowKev, highKev);
```

with:

```ts
const [low, high] = await Promise.all([
  loadSeries(path, lowUid),
  loadSeries(path, highUid),
]);
// low.metadata, low.voxels  →  the 70 keV volume
// high.metadata, high.voxels →  the 140 keV volume
```

Downstream code that used `de.shape`, `de.spacing`, `de.low_kev`, `de.high_kev`, etc. will now read:

```ts
const shape: [number, number, number] = [
  low.metadata.num_slices,
  low.metadata.rows,
  low.metadata.cols,
];
const spacing: [number, number, number] = [
  low.metadata.slice_spacing,
  low.metadata.pixel_spacing[0],
  low.metadata.pixel_spacing[1],
];
const lowKev = parseFloat(low.metadata.image_comments?.match(/([0-9.]+)\s*keV/i)?.[1] ?? '0');
const highKev = parseFloat(high.metadata.image_comments?.match(/([0-9.]+)\s*keV/i)?.[1] ?? '0');
```

If the MMD Rust command previously relied on `state.volume` / `state.dual_energy` being populated by `load_dual_energy`, update it to accept the voxel data (or file paths + uids) directly via its own Tauri command signature. If that requires a separate rework pass, add a TODO comment referencing this plan and hoist MMD integration into a follow-up task in the notes at the end of this plan — do not block DICOM loading on MMD.

- [ ] **Step 3: Type-check + build**

Run: `npx svelte-check --tsconfig tsconfig.json`
Expected: passes.

Run: `cargo build -p pcat-workstation-v2 --release`
Expected: passes.

- [ ] **Step 4: Commit**

```bash
git add -A src
git commit -m "refactor(mmd): dual-energy load = two parallel loadSeries calls"
```

---

## Task 18: End-to-end Playwright test

**Files:**
- Create: `tests/e2e/dicom-load.spec.ts` (if Playwright is already configured — check `playwright.config.ts` at repo root)

- [ ] **Step 1: Check if Playwright is configured**

Run: `ls playwright.config.* 2>/dev/null; ls tests/e2e/ 2>/dev/null`

- If neither exists, this task becomes "skip Playwright; add a short `vitest` integration test that mocks `@tauri-apps/api/core` with a fake `invoke`."
- If Playwright is configured, proceed to Step 2.

- [ ] **Step 2: Write the failing test**

Create `tests/e2e/dicom-load.spec.ts`:

```ts
import { test, expect } from '@playwright/test';

test('loading a mocked DICOM series renders the MPR viewport', async ({ page }) => {
  // Stub the Tauri invoke bridge so the test runs without a real backend.
  await page.addInitScript(() => {
    const fakeSeries = [{
      uid: '1.2.3.4',
      description: 'Fixture CT',
      image_comments: null,
      rows: 64, cols: 64, num_slices: 10,
      pixel_spacing: [0.5, 0.5],
      slice_spacing: 1.0,
      orientation: [1, 0, 0, 0, 1, 0],
      rescale_slope: 1.0, rescale_intercept: -1024.0,
      window_center: 40, window_width: 400,
      patient_name: 'TEST', study_description: 'FIX',
      file_paths: Array.from({ length: 10 }, (_, i) => `/tmp/s${i}.dcm`),
      slice_positions_z: Array.from({ length: 10 }, (_, i) => i),
    }];

    // Build a framed buffer: [u32 meta_len][meta_json][i16 voxel bytes]
    const meta = {
      series_uid: '1.2.3.4',
      series_description: 'Fixture CT',
      image_comments: null,
      rows: 64, cols: 64, num_slices: 10,
      pixel_spacing: [0.5, 0.5],
      slice_spacing: 1.0,
      orientation: [1, 0, 0, 0, 1, 0],
      window_center: 40, window_width: 400,
      patient_name: 'TEST', study_description: 'FIX',
      slice_positions_z: Array.from({ length: 10 }, (_, i) => i),
    };
    const metaJson = new TextEncoder().encode(JSON.stringify(meta));
    const sliceLen = 64 * 64;
    const voxels = new Int16Array(sliceLen * 10);
    for (let i = 0; i < voxels.length; i++) voxels[i] = (i % 1024) - 512;
    const voxelBytes = new Uint8Array(voxels.buffer);

    const buf = new ArrayBuffer(4 + metaJson.length + voxelBytes.length);
    const view = new DataView(buf);
    view.setUint32(0, metaJson.length, true);
    new Uint8Array(buf, 4, metaJson.length).set(metaJson);
    new Uint8Array(buf, 4 + metaJson.length).set(voxelBytes);

    // @ts-expect-error window.__TAURI_INTERNALS__ is a v2 internal hook
    window.__TAURI_INTERNALS__ = {
      invoke: async (cmd: string, args: unknown) => {
        if (cmd === 'open_dicom_dialog') return '/fake/dir';
        if (cmd === 'scan_series') return fakeSeries;
        if (cmd === 'load_series') return buf;
        if (cmd === 'get_recent_dicoms') return [];
        return null;
      },
    };
  });

  await page.goto('/');

  // Trigger the loader via the app's UI button.
  await page.getByRole('button', { name: /open dicom/i }).click();

  // The MPR viewport canvas should be non-blank within 3 s.
  const canvas = page.locator('canvas').first();
  await expect(canvas).toBeVisible({ timeout: 3000 });
});
```

- [ ] **Step 3: Run it**

Run: `npx playwright test tests/e2e/dicom-load.spec.ts`
Expected: PASS (or informative failure that identifies a UI-wiring regression).

- [ ] **Step 4: Commit**

```bash
git add tests/e2e/dicom-load.spec.ts
git commit -m "test(e2e): Playwright test for mocked DICOM load flow"
```

---

## Task 19: Manual benchmark on real SMB data

**Files:**
- Update: `docs/specs/2026-04-20-fast-dicom-loading-design.md` (append a "Measurements" appendix)

- [ ] **Step 1: Mount the SMB share**

Verify the SMB share is mounted at `/Volumes/labshare/Shu Nie/...` (per memory).

- [ ] **Step 2: Run the benchmark against a real patient**

Run: `BENCH_DICOM_DIR="/Volumes/labshare/Shu Nie/<pick-a-patient>" cargo bench -p pcat-pipeline --bench dicom_load_bench 2>&1 | tee /tmp/bench-new.txt`

Record `scan_series`, `load_series`, and TOTAL values.

- [ ] **Step 3: Run the dev app against the same patient**

Run: `npm run tauri dev` (or the repo's dev command).
Open the patient folder via the UI, time with a stopwatch the interval between "click open" and "first MPR frame visible."

Record the number.

- [ ] **Step 4: Append measurements to the spec**

Append to `docs/specs/2026-04-20-fast-dicom-loading-design.md`:

```markdown
## Measurements (2026-XX-XX)

Patient: `<folder name>` on `/Volumes/labshare/Shu Nie/...`, <N> slices, <rows>×<cols>.

| Phase | Time |
|---|---|
| scan_series | X.X s |
| load_series | X.X s |
| Rust cold total | X.X s |
| UI "click open" → first MPR frame | X.X s |

Goal was ≤ 10 s; stretch ≤ 5 s. **Actual: X.X s.**
```

- [ ] **Step 5: Commit the measurement record**

```bash
git add docs/specs/2026-04-20-fast-dicom-loading-design.md
git commit -m "docs(spec): record real-SMB benchmark after rewrite"
```

---

## Task 20: Cleanup pass — verify nothing else references dead code

**Files:**
- Various

- [ ] **Step 1: Sweep for dead references**

Run:
```bash
rg -n "getSlice|load_dicom\b|load_dual_energy|DualEnergyInfo|DualEnergyVolume|VolumeInfo\b|SeriesInfo\b|scan_dicom_series" src-tauri/src src crates/pcat-pipeline/src
```
Expected: no hits (other than in tests, plan docs, or legacy comments).

If there are hits, remove or rewire.

- [ ] **Step 2: Full workspace build + test**

Run: `cargo build --workspace --release` — expect success.
Run: `cargo test --workspace` — expect all tests pass.
Run: `npx svelte-check --tsconfig tsconfig.json` — expect clean.

- [ ] **Step 3: Manual verification checklist**

Go through the 5 steps in the spec's Manual Verification section:

1. Load a real patient from the SMB share. Measure wall clock. Confirm < 10 s for a ~300-slice study.
2. Load a dual-energy folder. Both series load in parallel via Promise.all. MMD panel works.
3. Pull the SMB cable mid-load (`diskutil unmount force /Volumes/labshare` or similar). Graceful error surfaces; UI doesn't crash.
4. Open a folder containing `.DS_Store` and a `README.txt`. They are silently skipped.
5. Open a MonoPlus folder. The keV picker reads `image_comments`, not `series_description`.

- [ ] **Step 4: Commit nothing (or a CHANGELOG entry)**

If a CHANGELOG exists at repo root:

```bash
echo "- DICOM loading is ~10-20x faster on SMB (stateless three-stage pipeline, bulk ArrayBuffer transport)" >> CHANGELOG.md
git add CHANGELOG.md
git commit -m "docs(changelog): record DICOM loader speedup"
```

---

## Notes

### Out-of-scope follow-ups (NOT part of this plan)

- **Multi-series picker UI.** If App.svelte currently auto-loads the first series, this plan preserves that. A proper "list series → click to load" UI is a follow-up.
- **MMD command rewire.** If MMD backend commands relied on backend state (`state.dual_energy`), migrating them to accept loaded voxel data is out of scope. The dual-energy volumes are now held in the frontend; MMD commands need to either accept them as args or be rewired to re-load from `load_series`. Keep them functional via a temporary bridge if needed; fully redesign in a follow-up.
- **Progressive / streaming loader.** Spec explicitly defers this. If the 5–10 s cold load proves to be a workflow blocker, add an Approach-1-style streaming protocol later.
- **Multi-frame / cine DICOMs.** Single-frame series only in this rewrite.

### Rollback strategy

If `tauri::ipc::Response` binary return path misbehaves at large payloads (≥ 500 MB), fall back to chunking: return `Vec<Vec<u8>>` of N ~50 MB chunks, still far fewer than 384. Keep `encode_frame` unchanged; adapt `load_series_v2` to split and adapt `api.ts::loadSeries` to concatenate. Chunking adds ~20 LOC.

### Success criteria recap

- Cold SMB load of a 300-slice 512² CT: **≤ 10 s** (stretch ≤ 5 s).
- Total LOC delta: **negative** (fewer lines after rewrite).
- IPC roundtrips per load: **2** (from ~390).
- No new persistent state on disk or in memory beyond the active volume.
- All tests green: unit, integration, e2e (Playwright), full `cargo test --workspace`.
