# Fast DICOM Loading — Design

**Date:** 2026-04-20
**Status:** Approved, awaiting implementation plan
**Scope:** PCAT Workstation v2 (Tauri + Svelte + Rust)

## Problem

Loading a DICOM series from the SMB-mounted lab share currently takes **40–85 s** for a typical 300-slice 512×512 cardiac CT. User reports it as "so long." The slowness has three structural root causes, not tunable parameters:

1. **Double directory walk** — `scan_dicom_series` walks the folder and opens every file to extract metadata; then `load_dicom_directory` walks the *same* folder and opens every file *again* to read pixel data.
2. **Serial file I/O** — the `rayon` crate is imported but never used; all `open_file` calls are in a single-threaded `WalkDir` iterator. On a ~50 ms-latency SMB link, 300 files × 50 ms = 15 s just for opens, before any parsing.
3. **384-roundtrip IPC** — after pixels decode, the frontend invokes `get_slice` once per slice, each call paying JSON serialization (`Vec<u8>` → `number[]` → `Uint8Array` → `Int16Array`). For 384 slices this is the single largest post-decode cost.

Additional friction: no header-only read path (pydicom-v1-equivalent); full-file fetch for every indexing operation; pixel data transferred even when only `SeriesInstanceUID` is needed.

## Goals and non-goals

**Goals**
- Cold first-load of a 300-slice CT from SMB in **≤ 10 s** (vs current 40–85 s). Target is a **5–10× speedup**.
- Robust to downstream change: DICOM transfer-syntax changes, folder layout changes, cornerstone3D API churn, Tauri version upgrades.
- Deterministic: same folder → same load time every time; no persistent cache, no staleness bugs.
- Simple code: fewer moving parts than the current implementation.

**Non-goals**
- No persistent on-disk cache. (User prefers deterministic cold path over cache complexity.)
- No progressive / first-pixel-ASAP UX in this revision. The user's workflow is "open one patient and stay" — the spinner is paid once per session.
- No DICOMweb / WADO-RS support. Local filesystem and SMB only.
- No progressive JPEG 2000 decoding. Legacy J2K handled by `dicom-rs`'s full-frame decoder; HTJ2K deferred.

## Architecture

Three stages, each independently replaceable:

```
┌──────────────────────────────────────────────────────────────┐
│  Svelte UI                                                    │
│  PatientBrowser → SeriesPicker → cornerstone3D MPR            │
└────────────┬────────────────────┬────────────────────────────┘
             │ scan_series        │ load_series
             ▼                    ▼
┌──────────────────────────────────────────────────────────────┐
│  Rust (pcat-pipeline crate)                                   │
│                                                               │
│  STAGE 1 — SCAN (header-only, parallel async)                 │
│    tokio::fs::read_dir + buffer_unordered(48)                 │
│    dicom-rs OpenFileOptions::read_until(tags::PIXEL_DATA)     │
│    group by SeriesInstanceUID, sort by IPP[2]                 │
│    → Vec<SeriesDescriptor>                                    │
│                                                               │
│  STAGE 2 — LOAD (pixel decode, parallel rayon)                │
│    rayon::par_iter over the selected series' file paths       │
│    dicom-rs decode_pixel_data + rescale slope/intercept       │
│    stack into Vec<i16> (z-major, slices sorted by IPP[2])     │
│                                                               │
│  STAGE 3 — TRANSFER (single framed IPC response)              │
│    tauri::ipc::Response with [u32 metadata_len][metadata_json]│
│    [i16 voxels]  — no per-voxel JSON encoding                 │
└──────────────────────────────────────────────────────────────┘
```

### Architectural properties

- **Narrow contracts between stages.** Stage 1 outputs `Vec<SeriesDescriptor>`. Stage 2 takes `(dir, series_uid)`, outputs `LoadedVolume`. Stage 3 is a Tauri command signature. Each stage can be rewritten without touching the others.
- **Stateless commands.** Backend holds no session state. User cancels mid-load by dropping the future; `tokio` and `rayon` clean up naturally. There is no `cancel_load` API because there is no state to cancel.
- **Concurrency bounded at one place each.** `Semaphore(48)` for async scan; `rayon` default thread pool for CPU-bound decode. No risk of fan-out storming the SMB server or the CPU.
- **No cache.** Every load is a cold load. Source of truth is the filesystem.

## Components

### `crates/pcat-pipeline/src/dicom_loader.rs` — rewrite

Public API (the only surface downstream consumers depend on):

```rust
pub struct SeriesDescriptor {
    pub uid: String,
    pub description: String,
    pub image_comments: Option<String>,  // MonoPlus keV truth
    pub rows: u32,
    pub cols: u32,
    pub num_slices: usize,
    pub pixel_spacing: [f64; 2],
    pub slice_spacing: f64,
    pub orientation: [f64; 6],
    pub rescale_slope: f64,
    pub rescale_intercept: f64,
    pub file_paths: Vec<PathBuf>,        // sorted by IPP[2]
}

pub struct LoadedVolume {
    pub metadata: VolumeMetadata,        // subset of SeriesDescriptor
    pub voxels_i16: Vec<i16>,            // tightly packed, z-major
}

pub async fn scan_series(dir: &Path) -> Result<Vec<SeriesDescriptor>, DicomLoadError>;
pub async fn load_series(dir: &Path, uid: &str) -> Result<LoadedVolume, DicomLoadError>;
```

Internal helpers:

- `read_header(path)` — `OpenFileOptions::new().read_until(tags::PIXEL_DATA).open_file(...)` returning a `SliceHeader` struct with only the tags we need.
- `decode_slice_i16(path, rescale)` — full open + `decode_pixel_data` + HU rescale + clamp to i16.
- `group_by_series(headers)` — `HashMap<String, Vec<SliceHeader>>`, each bucket sorted by `IPP[2]` (falling back to `InstanceNumber` when IPP is absent).

**Concurrency mix:** `tokio::fs` for scan (I/O-bound, latency-sensitive), `rayon` for load (CPU-bound decode mixed with I/O that benefits from parallelism). Each stage picks the right tool.

Estimated size: ~350 LOC (current: ~500 LOC).

### `src-tauri/src/commands/dicom.rs` — slim down

Two Tauri commands only:

```rust
#[tauri::command]
async fn scan_series(dir: String) -> Result<Vec<SeriesDto>, String>;

#[tauri::command]
async fn load_series(dir: String, uid: String) -> Result<tauri::ipc::Response, String>;
```

**Binary transport detail:** `load_series` returns a framed buffer via `tauri::ipc::Response`, which in Tauri v2 is passed to the frontend as a raw `ArrayBuffer` — no JSON encoding of the ~200 MB payload.

Frame layout (single buffer, no second round-trip for metadata):

```
offset 0:              u32 little-endian   metadata_json_length
offset 4:              UTF-8 JSON bytes    metadata_json
offset 4+meta_len:     i16 little-endian   voxel buffer (rows × cols × num_slices)
```

Dual-energy: no special command. The frontend invokes `load_series` twice in parallel. Rayon saturates across both loads.

Deleted commands:
- `get_slice` — redundant after bulk transfer.
- `load_dual_energy` — frontend orchestrates two `load_series` calls.
- Ad-hoc scan helpers inside `list_patients` paths that duplicate `scan_series` logic.

Estimated size: ~120 LOC (current: ~300 LOC).

### `src/lib/cornerstone/volumeLoader.ts` — rewrite

```ts
export async function loadVolume(dir: string, seriesUid: string): Promise<Volume> {
    const { metadata, bytes } = await loadSeries(dir, seriesUid);
    const scalarData = new Int16Array(bytes);  // zero-copy view
    return volumeLoader.createLocalVolume(volumeId, {
        dimensions: [metadata.cols, metadata.rows, metadata.num_slices],
        spacing: [...metadata.pixel_spacing, metadata.slice_spacing],
        orientation: metadata.orientation,
        scalarData,
    });
}
```

Uses cornerstone3D's boring, stable `createLocalVolume` API. No streaming machinery.

Estimated size: ~60 LOC (current: ~100 LOC with 384-IPC loop).

### `src/lib/api.ts` — add two typed wrappers

```ts
export async function scanSeries(dir: string): Promise<SeriesDescriptor[]>;
export async function loadSeries(dir: string, uid: string): Promise<{
    metadata: VolumeMetadata;
    bytes: ArrayBuffer;
}>;
```

`loadSeries` parses the frame header (first 4 bytes = metadata length, next N bytes = JSON, rest = voxels) and returns the voxel buffer as a `Uint8Array` / `ArrayBuffer` suitable for `Int16Array` construction.

### UI wiring

`DicomBrowser.svelte` and surrounding components replace the current "pick folder → load everything" flow with:

1. Pick folder → `scanSeries` → render a list of series with descriptions, keV (from `image_comments`), slice counts.
2. User clicks a series → `loadVolume` → cornerstone3D renders MPR.
3. Dual-energy / MMD: pick two series → `await Promise.all([loadSeries(lo), loadSeries(hi)])`.

### Deletions

- `src-tauri/src/commands/volume.rs::get_slice` — dead after rewrite.
- The 6-stage format-conversion pipeline in `volume.rs` (f32→i16→bytes→JSON→Uint8→Int16).
- Duplicate directory-walk logic in `dicom_loader::collect_slices` and `scan_dicom_series`.

## Data flow — timings

### Scenario A: Single-series load (300-slice CT, ~50 ms SMB latency)

| Time | Stage | What happens |
|---|---|---|
| 0.00 s | user action | Picks folder, frontend invokes `scan_series` |
| 0.02 s | scan | `tokio::fs::read_dir` — 1 SMB round-trip to list 384 entries |
| 0.42 s | scan | `buffer_unordered(48)` header reads; each reads ~4 KB (not 512 KB) |
| 0.43 s | scan | Group by UID, sort by IPP[2], return `Vec<SeriesDescriptor>` |
| 0.43 s | user action | Frontend renders picker, user clicks series, invokes `load_series` |
| 3.50 s | load | `rayon::par_iter` over 384 files: open + decode + rescale → `Vec<i16>` |
| 3.60 s | transfer | Framed `ArrayBuffer` returned; frontend parses metadata and builds `Int16Array` |
| 3.70 s | render | cornerstone3D constructs volume and renders first MPR frame |

**Total cold load: ~3.7 s.** Current: ~40–85 s. Speedup: **10–20×**.

### Scenario B: Dual-energy / MMD load

```js
const series = await scanSeries(dir);                         // ~0.4 s (shared)
const [low, high] = await Promise.all([
    loadSeries(dir, loSeries.uid),                             // ~3.3 s
    loadSeries(dir, hiSeries.uid),                             // ~3.3 s (parallel)
]);
```

**Total: ~3.7 s** — the two loads share the rayon pool and run in parallel. Current `load_dual_energy`: serial, ~20+ s.

### Scenario C: Cancellation

User picks series A, then picks series B mid-load. The Promise for A is dropped by the frontend; Tauri drops the future; `tokio` cleans up the task; `rayon` in-flight jobs finish naturally. No state leaks.

### Scenario D: First SMB touch after mount

First `read_dir` can stall 1–3 s for SMB negotiation. We surface a "Connecting to server…" message if `scan_series` has not produced a result within 500 ms.

### Bytes on the wire

| Call | Request | Response | Roundtrips |
|---|---|---|---|
| `scan_series` | JSON `{dir}` | JSON `Vec<SeriesDto>`, <10 KB | 1 |
| `load_series` | JSON `{dir, uid}` | Framed binary `ArrayBuffer`, ~200 MB for 300×512² | 1 |

**Total: 2 IPC roundtrips per patient load.** Current: ~390 (`list_patients` + `scan_series` + `load_dicom` + 384 × `get_slice`).

## Error handling

**Principle:** fail fast on structural errors (unrecoverable); skip-and-warn on per-file errors during scan (one bad DICOM in a folder shouldn't kill the whole scan); fail fast during load (can't have holes in a volume); structured errors returned to UI.

### Error taxonomy

| Failure | Layer | Strategy | User sees |
|---|---|---|---|
| Directory unreadable | scan | fail fast | "Can't read folder: <reason>" |
| SMB disconnected mid-scan | scan | fail fast with partial count | "Lost connection after N files" |
| All files fail header parse | scan | fail with hint | "No DICOM files found in folder" |
| Single header parse fails | scan | **skip + warn** | scan continues; UI shows `skipped: N` |
| SMB read error mid-load | load | fail fast | "Read error at slice N: <reason>" |
| Single file fails pixel decode | load | **fail fast** | "Failed to decode slice N" |
| Inconsistent dims | load | fail fast | "Slice N has mismatched dimensions" |
| Missing IPP | scan | fallback to `InstanceNumber`; warn | loads with warn indicator |
| Missing rescale tags | load | default slope=1, intercept=0; warn | loads with warn icon |
| Unsupported transfer syntax | load | fail fast with syntax name | "Unsupported compression: <UID>" |
| User cancels mid-load | load | drop future, silent | (silent) |
| Volume too large | load | pre-check dims × 2 bytes > budget | "Volume too large: N GB > limit" |

### Rust error type

```rust
#[derive(thiserror::Error, Debug)]
pub enum DicomLoadError {
    #[error("folder not readable: {0}")]
    IoError(#[from] std::io::Error),

    #[error("no DICOM files found ({scanned} scanned, {skipped} failed header parse)")]
    NoDicoms { scanned: usize, skipped: usize },

    #[error("series {uid} not found in folder")]
    SeriesNotFound { uid: String },

    #[error("slice {path:?} has dims {got:?}, expected {expected:?}")]
    InconsistentDims { path: PathBuf, got: [u32; 2], expected: [u32; 2] },

    #[error("unsupported transfer syntax: {0}")]
    UnsupportedTransferSyntax(String),

    #[error("volume too large: {requested_mb} MB exceeds {limit_mb} MB limit")]
    VolumeTooLarge { requested_mb: usize, limit_mb: usize },

    #[error("pixel decode failed at {path:?}: {source}")]
    DecodeFailed { path: PathBuf, source: dicom_pixeldata::Error },
}
```

Tauri commands map these to `Result<_, String>` via `Display`. Frontend renders the string directly.

### What is *not* an error

- A file without a `.dcm` extension — content-sniff for DICM magic at byte 128, not extension.
- A file that is not a DICOM at all (`.DS_Store`, `README.txt`) — skipped silently at magic check, not counted toward `skipped`.
- Non-image DICOMs (SR, KO, PR) — skipped at header parse (no `PixelData` tag), not logged as failure.

## Testing

### Unit tests (Rust)

In `crates/pcat-pipeline/src/dicom_loader.rs`:

- `test_read_header_skips_pixel_data` — verify `read_until(PIXEL_DATA)` does not read pixel bytes. Assert on byte count (via a wrapping `Read` adapter).
- `test_group_by_series` — mixed-UID header list produces correct buckets.
- `test_sort_by_ipp` — slices sorted by `ImagePositionPatient[2]`.
- `test_sort_fallback_to_instance_number` — IPP absent → sort by `InstanceNumber`.
- `test_missing_rescale_defaults_to_identity` — slope=1, intercept=0.
- `test_magic_byte_sniff` — non-DICOM file skipped without error.
- `test_inconsistent_dims_rejected` — load fails with `InconsistentDims`.

### Integration tests (Rust)

Fixture: 30-slice mini-CT committed at `crates/pcat-pipeline/tests/fixtures/mini-ct/` (~4 MB, generated from a synthetic phantom with `dcmtk`).

- `test_scan_fixture_returns_one_series` — scan end-to-end.
- `test_load_fixture_returns_correct_volume` — byte-level compare loaded `i16` against a reference `.raw` in the fixture.
- `test_load_with_corrupted_file_in_folder` — non-DICOM file present; verify silent skip.
- `test_dual_series_fixture` — two series in same folder, each loads independently.

### Benchmark harness

`cargo bench` (or `#[test]` guarded by `#[cfg(feature = "bench")]`) that:

- Accepts a folder via `BENCH_DICOM_DIR` env var (default: fixture).
- Times `scan_series`, `load_series`, total cold path.
- Compares against pre-rewrite code paths when invoked with `--baseline`.
- Not run in CI. Run manually pre-merge against a real SMB folder.

**Success criterion:** ≥ 10× speedup on real SMB data.

### Frontend test

One Playwright test at `tests/e2e/dicom-load.spec.ts`:

- Mocks Tauri `load_series` to return a synthetic 10-slice framed buffer.
- Verifies cornerstone3D renders and MPR viewport is non-blank.
- Verifies loading state is shown while the promise is pending.

### Not tested

- Pixel-perfect DICOM round-trip across every transfer syntax (trust `dicom-rs`).
- macOS SMB stack behavior (environmental, not our code).
- cornerstone3D MPR correctness (trust their tests).

### Manual verification before claiming done

1. Load a real patient folder from the SMB share. Measure wall clock. Confirm < 10 s for a 300-slice study.
2. Load a dual-energy folder. Both series load in parallel. MMD panel works.
3. Pull SMB cable mid-load. Graceful error; UI does not crash.
4. Open a folder with `.DS_Store` and `README.txt`. Silently skipped.
5. Open a MonoPlus folder. Confirm `image_comments` keV is captured, not `SeriesDescription` (per prior finding).

## Open questions

- Actual SMB latency on the lab network — the 40–85 s estimate is extrapolated from `WalkDir` serial cost at ~50 ms/op, not measured. Benchmark on day 1 to calibrate the target.
- Volume-too-large budget — pick a number that covers the largest expected cardiac CT (~1.5 GB for a 1000-slice 512² i16) but rejects accidental whole-study loads. Propose 4 GB soft limit.
- Whether to keep the existing `list_patients` command or fold it into a unified scan API. Tentative: keep for now; revisit after rewrite.

## Risks and mitigations

- **Risk:** `tauri::ipc::Response` binary return path misbehaves for very large (~1 GB) payloads. **Mitigation:** benchmark early. Fallback is to chunk into a handful of bulk calls (~4× 50 MB), still far fewer than the current 384.
- **Risk:** `dicom-rs`'s `read_until` interacts badly with some encapsulated transfer syntaxes. **Mitigation:** integration fixture covers at least one JPEG-2000-compressed series. If it breaks, the fallback is a two-pass design where we read full files in scan but still skip pixel decode — same speedup minus one factor.
- **Risk:** Cornerstone3D's `createLocalVolume` requires a specific orientation convention. **Mitigation:** cover with the Playwright test; the current v2 already does orientation math correctly, carry it over.

## Success metrics

- Cold load time for the reference 300-slice SMB study: **≤ 10 s** (stretch: ≤ 5 s).
- Total LOC delta: negative (fewer lines than current).
- IPC roundtrips per load: 2 (from ~390).
- No new persistent state on disk or in memory outside the active volume.
