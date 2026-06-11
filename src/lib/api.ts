/**
 * API layer — Tauri IPC commands to the Rust backend.
 * No HTTP, no Python sidecar. All calls go through `invoke()`.
 */
import { invoke } from '@tauri-apps/api/core';

/** Open native folder picker. Returns path or null if cancelled. */
export async function openDicomDialog(): Promise<string | null> {
  return invoke<string | null>('open_dicom_dialog');
}

/** Get list of recently opened DICOM folder paths. */
export async function getRecentDicoms(): Promise<string[]> {
  return invoke<string[]>('get_recent_dicoms');
}

/** Save seeds JSON to app data directory, keyed by DICOM path. Returns the file path. */
export async function saveSeeds(seedsJson: string, dicomPath: string): Promise<string> {
  return invoke<string>('save_seeds', { seedsJson, dicomPath });
}

/** Load seeds JSON from app data directory, keyed by DICOM path. Returns null if no file. */
export async function loadSeeds(dicomPath: string): Promise<string | null> {
  return invoke<string | null>('load_seeds', { dicomPath });
}

/* ── DICOM scan & bulk load (fast path) ────────────────────── */

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
  /** ImagePositionPatient of first slice in patient LPS mm: [x, y, z]. */
  image_position_patient: [number, number, number];
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
  /** ImagePositionPatient of first slice in patient LPS mm: [x, y, z]. */
  image_position_patient: [number, number, number];
}

export type DicomLoadPhase =
  | 'scanning'
  | 'scanned'
  | 'decoding'
  | 'patient_series'
  | 'done';

export interface DicomLoadProgress {
  phase: DicomLoadPhase;
  done: number;
  total: number;
  /** Optional detail string — e.g. the current series folder name during a
   *  patient-wide load. Omitted for phases that don't carry one. */
  detail?: string;
}

/** Scan a DICOM folder (header-only). Header-only; no pixel data decoded. */
export async function scanSeries(dir: string): Promise<SeriesDescriptor[]> {
  return invoke<SeriesDescriptor[]>('scan_series', { path: dir });
}

/**
 * Load one series as a single bulk transfer: returns the parsed metadata and
 * a view over the decoded i16 HU voxels (z-major).
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

  const voxelOffset = 4 + metaLen;
  let voxels: Int16Array;
  if (voxelOffset % 2 === 0) {
    voxels = new Int16Array(buf, voxelOffset);
  } else {
    // Int16Array requires 2-byte alignment. If metaLen is odd, copy the bytes.
    const copy = new Uint8Array(buf.byteLength - voxelOffset);
    copy.set(new Uint8Array(buf, voxelOffset));
    voxels = new Int16Array(copy.buffer);
  }
  return { metadata, voxels };
}

/**
 * Load low-energy + high-energy DICOM series in parallel for dual-energy MMD.
 *
 * Folder names MUST contain a keV label (e.g. `MonoPlus_70keV`,
 * `MonoPlus_150keV`) — lab-internal data has `ImageComments` stripped and
 * `SeriesDescription` mislabeled, so the folder name is the only reliable
 * source of keV.
 *
 * Populates both `state.dual_energy` (for MMD) and `state.volume` (low),
 * and returns the low-energy framed bundle so the caller can build the
 * cornerstone volume exactly as with `loadSeries`.
 */
export type LoadedSeriesDescriptor = {
  name: string;
  path: string;
  uid: string;
  series_description: string;
  kev: number | null;
  num_slices: number;
  rows: number;
  cols: number;
};

export type PatientLoadResult = {
  series: LoadedSeriesDescriptor[];
  active_index: number;
  failures: string[];
};

/** Load every DICOM series under the patient directory into the Rust-side
 *  volume cache. After this returns, `setActiveVolume` can switch between
 *  them instantly without redecoding. Auto-pairs MonoPlus keV series for MMD. */
export async function loadPatientAll(
  patientDir: string,
): Promise<PatientLoadResult> {
  return invoke<PatientLoadResult>('load_patient_all', { patientDir });
}

/** Switch the active volume to a series previously loaded into the cache.
 *  Returns the framed binary bundle so the caller can rebuild the
 *  cornerstone3D volume exactly as with `loadSeries`. */
export async function setActiveVolume(
  dir: string,
  uid: string,
): Promise<{ metadata: VolumeMetadata; voxels: Int16Array }> {
  const buf = await invoke<ArrayBuffer>('set_active_volume', { dir, uid });
  const view = new DataView(buf);
  const metaLen = view.getUint32(0, true);
  const metaBytes = new Uint8Array(buf, 4, metaLen);
  const metadata = JSON.parse(new TextDecoder().decode(metaBytes)) as VolumeMetadata;

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

export async function loadDualEnergy(
  lowDir: string,
  highDir: string,
): Promise<{ metadata: VolumeMetadata; voxels: Int16Array }> {
  const buf = await invoke<ArrayBuffer>('load_dual_energy', { lowDir, highDir });
  const view = new DataView(buf);
  const metaLen = view.getUint32(0, true);
  const metaBytes = new Uint8Array(buf, 4, metaLen);
  const metadata = JSON.parse(new TextDecoder().decode(metaBytes)) as VolumeMetadata;

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

/**
 * Query whether the Rust AppState currently holds the volume at (dir, uid).
 * If so, returns the cached pipeline metadata (same shape as `loadSeries`'s
 * `metadata` tuple field) so the frontend can skip the heavy load. Returns
 * null when the Rust side doesn't currently hold this volume.
 */
export async function reuseLoadedVolume(
  dir: string,
  uid: string,
): Promise<VolumeMetadata | null> {
  return invoke<VolumeMetadata | null>('reuse_loaded_volume', { dir, uid });
}

/**
 * Subscribe to DICOM load progress events emitted by scan_series and
 * load_series. Returns an unlisten function — call it when the consumer
 * is finished to avoid leaks.
 */
export async function onDicomLoadProgress(
  handler: (progress: DicomLoadProgress) => void,
): Promise<() => void> {
  const { listen } = await import('@tauri-apps/api/event');
  const unlisten = await listen<DicomLoadProgress>('dicom_load_progress', (e) => {
    handler(e.payload);
  });
  return unlisten;
}

/* ── Annotation + Snake commands ─────────────────────────── */

export type AnnotationTarget = {
  image: number[];         // HU values, row-major, pixels*pixels
  pixels: number;
  width_mm: number;
  arc_mm: number;
  frame_index: number;
  vessel_wall: [number, number][];  // [x,y] pixel coords
  vessel_radius_mm: number;
  init_boundary: [number, number][]; // [x,y] pixel coords
};

export type SnakeResult = {
  points: [number, number][];
  iterations: number;
  max_displacement: number;
  converged: boolean;
};

export type MmdSummary = {
  method: string;
  iterations: number;
  converged: boolean;
  n_voxels: number;
  mean_water_frac: number;
  mean_lipid_frac: number;
};

/** Generate annotation targets for all cross-section frames along a centerline.
 *
 * `ostiumMm` (optional, in `[z, y, x]` pipeline order) shifts the first
 * cross-section to start at the coronary ostium rather than the first
 * centerline waypoint. Pass `null` to keep legacy behaviour. */
export async function generateAnnotationTargets(
  centerlineMm: [number, number, number][],
  ostiumMm: [number, number, number] | null = null,
): Promise<AnnotationTarget[]> {
  return invoke<AnnotationTarget[]>('generate_annotation_targets', {
    centerlineMm,
    ostiumMm,
  });
}

/** Initialize a circular snake contour for a given annotation target. */
export async function initSnake(
  targetIndex: number,
  initRadiusMm?: number,
): Promise<[number, number][]> {
  return invoke<[number, number][]>('init_snake', {
    targetIndex,
    initRadiusMm: initRadiusMm ?? null,
  });
}

/** Evolve the active contour (snake) for a target. */
export async function evolveSnake(
  targetIndex: number,
  nIterations: number = 200,
): Promise<SnakeResult> {
  return invoke<SnakeResult>('evolve_snake', {
    targetIndex,
    nIterations,
  });
}

/** Replace the snake control points for a target (after manual drag). */
export async function updateSnakePoints(
  targetIndex: number,
  points: [number, number][],
): Promise<void> {
  return invoke<void>('update_snake_points', {
    targetIndex,
    points,
  });
}

/** Insert a new control point on the snake contour at a given position. */
export async function addSnakePoint(
  targetIndex: number,
  position: [number, number],
): Promise<number> {
  return invoke<number>('add_snake_point', {
    targetIndex,
    position,
  });
}

/** Finalize the contour for a target (marks it as done). */
export async function finalizeContour(
  targetIndex: number,
): Promise<void> {
  return invoke<void>('finalize_contour', {
    targetIndex,
  });
}

export type AdoptedContour = {
  target_index: number;
  points: [number, number][];
};

/** Adopt the auto-detected vessel wall as the finalized contour, skipping
 *  the snake evolve step. Pass `all: true` to finalize every target whose
 *  vessel wall is non-empty. Returns the resampled contours the backend
 *  now holds so the UI stays in sync. */
export async function useVesselWallAsContour(
  opts: { targetIndex?: number; all?: boolean } = {},
): Promise<AdoptedContour[]> {
  return invoke<AdoptedContour[]>('use_vessel_wall_as_contour', {
    targetIndex: opts.targetIndex ?? null,
    all: opts.all ?? false,
  });
}

/** Run noise-aware water/lipid (GLS) decomposition on the annotated ROI.
 *  'gls' is the only supported method (the 3-material solvers were dropped). */
export async function runMmdOnRoi(
  method: string = 'gls',
): Promise<MmdSummary> {
  return invoke<MmdSummary>('run_mmd_on_roi', { method });
}
/* ── Surface sampling + MMD overlay ─────────────────────── */

export type CrossSectionSurface = {
  arc_mm: number;
  theta_deg: number[];
  r_mm: number[];
  surface: number[];
  n_theta: number;
  n_radial: number;
  max_r_per_theta: number[];
};

/** Sample radial-angular surface data from MMD result for all finalized cross-sections. */
export async function sampleSurfaces(
  material: string,
  unit: string,
): Promise<CrossSectionSurface[]> {
  return invoke<CrossSectionSurface[]>('sample_surfaces', { material, unit });
}

/** Get MMD material overlay for a single cross-section (flat f32 array). */
export async function getMmdOverlay(
  targetIndex: number,
  material: string,
  unit: string,
): Promise<number[]> {
  return invoke<number[]>('get_mmd_overlay', { targetIndex, material, unit });
}

/* ── Whole-volume noise-aware GLS water/lipid decomposition ─────────────────
 *
 * The latest method (ported from wl-noise-aware-mmd): self-calibrated, runs on
 * the whole dual-energy volume, renders axial f_w/f_l maps like
 * fw_fl_maps_theolipid_57955439.png. Calibration is the single source of truth;
 * slices are derived on demand by `getWlSlice`. */

export type WlAnchor = 'adipose' | 'theoretical';

/** Self-measured calibration returned by `runWaterLipid`. Mirrors the Rust
 *  `WlCalibration` struct (snake_case fields). All endpoints/noise lines are
 *  measured from the patient's own fat + muscle — no external phantom. */
export type WlCalibration = {
  low_kev: number;
  high_kev: number;
  hu_w: [number, number];
  hu_l_adipose: [number, number];
  hu_l_theo: [number, number];
  hu_mus: [number, number];
  ab_low: [number, number];
  ab_high: [number, number];
  sigma_fat: [number, number];
  sigma_mus: [number, number];
  rho: number;
  gate: [number, number];
  n_fat: number;
  n_mus: number;
  fraction_kept: number;
  /** [nz, ny, nx] of the dual-energy grid. */
  dims: [number, number, number];
};

/** One decoded axial slice: CT HU (i16) + water fraction + σ_f (both f32, NaN
 *  outside the soft-tissue gate). f_l = 1 − f_w is derived on the client. */
export type WlSlice = {
  z: number;
  ny: number;
  nx: number;
  anchor: WlAnchor;
  ct: Int16Array;
  fw: Float32Array;
  sf: Float32Array;
};

/** Self-calibrate + register the whole-volume GLS water/lipid solver for the
 *  loaded dual-energy volume. Returns the measured calibration. Run once. */
export async function runWaterLipid(): Promise<WlCalibration> {
  return invoke<WlCalibration>('run_water_lipid');
}

/** Derive one axial slice's CT + f_w + σ_f maps from the stored calibration. */
export async function getWlSlice(z: number, anchor: WlAnchor): Promise<WlSlice> {
  const buf = await invoke<ArrayBuffer>('get_wl_slice', { z, anchor });
  const view = new DataView(buf);
  const metaLen = view.getUint32(0, true);
  const metaBytes = new Uint8Array(buf, 4, metaLen);
  const meta = JSON.parse(new TextDecoder().decode(metaBytes)) as {
    z: number;
    ny: number;
    nx: number;
    anchor: WlAnchor;
  };

  const plane = meta.ny * meta.nx;
  let off = 4 + metaLen;
  // Copy each region into its own buffer — guarantees the typed-array
  // alignment Int16Array/Float32Array require regardless of the JSON length.
  const ct = new Int16Array(buf.slice(off, off + plane * 2));
  off += plane * 2;
  const fw = new Float32Array(buf.slice(off, off + plane * 4));
  off += plane * 4;
  const sf = new Float32Array(buf.slice(off, off + plane * 4));

  return { z: meta.z, ny: meta.ny, nx: meta.nx, anchor: meta.anchor, ct, fw, sf };
}

/* ── Save/Load annotations + CSV export ───────────────── */

export type AnnotationStateJson = {
  centerline_mm: [number, number, number][];
  snake_contours: Record<number, [number, number][]>;
  finalized: Record<number, boolean>;
  mmd_method: string | null;
  mmd_iterations: number | null;
  mmd_converged: boolean | null;
};

/** Save the current annotation state for the given patient. Returns the file path. */
export async function saveAnnotations(
  dicomPath: string,
  centerlineMm: [number, number, number][],
): Promise<string> {
  return invoke<string>('save_annotations', { dicomPath, centerlineMm });
}

/** Load saved annotation state for the given patient. Returns null if no save exists. */
export async function loadAnnotations(
  dicomPath: string,
): Promise<AnnotationStateJson | null> {
  return invoke<AnnotationStateJson | null>('load_annotations', { dicomPath });
}

/** Export current MMD surface data as a CSV string. */
export async function exportMmdCsv(
  patientId: string,
): Promise<string> {
  return invoke<string>('export_mmd_csv', { patientId });
}

/* ── Patient browser ──────────────────────────────────── */

export type PatientStatus = 'not_started' | 'in_progress' | 'complete';

export type PatientInfo = {
  /** Folder name (stable patient ID). */
  id: string;
  /** Absolute path to the patient's DICOM folder. */
  path: string;
  status: PatientStatus;
  /** Number of cross-sections marked finalized in saved annotations. */
  finalized_count: number;
  /** Whether MMD has been run and stored in saved annotations. */
  has_mmd: boolean;
};

/**
 * Walk `rootDir` for patient subfolders and return a sorted list with status
 * badges derived from each patient's saved annotation JSON.
 */
export async function listPatients(rootDir: string): Promise<PatientInfo[]> {
  return invoke<PatientInfo[]>('list_patients', { rootDir });
}

export type SeriesDirInfo = {
  /** Folder name (e.g. `MonoPlus_70keV`). */
  name: string;
  /** Absolute path to the series folder. */
  path: string;
  /** File count in that folder (≈ DICOM slices). */
  num_files: number;
};

/**
 * List immediate subdirectories of a patient folder. Each subdirectory is
 * typically one DICOM series. No DICOM headers are parsed — fast.
 */
export async function listSeriesDirs(patientPath: string): Promise<SeriesDirInfo[]> {
  return invoke<SeriesDirInfo[]>('list_series_dirs', { patientPath });
}
