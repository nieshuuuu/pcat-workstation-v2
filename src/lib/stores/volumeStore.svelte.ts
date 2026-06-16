/**
 * Reactive store for the currently loaded CT volume.
 *
 * Uses Svelte 5 runes (`$state`) for fine-grained reactivity.
 * Import `volumeStore` from any `.svelte` or `.svelte.ts` file.
 */

export type VolumeMetadata = {
  volumeId: string;
  /** Dimensions in Python/NumPy order: [Z, Y, X] (slices, rows, columns). */
  shape: [number, number, number];
  /** Voxel spacing in Python order: [sz, sy, sx]. */
  spacing: [number, number, number];
  origin: [number, number, number];
  /** Row-major 3x3 direction cosine matrix (9 elements). */
  direction: number[];
  windowCenter: number;
  windowWidth: number;
  patientName: string;
  studyDescription: string;
  dicomPath: string;
};

/** One volume resident in the Rust-side volume cache, usable as an
 *  immediate switch target without another decode. */
export type LoadedSeriesEntry = {
  name: string;
  path: string;
  uid: string;
  seriesDescription: string;
  kev: number | null;
  numSlices: number;
  rows: number;
  cols: number;
};

let currentVolume = $state<VolumeMetadata | null>(null);
let cornerstoneVolumeId = $state<string | null>(null);
let loading = $state(false);
let loadProgress = $state(0);
/** Optional human-readable detail shown alongside the progress bar
 *  (e.g. "Loading 3/6: MonoPlus_100keV"). Empty string hides the suffix. */
let loadMessage = $state('');
/** Volumes currently loaded in the Rust cache for the active patient. */
let loaded = $state<LoadedSeriesEntry[]>([]);

export const volumeStore = {
  get current() {
    return currentVolume;
  },
  get dicomPath(): string | null {
    return currentVolume?.dicomPath ?? null;
  },
  get cornerstoneVolumeId() {
    return cornerstoneVolumeId;
  },
  get loading() {
    return loading;
  },
  /** Progress percentage 0-100. */
  get loadProgress() {
    return loadProgress;
  },
  get loaded() {
    return loaded;
  },
  get loadMessage() {
    return loadMessage;
  },

  set(vol: VolumeMetadata) {
    currentVolume = vol;
  },
  setCornerstoneVolumeId(id: string) {
    cornerstoneVolumeId = id;
  },
  setLoading(v: boolean) {
    loading = v;
  },
  setLoadProgress(v: number) {
    loadProgress = v;
  },
  setLoadMessage(msg: string) {
    loadMessage = msg;
  },
  setLoaded(entries: LoadedSeriesEntry[]) {
    loaded = entries;
  },
  clear() {
    // Null the stored references only. Cornerstone's own LRU cache handles
    // eviction if memory pressure arises; keeping the volume cached lets the
    // A→B→A fast-reload path short-circuit via cache.getVolume(csId).
    currentVolume = null;
    cornerstoneVolumeId = null;
    loading = false;
    loadProgress = 0;
    loadMessage = '';
    loaded = [];
  },
};
