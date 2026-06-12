/**
 * Reactive store for the FAI analysis pipeline.
 *
 * Manages pipeline execution state, listens for Tauri progress events,
 * and stores per-vessel FAI statistics results.
 *
 * Uses Svelte 5 runes (`$state`) for fine-grained reactivity.
 */

import { invoke } from '@tauri-apps/api/core';
import { listen } from '@tauri-apps/api/event';
import { seedStore, type Vessel } from './seedStore.svelte';

// ---- Types ----

export type PipelineStatus = 'idle' | 'running' | 'complete' | 'error';

export type VesselProgress = {
  stage: string;
  progress: number; // 0.0 - 1.0
};

export type FaiStats = {
  vessel: string;
  n_voi_voxels: number;
  n_fat_voxels: number;
  fat_fraction: number;
  hu_mean: number;
  hu_std: number;
  hu_median: number;
  fai_risk: string; // "HIGH" or "LOW"
  histogram_bins: number[];
  histogram_counts: number[];
  radial_profile: {
    distances_mm: number[];
    mean_hu: number[];
    std_hu: number[];
  } | null;
  angular_asymmetry: {
    sectors: {
      label: string;
      angle_deg: number;
      hu_mean: number;
      hu_std: number;
      n_voxels: number;
      fai_risk: string;
    }[];
  } | null;
};

type PipelineProgressPayload = {
  vessel: string;
  stage: string;
  progress: number;
};

// Seed payload format expected by the Rust backend
type SeedPayload = {
  ostium_mm: [number, number, number];
  waypoints_mm: [number, number, number][];
  segment_start_mm: number;
  segment_length_mm: number;
};

// ---- Reactive state ----

let status = $state<PipelineStatus>('idle');
let progress = $state<Record<string, VesselProgress>>({});
let results = $state<Record<string, FaiStats> | null>(null);
let errorMsg = $state('');

// ---- Helpers ----

const ALL_VESSELS: Vessel[] = ['LAD', 'LCx', 'RCA'];

/**
 * Build the seeds payload from the seed store.
 * Only includes vessels that have at least an ostium + 1 waypoint.
 */
/** Swap cornerstone [x,y,z] to pipeline [z,y,x] ordering. */
function toZyx(p: [number, number, number]): [number, number, number] {
  return [p[2], p[1], p[0]];
}

function buildSeedsPayload(): Record<string, SeedPayload> {
  const seeds: Record<string, SeedPayload> = {};

  for (const vessel of ALL_VESSELS) {
    const data = seedStore.vessels[vessel];
    if (data.seeds.length < 2) continue;

    // Use ostiumFraction position if set, otherwise fall back to first seed
    const ostiumXyz = seedStore.getOstiumWorldPosForVessel(vessel) ?? data.seeds[0].position;

    seeds[vessel] = {
      ostium_mm: toZyx(ostiumXyz),
      waypoints_mm: data.seeds.map((s) => toZyx(s.position)),
      segment_start_mm: 0,
      segment_length_mm: 40,
    };
  }

  return seeds;
}

// ---- Store ----

export const pipelineStore = {
  // --- Getters ---

  get status(): PipelineStatus {
    return status;
  },
  get progress(): Record<string, VesselProgress> {
    return progress;
  },
  get results(): Record<string, FaiStats> | null {
    return results;
  },
  get error(): string {
    return errorMsg;
  },

  /** True if at least one vessel has enough seeds (2+) to run the pipeline. */
  get canRun(): boolean {
    return ALL_VESSELS.some((v) => seedStore.vessels[v].seeds.length >= 2);
  },

  // --- Actions ---

  async run() {
    if (status === 'running') return;

    // Reset state
    status = 'running';
    errorMsg = '';
    results = null;
    progress = {};

    // Initialize progress for each vessel that will be analyzed
    const seeds = buildSeedsPayload();
    for (const vessel of Object.keys(seeds)) {
      progress[vessel] = { stage: 'Queued', progress: 0 };
    }

    // Listen for progress events
    let unlisten: (() => void) | undefined;

    try {
      unlisten = await listen<PipelineProgressPayload>(
        'pipeline-progress',
        (event) => {
          const { vessel, stage, progress: prog } = event.payload;
          // Create a new object to trigger Svelte reactivity
          progress = {
            ...progress,
            [vessel]: { stage, progress: prog },
          };
        },
      );

      // Invoke the Rust pipeline command
      const result = await invoke<Record<string, FaiStats>>('run_pipeline', {
        seeds,
      });

      // Store results and mark complete
      results = result;
      status = 'complete';
    } catch (e) {
      errorMsg = e instanceof Error ? e.message : String(e);
      status = 'error';
      console.error('Pipeline failed:', e);
    } finally {
      // Clean up event listener
      unlisten?.();
    }
  },

  /** Reset all pipeline state back to idle. */
  reset() {
    status = 'idle';
    progress = {};
    results = null;
    errorMsg = '';
  },

  /** Restore saved pipeline results (from loaded file). */
  restoreResults(savedResults: Record<string, FaiStats>) {
    results = savedResults;
    status = 'complete';
    errorMsg = '';
  },
};
