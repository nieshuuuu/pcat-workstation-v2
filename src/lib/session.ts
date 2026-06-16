/**
 * Unified per-patient session save/load — "save everything" in one file.
 *
 * Bundles the work that previously lived in separate, fragmented saves:
 *   - seeds + ostium            (seedStore)
 *   - FAI analysis results      (pipelineStore)
 *   - water/lipid quantification: summary + 3D surfaces (mmdStore)
 *
 * The backend stores the opaque JSON blob; this module owns the schema, so
 * adding a field here never needs a Rust change. One saveSession() captures
 * everything; one loadSession() restores it into the stores.
 */

import { saveSession as saveSessionCmd, loadSession as loadSessionCmd, restoreWlCalibration } from '$lib/api';
import { seedStore } from '$lib/stores/seedStore.svelte';
import { pipelineStore } from '$lib/stores/pipelineStore.svelte';
import { mmdStore } from '$lib/stores/mmdStore.svelte';
import { wlStore } from '$lib/stores/wlStore.svelte';

const SESSION_VERSION = 1;

export type SessionMeta = { savedAt: string };

/** Persist seeds + FAI + water/lipid quantification for this patient. */
export async function saveSession(dicomPath: string): Promise<void> {
  if (!dicomPath) throw new Error('no patient loaded');
  const bundle = {
    version: SESSION_VERSION,
    savedAt: new Date().toISOString(),
    seeds: JSON.parse(seedStore.exportJson()),
    fai: pipelineStore.results,
    mmd: { summary: mmdStore.summary, surfaces: mmdStore.surfaces },
    // Whole-volume water/lipid calibration (the SSoT; f_w/σ_f maps derive from it).
    wl: wlStore.calibration,
  };
  await saveSessionCmd(dicomPath, JSON.stringify(bundle));
}

/** Restore everything saved by saveSession into the stores. Returns the save
 *  metadata, or null if no session exists for this patient. */
export async function loadSession(dicomPath: string): Promise<SessionMeta | null> {
  if (!dicomPath) return null;
  const json = await loadSessionCmd(dicomPath);
  if (!json) return null;
  const bundle = JSON.parse(json);
  if (bundle.seeds) seedStore.importJson(JSON.stringify(bundle.seeds));
  if (bundle.fai) pipelineStore.restoreResults(bundle.fai);
  if (bundle.mmd) mmdStore.restore(bundle.mmd);
  if (bundle.wl) {
    wlStore.restore(bundle.wl);
    // Push the calibration back into backend state so getWlSlice can derive maps.
    // Best-effort: requires the dual-energy volume to be loaded (it is after a
    // patient load with a keV pair); a single-series open has no WL to restore.
    try { await restoreWlCalibration(bundle.wl); } catch { /* dual-energy not loaded */ }
  } else {
    wlStore.reset();
  }
  return { savedAt: bundle.savedAt ?? '' };
}
