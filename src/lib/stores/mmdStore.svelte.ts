/**
 * Shared store for the water/lipid (MMD) analysis results.
 *
 * The results live here, not in MmdAnalysisView's local state, so the unified
 * session save/load ($lib/session) can snapshot and restore them from anywhere
 * — the same way pipelineStore holds the FAI results. The view writes its
 * computed summary/surfaces here; SurfacePlotPanel and the session reader them.
 */

import type { MmdSummary, CrossSectionSurface } from '$lib/api';

export type MmdSnapshot = {
  summary: MmdSummary | null;
  surfaces: CrossSectionSurface[];
};

let summary = $state<MmdSummary | null>(null);
let surfaces = $state<CrossSectionSurface[]>([]);

export const mmdStore = {
  get summary(): MmdSummary | null {
    return summary;
  },
  set summary(v: MmdSummary | null) {
    summary = v;
  },
  get surfaces(): CrossSectionSurface[] {
    return surfaces;
  },
  set surfaces(v: CrossSectionSurface[]) {
    surfaces = v;
  },

  /** Clear results (new patient / re-open before a fresh run). */
  clear() {
    summary = null;
    surfaces = [];
  },

  /** Restore a saved snapshot (from the unified session loader). */
  restore(snap: Partial<MmdSnapshot> | null | undefined) {
    if (!snap) return;
    summary = snap.summary ?? null;
    surfaces = snap.surfaces ?? [];
  },
};
