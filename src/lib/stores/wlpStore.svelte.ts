/**
 * Holds the whole-volume water/lipid/protein (poly2 surface) model so it can be
 * snapshotted and restored by the unified session loader ($lib/session).
 * WaterLipidView is the producer (sets it on Run). The model is the single
 * source of truth — per-slice f_w/f_l/f_p/σ_f maps are derived from it on demand
 * by `getWlpSlice`. Mirrors wlStore, for the 3-material path.
 */
import type { WlpModel } from '$lib/api';

let model = $state<WlpModel | null>(null);

export const wlpStore = {
  get model(): WlpModel | null {
    return model;
  },
  /** Record the model produced by a Run. */
  set(m: WlpModel | null) {
    model = m;
  },
  /** Clear (new patient load / failed run). */
  reset() {
    model = null;
  },
  /** Restore a saved model (from the unified session loader). */
  restore(m: WlpModel | null | undefined) {
    model = m ?? null;
  },
};
