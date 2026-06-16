/**
 * Holds the whole-volume water/lipid (GLS) calibration so it can be snapshotted
 * and restored by the unified session loader ($lib/session). WaterLipidView is
 * the producer (sets it on Run); SurfacePlotPanel-style consumers and the
 * session reader it. The calibration is the single source of truth — per-slice
 * f_w/σ_f maps are derived from it on demand by `getWlSlice`.
 */
import type { WlCalibration } from '$lib/api';

let calibration = $state<WlCalibration | null>(null);

export const wlStore = {
  get calibration(): WlCalibration | null {
    return calibration;
  },
  /** Record the calibration produced by a Run. */
  set(c: WlCalibration | null) {
    calibration = c;
  },
  /** Clear (new patient load / failed run). */
  reset() {
    calibration = null;
  },
  /** Restore a saved calibration (from the unified session loader). */
  restore(c: WlCalibration | null | undefined) {
    calibration = c ?? null;
  },
};
