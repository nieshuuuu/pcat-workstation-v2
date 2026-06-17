/**
 * Holds the currently-loaded patient's data-quality flag. The footer reads it;
 * the flag dialog sets it. Patient identity is NOT duplicated here — the caller
 * passes the series `dicomPath` (volumeStore.dicomPath) and the backend derives
 * the patient folder. `flags/{patient}.json` on disk is the single source of
 * truth; this store is the in-memory view for the loaded patient.
 */
import type { PatientFlag } from '$lib/api';
import { getPatientFlag, setPatientFlag } from '$lib/api';

let current = $state<PatientFlag | null>(null);

export const flagStore = {
  get current(): PatientFlag | null {
    return current;
  },
  /** Load the flag for the loaded patient (or clear if none / no patient). */
  async load(dicomPath: string | null): Promise<void> {
    if (!dicomPath) {
      current = null;
      return;
    }
    try {
      current = await getPatientFlag(dicomPath);
    } catch (e) {
      console.error('load patient flag failed:', e);
      current = null;
    }
  },
  /** Set (flagged=true) or clear (flagged=false) the flag and update the view. */
  async set(dicomPath: string, flagged: boolean, note: string): Promise<void> {
    const flaggedAt = new Date().toISOString();
    await setPatientFlag(dicomPath, flagged, note, flaggedAt);
    current = flagged ? { flagged: true, note, flagged_at: flaggedAt } : null;
  },
};
