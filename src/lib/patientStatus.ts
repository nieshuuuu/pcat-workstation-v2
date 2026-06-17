/**
 * Live workflow-status derivation for the footer. This MIRRORS the canonical
 * Rust rule (src-tauri/src/commands/dicom.rs `status_for` + `session_flags`,
 * which is unit-tested). Kept trivially simple so the two cannot drift. The Rust
 * version drives the patient-list badges (from disk); this drives the footer of
 * the currently-loaded patient (live, from the in-memory stores).
 */
export type WorkflowStatus = 'not_started' | 'in_progress' | 'complete';

export function derivePatientStatus(f: {
  hasSeeds: boolean;
  hasFai: boolean;
  hasMmd: boolean;
}): WorkflowStatus {
  if (f.hasFai && f.hasMmd) return 'complete';
  if (f.hasSeeds) return 'in_progress';
  return 'not_started';
}

/** Patient id (folder name) for display, from a *series* dicomPath: the patient
 *  folder is the series' parent, so the id is the second-to-last path segment. */
export function patientIdOf(dicomPath: string | null): string {
  if (!dicomPath) return '';
  const parts = dicomPath.replace(/\/+$/, '').split('/');
  return parts.length >= 2 ? parts[parts.length - 2] : (parts[parts.length - 1] ?? '');
}
