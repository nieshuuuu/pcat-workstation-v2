# Patient Status Display & Data-Quality Flagging — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Show which patient is loaded + its workflow status in the footer, highlight the loaded patient in the list, add a patient-level "Flagged" data-quality mark, and fix the pre-existing bug where `complete` never auto-shows.

**Architecture:** The only new *stored* fact is the human "Flagged" judgement (a per-patient sidecar `flags/{key}.json`). Everything else is *derived*: the patient list's status is derived on the Rust side from the session bundle the frontend actually writes (fixing a producer/consumer Single-Source-of-Truth violation), and the footer's live status is derived from the in-memory stores. Auto-save on FAI/MMD completion makes the derived disk status appear without a manual Save.

**Tech Stack:** Tauri 2 (Rust) backend, Svelte 5 runes frontend, serde_json sidecar persistence.

## Global Constraints

- **SSoT / derive-don't-sync:** never store a value that can be derived; never write the same fact to two files and "keep them in sync". (This bug was exactly that.)
- **YAGNI:** no flag categories, no series-level flag, no manual "mark complete", no flag-from-list. Patient-level boolean flag + optional note only.
- **`dicomPath` is a SERIES folder path**, not the patient folder. Patient folder = its parent dir. Patient id = `basename(parent(dicomPath))`.
- **Flag scope is patient-level.** Key flags by `patient_file_key(parent_of(dicomPath))` so every series of one patient shares one flag file.
- **`complete` = FAI results present AND MMD summary present** (in the same session bundle).
- **UI label for the flag is exactly `Flagged`** (red, with a `⚠` icon).
- **No new dependencies.** No frontend test runner exists — automated tests are Rust-side (`#[cfg(test)]` in the relevant file); frontend changes are verified by `npm run build` + manual runtime check.
- **Fail loud:** no silent `try/except: pass`. Auto-save failures are `console.error`-logged (acceptable for a background convenience save; the manual Save remains).
- **Rust test command:** `cargo test --manifest-path src-tauri/Cargo.toml --lib <filter>` (the src-tauri crate builds as a lib; existing unit tests run this way).

---

### Task 1: Rust — fix patient status derivation (read `sessions/`, `complete` = FAI && MMD)

Replaces the broken `read_annotation_summary` path (reads `annotations/{sanitize(folder_name)}.json` — a key nothing writes) with derivation from the per-series session bundles that `save_session` actually writes. Adds two pure, unit-tested functions.

**Files:**
- Modify: `src-tauri/src/commands/dicom.rs` (status logic in `list_patients` ~223-278; remove now-dead `read_annotation_summary` ~196-217; add pure helpers + tests in the existing `#[cfg(test)] mod tests` ~1323)

**Interfaces:**
- Produces: `fn session_flags(bundle: &serde_json::Value) -> (bool, bool)` — `(has_fai, has_mmd)`; `fn status_for(complete: bool, has_seeds: bool) -> &'static str`; `fn patient_session_summary(app: &tauri::AppHandle, patient_path: &str) -> (bool, bool)` — `(complete, has_mmd)`.
- Consumes: existing `sanitize_for_filename`, `patient_has_seeds`, `PatientInfo`.

- [ ] **Step 1: Write the failing tests**

Add to the `#[cfg(test)] mod tests` block in `src-tauri/src/commands/dicom.rs` (after the existing `patient_file_key_unique_across_patients_with_same_series_name` test):

```rust
    #[test]
    fn session_flags_complete_bundle() {
        // Mirrors src/lib/session.ts saveSession bundle shape.
        let bundle: serde_json::Value = serde_json::from_str(
            r#"{ "version":1, "savedAt":"t",
                 "seeds":{"activeVessel":"RCA","vessels":{}},
                 "fai":{"RCA":{"fai_mean_hu":-75.0}},
                 "mmd":{"summary":{"n":3},"surfaces":[]},
                 "wl":null }"#,
        )
        .unwrap();
        assert_eq!(session_flags(&bundle), (true, true));
    }

    #[test]
    fn session_flags_fai_only_is_not_complete() {
        let bundle: serde_json::Value =
            serde_json::from_str(r#"{"fai":{"RCA":{}},"mmd":{"summary":null}}"#).unwrap();
        assert_eq!(session_flags(&bundle), (true, false));
    }

    #[test]
    fn session_flags_empty_bundle() {
        let bundle: serde_json::Value =
            serde_json::from_str(r#"{"fai":null,"mmd":{"summary":null}}"#).unwrap();
        assert_eq!(session_flags(&bundle), (false, false));
    }

    #[test]
    fn status_for_rules() {
        // (complete, has_seeds) -> label. This is the bug regression: a bundle
        // with FAI+MMD must yield "complete".
        assert_eq!(status_for(true, true), "complete");
        assert_eq!(status_for(true, false), "complete");
        assert_eq!(status_for(false, true), "in_progress");
        assert_eq!(status_for(false, false), "not_started");
    }
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cargo test --manifest-path src-tauri/Cargo.toml --lib session_flags status_for`
Expected: FAIL — `cannot find function 'session_flags'` / `'status_for'`.

- [ ] **Step 3: Add the pure helpers**

In `src-tauri/src/commands/dicom.rs`, replace the entire `read_annotation_summary` function (lines ~196-217) with these three functions:

```rust
/// Whether a parsed session bundle contains FAI results and an MMD summary.
/// Contract with src/lib/session.ts `saveSession`: the bundle has
/// `fai: object|null` (a non-empty object means FAI ran) and
/// `mmd: { summary: any|null }` (non-null summary means MMD ran). These three
/// field names are the contract between session.ts and this reader; the tests
/// above pin it.
fn session_flags(bundle: &serde_json::Value) -> (bool, bool) {
    let has_fai = bundle
        .get("fai")
        .and_then(|v| v.as_object())
        .map(|o| !o.is_empty())
        .unwrap_or(false);
    let has_mmd = bundle
        .get("mmd")
        .and_then(|m| m.get("summary"))
        .map(|s| !s.is_null())
        .unwrap_or(false);
    (has_fai, has_mmd)
}

/// Derived patient status label. Mirrored (kept trivially simple) in
/// src/lib/patientStatus.ts `derivePatientStatus` for the live footer; this Rust
/// version is the canonical, tested one driving the patient-list badges.
fn status_for(complete: bool, has_seeds: bool) -> &'static str {
    if complete {
        "complete"
    } else if has_seeds {
        "in_progress"
    } else {
        "not_started"
    }
}

/// Scan saved session bundles for this patient. Sessions are keyed per *series*
/// (`patient_file_key` of a series subfolder) and a patient folder holds several
/// series, so we substring-match the sanitized patient path exactly like
/// `patient_has_seeds`. Returns (complete, has_mmd): `complete` iff ANY single
/// session has both FAI and MMD (that series is fully analyzed).
fn patient_session_summary(app: &tauri::AppHandle, patient_path: &str) -> (bool, bool) {
    let dir = app
        .path()
        .app_data_dir()
        .expect("app data dir")
        .join("sessions");
    let needle = sanitize_for_filename(patient_path);
    if needle.is_empty() {
        return (false, false);
    }
    let Ok(read) = std::fs::read_dir(&dir) else {
        return (false, false);
    };
    let mut complete = false;
    let mut has_mmd = false;
    for entry in read.flatten() {
        let name = entry.file_name().to_string_lossy().to_string();
        if !name.contains(&needle) {
            continue;
        }
        let Ok(data) = std::fs::read_to_string(entry.path()) else {
            continue;
        };
        let Ok(json) = serde_json::from_str::<serde_json::Value>(&data) else {
            continue;
        };
        let (f, m) = session_flags(&json);
        complete |= f && m;
        has_mmd |= m;
    }
    (complete, has_mmd)
}
```

- [ ] **Step 4: Rewrite the status block in `list_patients`**

In `list_patients` (lines ~256-272), replace the loop body that computes status:

```rust
    for (id, path) in entries {
        let path_str = path.to_string_lossy().to_string();
        let (finalized_count, has_mmd) = read_annotation_summary(&app, &id);
        let has_seeds = patient_has_seeds(&app, &path_str);
        let status = if has_mmd {
            "complete"
        } else if has_seeds || finalized_count > 0 {
            "in_progress"
        } else {
            "not_started"
        };
        patients.push(PatientInfo {
            id,
            path: path_str,
            status: status.to_string(),
            finalized_count,
            has_mmd,
        });
    }
```

with:

```rust
    for (id, path) in entries {
        let path_str = path.to_string_lossy().to_string();
        let has_seeds = patient_has_seeds(&app, &path_str);
        let (complete, has_mmd) = patient_session_summary(&app, &path_str);
        let status = status_for(complete, has_seeds);
        patients.push(PatientInfo {
            id,
            path: path_str,
            status: status.to_string(),
            // Finalized contour count is not persisted in the session bundle;
            // the old annotations-based count was already always 0. Kept for API
            // stability.
            finalized_count: 0,
            has_mmd,
        });
    }
```

- [ ] **Step 5: Run tests to verify they pass + the whole lib still compiles**

Run: `cargo test --manifest-path src-tauri/Cargo.toml --lib session_flags status_for patient_file_key`
Expected: PASS (all named tests, including the pre-existing `patient_file_key` test).

- [ ] **Step 6: Commit**

```bash
git add src-tauri/src/commands/dicom.rs
git commit -m "fix: derive patient status from session bundle (complete=FAI&&MMD)"
```

---

### Task 2: Frontend — auto-save session on FAI / MMD / Water-Lipid completion

The other half of the bug fix: persist the session automatically when an analysis finishes, so the disk (and thus the list badge) reflects `complete` without a manual Save.

**Files:**
- Modify: `src/App.svelte` (import line 32; FAI run call sites ~607, ~615)
- Modify: `src/components/MmdAnalysisView.svelte` (`handleRunMmd` ~374-388)
- Modify: `src/components/WaterLipidView.svelte` (imports; `handleRun` ~111-126)

**Interfaces:**
- Consumes: `saveSession(dicomPath)` from `$lib/session`, `volumeStore.dicomPath`, `pipelineStore.run()` (async, sets `status='complete'` on success).

- [ ] **Step 1: App.svelte — import `saveSession` and wrap the FAI run**

Change the import on line 32 from:

```ts
  import { loadSession, clearSession } from '$lib/session';
```

to:

```ts
  import { loadSession, clearSession, saveSession } from '$lib/session';
```

Add this function in the `<script>` block (near the other handlers, e.g. just after `autoRestoreSession`, ~line 73):

```ts
  /** Run FAI, then auto-persist the session so the patient list shows progress
   *  without a manual Save. Footer status updates live regardless. */
  async function runFaiAndSave() {
    await pipelineStore.run();
    if (pipelineStore.status === 'complete' && volumeStore.dicomPath) {
      try {
        await saveSession(volumeStore.dicomPath);
      } catch (e) {
        console.error('auto-save after FAI failed:', e);
      }
    }
  }
```

Change the two FAI run call sites. Line ~607 (re-run button):

```svelte
          onclick={() => { pipelineStore.run(); }}
```

→

```svelte
          onclick={() => { runFaiAndSave(); }}
```

Line ~615 (Analyze button):

```svelte
          onclick={() => pipelineStore.run()}
```

→

```svelte
          onclick={() => runFaiAndSave()}
```

- [ ] **Step 2: MmdAnalysisView.svelte — auto-save after MMD**

In `handleRunMmd` (lines ~374-388), after `await refreshSurfaces();`, add the auto-save. The function becomes:

```ts
  async function handleRunMmd() {
    if (mmdBusy) return;
    mmdBusy = true;
    mmdError = '';
    overlayCache = {}; // stale now
    try {
      mmdStore.summary = await runMmdOnRoi('gls');
      await refreshSurfaces();
      if (dicomPath) {
        try {
          await saveSession(dicomPath);
        } catch (e) {
          console.error('auto-save after MMD failed:', e);
        }
      }
    } catch (err) {
      mmdError = err instanceof Error ? err.message : String(err);
      console.error('MMD failed:', err);
    } finally {
      mmdBusy = false;
    }
  }
```

(`saveSession` and `dicomPath` are already imported/derived in this file at lines 34 and 50.)

- [ ] **Step 3: WaterLipidView.svelte — add imports and auto-save after the run**

Add these imports to the `<script>` block (with the other `$lib` imports):

```ts
  import { volumeStore } from '$lib/stores/volumeStore.svelte';
  import { saveSession } from '$lib/session';
```

In `handleRun` (lines ~111-126), after `wlStore.set(await runWaterLipid());`, add:

```ts
      if (volumeStore.dicomPath) {
        try {
          await saveSession(volumeStore.dicomPath);
        } catch (e) {
          console.error('auto-save after water/lipid failed:', e);
        }
      }
```

- [ ] **Step 4: Build to verify it compiles**

Run: `npm run build`
Expected: build succeeds (no Svelte/TS errors).

- [ ] **Step 5: Manual verification**

Run `npm run tauri dev`. Load a patient, place ≥2 seeds, click **Analyze** (FAI). When it finishes, open the patient browser → the patient badge should read **In progress** (FAI only, no MMD yet). Then run MMD on a cross-section, reopen the browser → badge should read **Done** (complete). No manual Save was clicked.

- [ ] **Step 6: Commit**

```bash
git add src/App.svelte src/components/MmdAnalysisView.svelte src/components/WaterLipidView.svelte
git commit -m "fix: auto-save session on FAI/MMD/water-lipid completion"
```

---

### Task 3: Rust — patient-level flag storage, commands, and list integration

New `flags/{key}.json` sidecar keyed by the *patient* folder (parent of the series `dicomPath`), `set`/`get` commands, and `flagged`/`flag_note` folded into `list_patients`/`PatientInfo`.

**Files:**
- Create: `src-tauri/src/commands/flag.rs`
- Modify: `src-tauri/src/commands/mod.rs` (add `pub mod flag;`)
- Modify: `src-tauri/src/lib.rs` (register two commands in `generate_handler!`)
- Modify: `src-tauri/src/commands/dicom.rs` (`PatientInfo` struct + flag read in `list_patients`)

**Interfaces:**
- Produces: `PatientFlag { flagged: bool, note: String, flagged_at: String }`; commands `set_patient_flag(app, dicom_path, flagged, note, flagged_at) -> Result<(), String>` and `get_patient_flag(app, dicom_path) -> Result<Option<PatientFlag>, String>`; pure helpers `patient_dir_of(&str) -> String`, `flag_file_path(&Path, &str) -> PathBuf`, `read_flag_at(&Path) -> Result<Option<PatientFlag>, String>`.
- Consumes: `crate::commands::dicom::patient_file_key`.

- [ ] **Step 1: Create `flag.rs` with the failing tests**

Create `src-tauri/src/commands/flag.rs`:

```rust
//! Patient-level data-quality flag ("Flagged" — data has a problem / unusable).
//! The only human-authored fact in the patient-status feature, so it is the only
//! thing persisted: one tiny sidecar `flags/{patient_file_key}.json` per patient.
//!
//! `dicom_path` everywhere in the app is a *series* subfolder; flags are
//! patient-level, so we key off the parent (patient) folder — every series of a
//! patient shares one flag file. `list_patients` reads the same key.

use std::path::{Path, PathBuf};
use tauri::Manager;

use crate::commands::dicom::patient_file_key;

#[derive(serde::Serialize, serde::Deserialize, Clone)]
pub struct PatientFlag {
    pub flagged: bool,
    pub note: String,
    pub flagged_at: String,
}

/// The patient folder that owns a loaded *series* path (its parent dir).
pub(crate) fn patient_dir_of(series_path: &str) -> String {
    Path::new(series_path)
        .parent()
        .map(|p| p.to_string_lossy().to_string())
        .filter(|s| !s.is_empty())
        .unwrap_or_else(|| series_path.to_string())
}

/// `<base>/<patient_file_key(patient_dir)>` — the same key `list_patients` uses,
/// so a flag set while viewing any series is found when listing the patient.
pub(crate) fn flag_file_path(base_dir: &Path, patient_dir: &str) -> PathBuf {
    base_dir.join(patient_file_key(patient_dir))
}

/// Read a flag file, returning None if absent. Errors only on real IO/parse
/// failures (never swallows them into None).
pub(crate) fn read_flag_at(path: &Path) -> Result<Option<PatientFlag>, String> {
    match std::fs::read_to_string(path) {
        Ok(s) => {
            let f: PatientFlag =
                serde_json::from_str(&s).map_err(|e| format!("parse failed: {e}"))?;
            Ok(Some(f))
        }
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => Ok(None),
        Err(e) => Err(format!("read failed: {e}")),
    }
}

fn flags_dir(app: &tauri::AppHandle) -> PathBuf {
    let dir = app
        .path()
        .app_data_dir()
        .expect("app data dir")
        .join("flags");
    let _ = std::fs::create_dir_all(&dir);
    dir
}

#[tauri::command]
pub async fn set_patient_flag(
    app: tauri::AppHandle,
    dicom_path: String,
    flagged: bool,
    note: String,
    flagged_at: String,
) -> Result<(), String> {
    let patient_dir = patient_dir_of(&dicom_path);
    let path = flag_file_path(&flags_dir(&app), &patient_dir);
    if flagged {
        let flag = PatientFlag {
            flagged: true,
            note,
            flagged_at,
        };
        let json =
            serde_json::to_string_pretty(&flag).map_err(|e| format!("serialize failed: {e}"))?;
        std::fs::write(&path, json).map_err(|e| format!("write failed: {e}"))?;
    } else {
        match std::fs::remove_file(&path) {
            Ok(()) => {}
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => {}
            Err(e) => return Err(format!("remove failed: {e}")),
        }
    }
    Ok(())
}

#[tauri::command]
pub async fn get_patient_flag(
    app: tauri::AppHandle,
    dicom_path: String,
) -> Result<Option<PatientFlag>, String> {
    let patient_dir = patient_dir_of(&dicom_path);
    let path = flag_file_path(&flags_dir(&app), &patient_dir);
    read_flag_at(&path)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn flag_is_patient_level_two_series_one_file() {
        let base = Path::new("/base");
        let a = flag_file_path(base, &patient_dir_of("/data/P001/MonoPlus_70keV"));
        let b = flag_file_path(base, &patient_dir_of("/data/P001/CCTA"));
        assert_eq!(a, b, "two series of one patient must map to one flag file");
        assert!(a.starts_with("/base"));
    }

    #[test]
    fn flag_distinct_patients_distinct_files() {
        let base = Path::new("/base");
        let a = flag_file_path(base, &patient_dir_of("/data/P001/MonoPlus_70keV"));
        let b = flag_file_path(base, &patient_dir_of("/data/P002/MonoPlus_70keV"));
        assert_ne!(a, b);
    }

    #[test]
    fn flag_round_trip_write_read_delete() {
        let dir = std::env::temp_dir().join("pcat_flag_test_round_trip");
        let _ = std::fs::create_dir_all(&dir);
        let path = dir.join("p.json");
        let _ = std::fs::remove_file(&path);
        assert!(read_flag_at(&path).unwrap().is_none());

        let flag = PatientFlag {
            flagged: true,
            note: "vessel jump ~slice 40".into(),
            flagged_at: "2026-06-17T00:00:00Z".into(),
        };
        std::fs::write(&path, serde_json::to_string(&flag).unwrap()).unwrap();
        let got = read_flag_at(&path).unwrap().unwrap();
        assert!(got.flagged);
        assert_eq!(got.note, "vessel jump ~slice 40");

        std::fs::remove_file(&path).unwrap();
        assert!(read_flag_at(&path).unwrap().is_none());
    }
}
```

- [ ] **Step 2: Register the module**

In `src-tauri/src/commands/mod.rs`, add `pub mod flag;` (keep alphabetical):

```rust
pub mod annotation;
pub mod cpr;
pub mod dicom;
pub mod flag;
pub mod framed;
pub mod pipeline;
pub mod water_lipid;
```

- [ ] **Step 3: Run the flag tests to verify they pass**

Run: `cargo test --manifest-path src-tauri/Cargo.toml --lib flag_`
Expected: PASS (`flag_is_patient_level_two_series_one_file`, `flag_distinct_patients_distinct_files`, `flag_round_trip_write_read_delete`).

- [ ] **Step 4: Register the commands in the invoke handler**

In `src-tauri/src/lib.rs`, inside `tauri::generate_handler![ ... ]` (after the `commands::water_lipid::*` entries, before the closing `]`), add:

```rust
        commands::flag::set_patient_flag,
        commands::flag::get_patient_flag,
```

- [ ] **Step 5: Add flag fields to `PatientInfo` and read them in `list_patients`**

In `src-tauri/src/commands/dicom.rs`, extend the `PatientInfo` struct (lines ~140-152) by adding two fields before the closing brace:

```rust
    /// Whether MMD has been run and stored in saved annotations.
    pub has_mmd: bool,
    /// Whether the patient's data has been flagged as problematic / unusable.
    pub flagged: bool,
    /// Optional note explaining the flag (None when not flagged).
    pub flag_note: Option<String>,
}
```

Add this import near the top of `dicom.rs` (with the other `use crate::...` or `use super::...` lines):

```rust
use crate::commands::flag::{flag_file_path, read_flag_at};
```

In `list_patients`, compute a `flags` base dir once before the `for (id, path) in entries` loop:

```rust
    let flags_base = app
        .path()
        .app_data_dir()
        .expect("app data dir")
        .join("flags");
```

Then inside the loop, after computing `status` and before `patients.push(...)`, read the flag and include it:

```rust
        let flag = read_flag_at(&flag_file_path(&flags_base, &path_str))
            .ok()
            .flatten();
        patients.push(PatientInfo {
            id,
            path: path_str,
            status: status.to_string(),
            finalized_count: 0,
            has_mmd,
            flagged: flag.as_ref().map(|f| f.flagged).unwrap_or(false),
            flag_note: flag.and_then(|f| if f.note.is_empty() { None } else { Some(f.note) }),
        });
```

(Replace the `patients.push(PatientInfo { ... })` you wrote in Task 1 Step 4 with this expanded version.)

- [ ] **Step 6: Verify the whole lib compiles and all Rust tests pass**

Run: `cargo test --manifest-path src-tauri/Cargo.toml --lib`
Expected: PASS — all tests, no warnings about unused `read_annotation_summary` (it was removed in Task 1).

- [ ] **Step 7: Commit**

```bash
git add src-tauri/src/commands/flag.rs src-tauri/src/commands/mod.rs src-tauri/src/lib.rs src-tauri/src/commands/dicom.rs
git commit -m "feat: patient-level data-quality flag (storage, commands, list integration)"
```

---

### Task 4: Frontend — flag API wrappers, `flagStore`, and `PatientInfo` type

**Files:**
- Modify: `src/lib/api.ts` (add `PatientFlag` type + 2 wrappers; extend `PatientInfo`)
- Create: `src/lib/stores/flagStore.svelte.ts`

**Interfaces:**
- Produces: `PatientFlag` type `{ flagged: boolean; note: string; flagged_at: string }`; `getPatientFlag(dicomPath)`, `setPatientFlag(dicomPath, flagged, note, flaggedAt)`; `flagStore` with `current`, `load(dicomPath)`, `set(dicomPath, flagged, note)`, `clear()`.
- Consumes: `invoke` (existing pattern: `invoke<Ret>('command', { camelCaseArgs })`).

- [ ] **Step 1: Add the `PatientFlag` type and wrappers to `api.ts`**

In `src/lib/api.ts`, add near the `PatientInfo` type / session wrappers:

```ts
/** Patient-level data-quality flag ("Flagged" — data has a problem / unusable). */
export type PatientFlag = { flagged: boolean; note: string; flagged_at: string };

/** Read a patient's flag (the patient folder is derived backend-side from the
 *  series `dicomPath`). Returns null if not flagged. */
export async function getPatientFlag(dicomPath: string): Promise<PatientFlag | null> {
  return invoke<PatientFlag | null>('get_patient_flag', { dicomPath });
}

/** Set or clear a patient's flag. `flagged=false` deletes the flag file. */
export async function setPatientFlag(
  dicomPath: string,
  flagged: boolean,
  note: string,
  flaggedAt: string,
): Promise<void> {
  await invoke('set_patient_flag', { dicomPath, flagged, note, flaggedAt });
}
```

- [ ] **Step 2: Extend the `PatientInfo` type**

In `src/lib/api.ts`, add two fields to the `PatientInfo` type (after `has_mmd`):

```ts
  /** Whether MMD has been run and stored in saved annotations. */
  has_mmd: boolean;
  /** Whether the patient's data is flagged as problematic / unusable. */
  flagged: boolean;
  /** Optional note explaining the flag (null when not flagged). */
  flag_note: string | null;
};
```

- [ ] **Step 3: Create `flagStore.svelte.ts`**

Create `src/lib/stores/flagStore.svelte.ts` (modeled on `wlStore.svelte.ts`):

```ts
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
  /** Clear on patient switch. */
  clear(): void {
    current = null;
  },
};
```

- [ ] **Step 4: Build to verify it compiles**

Run: `npm run build`
Expected: build succeeds.

- [ ] **Step 5: Commit**

```bash
git add src/lib/api.ts src/lib/stores/flagStore.svelte.ts
git commit -m "feat: flag API wrappers + flagStore"
```

---

### Task 5: Frontend — footer (patient id, live status, flag indicator + dialog)

**Files:**
- Create: `src/lib/patientStatus.ts` (pure `derivePatientStatus` + `patientIdOf`)
- Modify: `src/App.svelte` (imports; `autoRestoreSession` to load the flag; footer block ~764-796; flag-dialog state + markup)

**Interfaces:**
- Produces: `derivePatientStatus({ hasSeeds, hasFai, hasMmd }) -> 'not_started'|'in_progress'|'complete'`; `patientIdOf(dicomPath) -> string`.
- Consumes: `flagStore`, `mmdStore`, `pipelineStore.results`, `seedStore.vessels`, `volumeStore`.

- [ ] **Step 1: Create the pure helper module**

Create `src/lib/patientStatus.ts`:

```ts
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
```

- [ ] **Step 2: Add imports and flag-load to App.svelte**

Add to the `<script>` imports in `src/App.svelte`:

```ts
  import { mmdStore } from '$lib/stores/mmdStore.svelte';
  import { flagStore } from '$lib/stores/flagStore.svelte';
  import { derivePatientStatus, patientIdOf } from '$lib/patientStatus';
```

In `autoRestoreSession` (lines ~63-72), load the flag alongside the session. After the `try { ... } catch { ... }` block that restores the session, add:

```ts
    await flagStore.load(dicomPath);
```

(so the full tail of `autoRestoreSession` is the existing session restore followed by `await flagStore.load(dicomPath);`).

Add flag-dialog state near the other `$state` declarations (e.g. after `let errorMessage = $state('');`, ~line 49):

```ts
  let showFlagDialog = $state(false);
  let flagNoteDraft = $state('');
```

- [ ] **Step 3: Replace the footer "Volume loaded" branch**

In `src/App.svelte`, replace this branch (lines ~787-789):

```svelte
      {:else if volumeStore.current}
        <span class="h-1.5 w-1.5 rounded-full bg-success"></span>
        <span class="text-[11px] text-text-secondary">Volume loaded</span>
```

with:

```svelte
      {:else if volumeStore.current}
        {@const st = derivePatientStatus({
          hasSeeds: (['LAD', 'LCx', 'RCA'] as const).some(
            (v) => seedStore.vessels[v].seeds.length > 0,
          ),
          hasFai: pipelineStore.results !== null,
          hasMmd: mmdStore.summary !== null,
        })}
        <span class="h-1.5 w-1.5 rounded-full {st === 'complete' ? 'bg-success' : 'bg-warning'}"></span>
        <span class="text-[11px] font-medium text-text-primary">{patientIdOf(volumeStore.dicomPath)}</span>
        {#if volumeStore.current.studyDescription}
          <span class="truncate text-[11px] text-text-secondary">· {volumeStore.current.studyDescription}</span>
        {/if}
        <span
          class="rounded px-1.5 text-[10px] font-medium {st === 'complete'
            ? 'bg-success/15 text-success'
            : 'bg-warning/15 text-warning'}"
        >
          {st === 'complete' ? 'complete' : 'in progress'}
        </span>
        {#if flagStore.current?.flagged}
          <button
            class="rounded bg-error/15 px-1.5 text-[10px] font-medium text-error hover:bg-error/25"
            title={flagStore.current.note || 'Flagged: data problem / unusable'}
            onclick={() => { flagNoteDraft = flagStore.current?.note ?? ''; showFlagDialog = true; }}
          >
            ⚠ Flagged
          </button>
        {:else}
          <button
            class="rounded px-1.5 text-[10px] text-text-secondary hover:bg-surface-tertiary hover:text-text-primary"
            title="Flag this data as problematic / unusable"
            onclick={() => { flagNoteDraft = ''; showFlagDialog = true; }}
          >
            ⚑ Flag
          </button>
        {/if}
```

- [ ] **Step 4: Add the flag dialog**

In `src/App.svelte`, just after the patient-browser modal block (after line ~761, before the `<footer>`), add:

```svelte
  <!-- ===== Flag dialog ===== -->
  {#if showFlagDialog}
    <div
      class="absolute inset-0 z-50 flex items-center justify-center bg-black/40"
      onclick={(e) => { if (e.target === e.currentTarget) showFlagDialog = false; }}
    >
      <div class="w-80 rounded-lg border border-border bg-surface-secondary p-4 shadow-xl">
        <h3 class="mb-1 text-sm font-semibold text-text-primary">Flag this data</h3>
        <p class="mb-2 text-[11px] text-text-secondary">
          Mark this patient's data as problematic / unusable — e.g. vessel
          discontinuity, motion or recon artifact.
        </p>
        <textarea
          bind:value={flagNoteDraft}
          rows="3"
          placeholder="Optional note (what's wrong)…"
          class="w-full rounded border border-border bg-surface px-2 py-1 text-xs text-text-primary focus:border-accent focus:outline-none"
        ></textarea>
        <div class="mt-3 flex items-center justify-between gap-2">
          {#if flagStore.current?.flagged}
            <button
              class="rounded px-2 py-1 text-[11px] text-text-secondary hover:bg-surface-tertiary"
              onclick={async () => {
                if (volumeStore.dicomPath) await flagStore.set(volumeStore.dicomPath, false, '');
                showFlagDialog = false;
              }}
            >
              Remove flag
            </button>
          {:else}
            <span></span>
          {/if}
          <div class="flex gap-2">
            <button
              class="rounded px-2 py-1 text-[11px] text-text-secondary hover:bg-surface-tertiary"
              onclick={() => (showFlagDialog = false)}
            >
              Cancel
            </button>
            <button
              class="rounded bg-error/15 px-3 py-1 text-[11px] font-medium text-error hover:bg-error/25"
              onclick={async () => {
                if (volumeStore.dicomPath)
                  await flagStore.set(volumeStore.dicomPath, true, flagNoteDraft.trim());
                showFlagDialog = false;
              }}
            >
              Flag
            </button>
          </div>
        </div>
      </div>
    </div>
  {/if}
```

- [ ] **Step 5: Build to verify it compiles**

Run: `npm run build`
Expected: build succeeds.

- [ ] **Step 6: Manual verification**

`npm run tauri dev`. Load a patient → footer shows the patient id (folder name) + study description + a status chip (`in progress` / `complete`). Click **⚑ Flag**, type a note, click **Flag** → footer shows red **⚠ Flagged**. Reload the app (or switch away and back) → the flag persists. Open the dialog again → **Remove flag** clears it.

- [ ] **Step 7: Commit**

```bash
git add src/lib/patientStatus.ts src/App.svelte
git commit -m "feat: footer shows patient id, live status, and flag control"
```

---

### Task 6: Frontend — highlight loaded patient + show flag in PatientBrowser

**Files:**
- Modify: `src/App.svelte` (pass `currentPath` prop to `<PatientBrowser>` ~752-760)
- Modify: `src/components/PatientBrowser.svelte` (props; row highlight; flag badge; optional `flagged` filter)

**Interfaces:**
- Consumes: `volumeStore.dicomPath` (series path), `PatientInfo.flagged`, `PatientInfo.flag_note`.

- [ ] **Step 1: Pass the loaded path into PatientBrowser**

In `src/App.svelte`, add the `currentPath` prop to the `<PatientBrowser ... />` element (lines ~752-760):

```svelte
    <PatientBrowser
      currentPath={volumeStore.dicomPath}
      onSelect={(path) => { showPatientBrowser = false; loadFromPath(path); }}
```

(keep the other props unchanged.)

- [ ] **Step 2: Accept the prop in PatientBrowser**

In `src/components/PatientBrowser.svelte`, add `currentPath` to the `Props` type and the `$props()` destructuring (lines ~17-39):

```ts
  type Props = {
    /** Initial root directory (the user can edit it before scanning). */
    initialRootDir?: string;
    /** The currently-loaded series path, to highlight its patient. */
    currentPath?: string | null;
    /** Called for a regular single-series pick. */
    onSelect: (path: string) => void;
```

and:

```ts
  let {
    initialRootDir = '/Volumes/Molloilab/Shu Nie/UCI NAEOTOM CCTA Data',
    currentPath = null,
    onSelect,
    onSelectDualEnergy,
    onSelectPatient,
    onClose,
  }: Props = $props();
```

- [ ] **Step 3: Highlight the loaded patient's row + show the flag**

In `src/components/PatientBrowser.svelte`, in the `{#each filtered as p (p.id)}` block (lines ~246-291), change the row container. Replace:

```svelte
        {#each filtered as p (p.id)}
          {@const isOpen = !!expanded[p.id]}
          <div class="border-b border-border/60">
            <!-- Patient row: click expands to show series subfolders -->
            <div
              class="flex w-full items-center justify-between gap-2 px-4 py-2 hover:bg-accent/5"
            >
```

with (adds `isLoaded` and a left-accent/tint highlight):

```svelte
        {#each filtered as p (p.id)}
          {@const isOpen = !!expanded[p.id]}
          {@const isLoaded = !!currentPath && currentPath.startsWith(p.path)}
          <div class="border-b border-border/60">
            <!-- Patient row: click expands to show series subfolders -->
            <div
              class="flex w-full items-center justify-between gap-2 px-4 py-2 hover:bg-accent/5 {isLoaded
                ? 'border-l-2 border-accent bg-accent/10'
                : ''}"
            >
```

Then add a flag badge in the status area. Replace the status badge block (lines ~282-290):

```svelte
                <span
                  class="rounded px-2 py-0.5 text-[10px] font-medium {p.status === 'complete'
                    ? 'bg-success/15 text-success'
                    : p.status === 'in_progress'
                      ? 'bg-warning/15 text-warning'
                      : 'bg-text-secondary/15 text-text-secondary'}"
                >
                  {statusLabel(p.status)}
                </span>
```

with (adds the red flag marker before the status badge):

```svelte
                {#if p.flagged}
                  <span
                    class="rounded bg-error/15 px-1.5 py-0.5 text-[10px] font-medium text-error"
                    title={p.flag_note || 'Flagged: data problem / unusable'}
                  >
                    ⚠ Flagged
                  </span>
                {/if}
                <span
                  class="rounded px-2 py-0.5 text-[10px] font-medium {p.status === 'complete'
                    ? 'bg-success/15 text-success'
                    : p.status === 'in_progress'
                      ? 'bg-warning/15 text-warning'
                      : 'bg-text-secondary/15 text-text-secondary'}"
                >
                  {statusLabel(p.status)}
                </span>
```

- [ ] **Step 4: Add a "Flagged" filter option**

In `src/components/PatientBrowser.svelte`, widen `statusFilter` (line ~52):

```ts
  let statusFilter = $state<'all' | 'not_started' | 'in_progress' | 'complete' | 'flagged'>('all');
```

Add the option to the `<select>` (after the "Done" option, ~line 227):

```svelte
        <option value="complete">Done ({counts.cp})</option>
        <option value="flagged">Flagged ({counts.fl})</option>
```

Find the `filtered` derivation (the `$derived` that filters by `query` + `statusFilter`) and the `counts` derivation. In the filter predicate, handle `flagged` specially (it filters on `p.flagged`, not `p.status`). The status-match clause currently looks like `statusFilter === 'all' || p.status === statusFilter`; change it to:

```ts
      const statusOk =
        statusFilter === 'all'
          ? true
          : statusFilter === 'flagged'
            ? p.flagged
            : p.status === statusFilter;
```

and use `statusOk` in the existing `filtered` predicate (combined with the existing text-query check). In the `counts` derivation, add `fl`:

```ts
    fl: patients.filter((p) => p.flagged).length,
```

(`patients` is the unfiltered list the `counts` derivation already iterates — match its existing variable name.)

- [ ] **Step 5: Build to verify it compiles**

Run: `npm run build`
Expected: build succeeds.

- [ ] **Step 6: Manual verification**

`npm run tauri dev`. Load a patient, open the patient browser → that patient's row is highlighted (accent tint + left bar). Flag the patient from the footer, reopen the browser → the row shows **⚠ Flagged** next to its status badge, and selecting **Flagged** in the filter dropdown lists only flagged patients.

- [ ] **Step 7: Commit**

```bash
git add src/App.svelte src/components/PatientBrowser.svelte
git commit -m "feat: highlight loaded patient + show flag in patient browser"
```

---

## Self-Review

**Spec coverage:**
- Footer shows which patient + live status → Task 5. ✓
- Highlight loaded patient in list → Task 6. ✓
- Patient-level boolean flag + optional note, set from footer → Tasks 3 (backend), 4 (store), 5 (dialog). ✓
- Flag displayed in list → Task 6. ✓
- Auto-`complete` bug root cause fix (derive from session bundle + auto-save) → Tasks 1 + 2. ✓
- Term `Flagged` → Tasks 5, 6 markup. ✓
- Non-goals (categories, series-level, manual complete, flag-from-list) → none implemented. ✓
- Verification items: dicomPath = series path resolved (Task 5 `patientIdOf`, Task 3 `patient_dir_of`, Task 6 prefix highlight); legacy annotations-only patients read as not_started/in_progress (acceptable — documented in spec §11); FAI-before-MMD ordering is handled correctly by the per-file `f && m` check regardless of order.

**Placeholder scan:** No TBD/TODO; every code step shows complete code; test code is concrete.

**Type consistency:** `PatientFlag { flagged, note, flagged_at }` identical in Rust (Task 3) and TS (Task 4). `set_patient_flag(dicom_path, flagged, note, flagged_at)` ↔ `setPatientFlag(dicomPath, flagged, note, flaggedAt)` (Tauri camelCase mapping). `PatientInfo.flagged/flag_note` added in both Rust (Task 3) and TS (Task 4). `derivePatientStatus` signature (Task 5) matches the Rust `status_for`/`session_flags` rule (Task 1). `flagStore.set/load/clear/current` consistent between Task 4 (definition) and Task 5 (use).

**Known coupling (intentional, documented in code):** `session_flags` (Rust) probes the field names `fai`, `mmd.summary` written by `saveSession` (session.ts) — pinned by the Task 1 unit tests. `derivePatientStatus` (TS) mirrors `status_for` (Rust) — both trivially simple, Rust side is the tested canonical one.
