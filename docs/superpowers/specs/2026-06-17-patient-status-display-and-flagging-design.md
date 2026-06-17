# Patient identity, workflow status, and data-quality flagging — Design

- **Date:** 2026-06-17
- **Status:** Draft (awaiting review)
- **Scope:** Single implementation plan (one feature, plus one root-cause bug fix it depends on)

## 1. Motivation

While reviewing a patient (e.g. an `r`-prefixed case), the operator noticed the
coronary vessel jump from one position to another between slices — a likely
acquisition problem (motion or reconstruction artifact) that makes the vessel
discontinuous and the dataset unreliable. Two gaps surfaced:

1. **You cannot tell which patient is currently loaded.** The footer shows only
   the generic text `Volume loaded`. The patient list does not highlight the
   loaded case.
2. **There is no way to mark a dataset as problematic.** The operator wants to
   record "this data has a problem / is unusable" so it is visible later.

A third, related problem was discovered during investigation:

3. **The `complete` status never appears automatically**, even after FAI and MMD
   have been run. This turned out to be a pre-existing bug (Section 3).

## 2. Goals and non-goals

### Goals
- Show **which patient** is loaded in the footer, plus its live workflow status.
- **Highlight** the currently-loaded patient in the patient list.
- Add a **per-patient `Flagged` mark** ("data has a problem / unusable") with an
  optional free-text note, settable from the footer while viewing the patient.
- **Fix** the auto-`complete` bug so status reflects analysis progress with no
  manual save.

### Non-goals (YAGNI — explicitly deferred)
- **No flag categories.** A fixed dropdown (motion / recon / …) has too many
  possibilities to enumerate; the optional free-text note covers "why" when
  wanted.
- **No series/volume-level flag.** Flags are **patient-level**, matching the
  existing storage grain of seeds and annotations. (Motion affects the whole
  acquisition anyway; recon-specific cases can be described in the note.)
- **No manual "mark complete" button.** Completion is fully auto-derived.
- **No flagging from the patient list.** v1 sets the flag from the footer (the
  currently-viewed patient); the list only *displays* it.

## 3. Root cause of the "`complete` never auto-shows" bug

The per-patient status in the list is **derived** on the Rust side by
`list_patients`, which reads sidecar files from disk. Tracing the write path vs.
the read path revealed a producer/consumer mismatch — a Single-Source-of-Truth
violation where one fact ("has MMD run?") is stored in two places and the writer
and reader disagree on which.

**Read path** (`src-tauri/src/commands/dicom.rs`):
- `list_patients` → `read_annotation_summary(app, id)` reads
  `app_data_dir/annotations/{key}.json` and sets `has_mmd` from the presence of
  the `mmd_method` field.
- Status rule: `has_mmd → "complete"`, else `has_seeds || finalized_count > 0 →
  "in_progress"`, else `"not_started"`.

**Write path** (frontend):
- The only save the UI performs is `saveSession(dicomPath)`
  (`src/lib/session.ts`), which calls the `save_session` command and writes
  `app_data_dir/sessions/{key}.json` — an **opaque bundle** `{ seeds, fai, mmd,
  wl }`. It **never** writes `annotations/{key}.json` and **never** writes a
  `mmd_method` field.
- `save_annotations` (the function that *does* write `mmd_method`) is **dead
  code** — no caller in the frontend.
- Worse, **no save happens automatically** after FAI/MMD complete; results live
  only in `AppState` (Rust) and the in-memory Svelte stores until the user
  clicks an explicit "Save".

**Consequence:** `annotations/{key}.json` never receives `mmd_method`, so
`has_mmd` is always `false`, so the badge never reaches `complete`.

### Fix direction (root cause, not symptom)
Consolidate to **one source of truth**: the session bundle the frontend actually
writes.

1. **Derive status from the session bundle.** `list_patients` reads
   `sessions/{key}.json` and derives:
   - `has_fai`  = bundle `fai` is a non-null object with ≥1 vessel entry
   - `has_mmd`  = bundle `mmd.summary` is non-null
   - `has_seeds`= bundle `seeds` has any vessel with ≥1 seed
   - status: `complete` iff `has_fai && has_mmd`; else `in_progress` iff
     `has_seeds`; else `not_started`.
   The `annotations/{key}.json` / `mmd_method` path is no longer consulted for
   status. We do **not** "fix" this by writing both files in sync — that is the
   anti-pattern that caused the bug.
2. **Auto-save on completion.** After FAI completes, after MMD completes, and
   after whole-volume Water/Lipid completes, call `saveSession(dicomPath)` so the
   disk reflects progress without a manual click. This also prevents silent loss
   of analysis work.

> Note on coupling: `list_patients` now probes three field names (`fai`,
> `mmd.summary`, `seeds`) inside the frontend-owned bundle. This is a small,
> deliberate contract between `session.ts` and `list_patients`, guarded by a test
> fixture (Section 8). It is preferable to storing a derived `status` field in
> the bundle (which would violate "derive, don't sync").

> Verification item for the plan: confirm whether the UI requires FAI before MMD.
> If it does, `complete` transitively requires both already; if not, a
> MMD-without-FAI case correctly stays `in_progress` under the rule above. Either
> way the derivation is correct.

## 4. Design principles

This feature is almost entirely **display + derivation**. The only genuinely new
*stored* fact is the human judgement "this data is problematic" — which cannot be
derived and therefore must be persisted. Everything else (the workflow status) is
derived from existing data, per the project's SSoT rules ("a fact lives in
exactly one place; derive, don't sync").

## 5. Data model — the one new stored fact

A patient-level flag, stored in a dedicated sidecar parallel to the existing
`seeds/`, `annotations/`, `sessions/` directories:

```
app_data_dir/flags/{patient_file_key}.json
```

```jsonc
// PatientFlag
{
  "flagged": true,
  "note": "proximal RCA jumps position ~slice 40 — suspect motion/recon",
  "flagged_at": "2026-06-17T21:30:00.000Z"  // ISO 8601, stamped by the frontend at flag time
}
```

- Keyed by `patient_file_key(dicom_path)` — the existing per-patient key
  (Section 3), so the grain matches the chosen patient-level scope.
- **Independent of analysis state.** A patient can be flagged before any seeds
  exist (the operator opens it, sees the discontinuity, flags it immediately).
  This is why the flag is *not* folded into `annotations/` or the session bundle.
- **Un-flagging deletes the file** (`flagged=false` ⇒ remove). Presence of the
  file with `flagged=true` is the flag. No stale residue, no `false` rows to
  reason about. (If flag history is wanted later, add it then.)

## 6. Backend changes (Rust)

New module `src-tauri/src/commands/flag.rs` (registered in `commands/mod.rs` and
the invoke handler in `lib.rs`):

```rust
pub struct PatientFlag {
    pub flagged: bool,
    pub note: String,
    pub flagged_at: String, // ISO 8601, supplied by the frontend
}

#[tauri::command]
pub async fn set_patient_flag(
    app: AppHandle,
    dicom_path: String,
    flagged: bool,
    note: String,
    flagged_at: String,
) -> Result<(), String>;     // flagged=true → write; false → delete file

#[tauri::command]
pub async fn get_patient_flag(
    app: AppHandle,
    dicom_path: String,
) -> Result<Option<PatientFlag>, String>;

// helper, used by list_patients (reads flags/{key}.json, returns None if absent)
fn read_patient_flag(app: &AppHandle, key: &str) -> Option<PatientFlag>;
```

`src-tauri/src/commands/dicom.rs`:
- Extend `PatientInfo` with `flagged: bool` and `flag_note: Option<String>`.
- Rewrite the status derivation in `list_patients` to read `sessions/{key}.json`
  (Section 3) and fold in `read_patient_flag` for the flag fields.
- Writes the flag timestamp as supplied by the caller (no new time dependency in
  Rust).

The flag is **orthogonal** to the progress status: `PatientInfo.status` stays
`not_started | in_progress | complete`; `flagged` is a separate field. A
`complete` patient can also be `flagged` (analysis done, but a data-quality
caveat is recorded).

## 7. Frontend changes (Svelte)

### 7.1 New store `src/lib/stores/flagStore.svelte.ts`
Owns only the **currently-loaded** patient's flag (focused unit; patient identity
is not duplicated — it is read from `volumeStore`).

```ts
let current = $state<PatientFlag | null>(null);
async function load(dicomPath: string): Promise<void>; // get_patient_flag
async function set(dicomPath, flagged, note): Promise<void>; // set_patient_flag; stamps flaggedAt
function clear(): void;  // on patient switch
```

`src/lib/api.ts`: add `PatientFlag` type, `getPatientFlag`, `setPatientFlag`
wrappers, and the two new `PatientInfo` fields. `flagStore.load(dicomPath)` is
invoked wherever the session is restored on patient load (alongside
`autoRestoreSession`).

### 7.2 Pure status helper (the deliberate rule mirror)
`derivePatientStatus({ hasSeeds, hasFai, hasMmd })` →
`'not_started' | 'in_progress' | 'complete'`, used by the footer with live store
values (`seedStore` has seeds; `pipelineStore.status === 'complete'`;
`mmdStore.summary !== null`).

> This rule exists in two runtimes by necessity: Rust `list_patients` (disk scan,
> for non-loaded patients) and this TS helper (in-memory, live, for the loaded
> patient). They are the persisted vs. live views of the same rule — not two
> sources for the same value at one time. The rule is kept trivially simple and
> is tested on both sides (Section 8) to prevent drift.

### 7.3 Footer (`src/App.svelte`, current `Volume loaded` branch)
When a volume is loaded, replace `Volume loaded` with:
- **Patient identifier** matching the list's `id` (the patient folder name; see
  verification note below) + study description (secondary, truncated).
- **Live status chip**: `in progress` / `complete`, from `derivePatientStatus`.
  Turns green the instant MMD finishes — no disk round-trip.
- **Flag area**: a `⚑` button that opens the flag dialog. When
  `flagStore.current?.flagged`, render a red `⚠ Flagged` with the note on hover.

During loading, include the patient identifier in the existing progress message
(currently only the series name is shown).

> Verification item for the plan: confirm whether `volumeStore.current.dicomPath`
> is the patient folder or a series subfolder. Show the patient folder name (same
> `id` as the list). If `dicomPath` is a series subpath, derive the patient id
> from its parent / the entry under the patients root.

### 7.4 Flag dialog
Compact popover/modal anchored to the footer `⚑`:
- A toggle: **Flagged** — "data has a problem / unusable".
- An **optional** one-line note (free text; leave blank to record just the flag).
- Actions: **Save** / **Remove flag**.
- On save → `flagStore.set(...)` → footer and list reflect it.

### 7.5 Patient list (`src/components/PatientBrowser.svelte`)
- Accept a `currentPath` prop (= `volumeStore.dicomPath`). **Highlight** the
  patient row whose `path` is a prefix of `currentPath` (prefix match handles a
  loaded series subfolder): background tint + left accent border.
- For `flagged` patients, render a red `⚠` next to the existing status badge
  (`title` = note). The flag is additive — it does not replace the
  `not_started/in_progress/complete` badge.
- (Optional, low-priority) Add `flagged` to the existing status filter so the
  operator can list all problematic datasets.

### 7.6 Auto-save wiring (the bug fix, frontend half)
Add `await saveSession(dicomPath)` after the run resolves successfully (in the
handler, not as a reactive effect, to avoid double-firing) in:
- the FAI run handler (after `run_pipeline` resolves and results are stored),
- `MmdAnalysisView.handleRunMmd` (after `mmdStore.summary` is set + surfaces
  refreshed),
- `WaterLipidView.handleRun` (after `wlStore.set(...)`).
Guarded so a save failure surfaces (no silent swallow, per the project's
fail-loud rule).

## 8. Terminology

UI label: **`Flagged`**. Considered and rejected: `Non-diagnostic` (accurate but
too clinical/jargon and visually heavy in a 10px badge), `Excluded`, `Suspect`,
`Limited`. `Flagged`'s usual weakness (vague valence — "flagged for what?") is
neutralized here: the flag has exactly one meaning ("data problem / unusable"),
reinforced by red colour + `⚠` icon + the optional note.

## 9. Status lifecycle (summary)

```
not_started ──(seeds placed)──▶ in_progress ──(FAI && MMD done)──▶ complete
     │                               │                                │
     └───────────── orthogonal ──────┴──────── Flagged (red ⚠) ───────┘
                          (set/cleared manually, any time)
```

- Progress axis: derived (live in footer; from disk in list).
- Flag axis: stored; the only manual input in the whole feature.

## 10. Testing strategy (TDD)

- **Regression test for the bug (write first, must fail before the fix):** given
  a `sessions/{key}.json` fixture containing non-empty `fai` and non-null
  `mmd.summary`, `list_patients` returns `status == "complete"`. Variants:
  fai-only → `in_progress`; seeds-only → `in_progress`; empty → `not_started`.
  This pins the producer/consumer contract (Section 3 coupling note).
- **Flag round-trip:** `set_patient_flag(flagged=true, note)` then
  `get_patient_flag` returns it; `list_patients` reports `flagged=true` +
  `flag_note`; `set_patient_flag(flagged=false)` deletes the file and
  `get_patient_flag` → `None`.
- **Status rule invariant (both sides):** a small table-driven test of
  `derivePatientStatus` (TS) and the Rust derivation over the same
  `(hasSeeds, hasFai, hasMmd)` inputs, asserting identical outputs — the
  anti-drift guard for the mirrored rule.
- **Manual verification:** load a patient → footer shows its id + status; run FAI
  then MMD → footer flips to `complete` and the list badge shows `complete` on
  reopen (no manual save); flag it → red `⚠` in footer and list; reload the app →
  flag persists.

## 11. Open / verification items (resolve during planning, not now)
1. Is FAI a prerequisite for MMD in the UI? (Affects nothing in the rule, but
   documents whether `complete` is reachable via MMD alone.)
2. Exact meaning of `volumeStore.current.dicomPath` (patient folder vs. series
   subfolder) for the footer patient id and the list prefix-highlight.
3. Any legacy patients whose status lives only in old `annotations/{key}.json`
   (no session bundle). They will read as `in_progress`/`not_started` until
   re-saved. Given the bug means nothing currently shows `complete`, no working
   data is lost; flag if a migration is wanted.
