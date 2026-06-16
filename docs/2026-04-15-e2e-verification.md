# End-to-End Verification: PVAT MMD Pipeline

**Date:** 2026-04-15
**Reference patient:** `57955439` at `/Volumes/Molloilab/Shu Nie/UCI NAEOTOM CCTA Data/57955439/`
**Available series:** `MonoPlus_70keV`, `MonoPlus_150keV` (dual-energy pair), `CCTA_Bv44__CORONARY  Bv44 Q2 75%` (single-energy reference), `CA_SCORING`.

## Purpose

Walk the full pipeline end-to-end on patient 57955439 to confirm Phases 0–5 fit
together, surface any UI/data path defects, and establish a baseline for the
remaining 65 patients.

## Pre-flight

- [ ] Latest `feature/pvat-mmd-pipeline` branch checked out.
- [ ] `cargo test --workspace` exit 0 (97 tests pass + 1 integration).
- [ ] `npm run build` exit 0.
- [ ] SMB share `/Volumes/Molloilab/Shu Nie/UCI NAEOTOM CCTA Data` mounted and readable.
- [ ] Run `npm run tauri dev` — app boots, header shows "Patients" and "Open DICOM" buttons.

## Scenario 1 — Single-energy editor flow

Goal: confirm seeds + centerline + CPR still work after Phase 5 changes.

1. Click **Patients** → modal opens, lists 28 folders, status `not_started`.
2. Filter by `5795`, click `57955439` → modal closes, volume loads.
3. Verify MPR panes render (axial + sagittal + coronal).
4. Click 3 seeds along the proximal RCA in the axial view.
5. Switch active vessel `1` → RCA, `2` → LAD, etc., with seeds preserved.
6. Press `Cmd+Z` once → last seed undone.
7. Re-open the patient via **Patients** modal → seeds auto-restore.

**Pass criteria:** MPR responsive, seeds persist, undo works, vessel switching
works, no console errors.

## Scenario 2 — CPR + cross-section

1. With ≥3 RCA seeds placed, the CPR view (Editor tab) should render the
   straightened RCA strip.
2. Drag the rotation slider — image rotates.
3. Set needle position A/B/C — cross-section renders with vessel diameter
   measurement.

**Pass criteria:** CPR strip is contiguous (no black bands), cross-sections
show vessel lumen as bright disk surrounded by dark fat ring.

## Scenario 3 — MMD analysis flow (Phase 2–5 path)

1. Switch to **MMD Analysis** tab.
2. Series selector should appear (or auto-pick the dual-energy pair).
   Confirm 70 keV + 150 keV are correctly identified (check `ImageComments`
   override — the `SeriesDescription` may be mislabeled).
3. **Generate annotation targets** — strip populates with 20 thumbnails along
   proximal 40 mm of RCA. Each shows vessel wall (red dashed) and 1-vessel-
   diameter init circle (blue dashed).
4. Pick target #0 — large view shows vessel cross-section. Initialize snake
   from the blue circle. Click **Evolve** → snake contracts toward fat–tissue
   boundary.
5. If the snake misses a triangular fat extension: click **Add Point**, drop
   a waypoint near the apex, drag to pull, **Evolve** again.
6. **Accept** → status badge turns green for that target.
7. Repeat for at least 5 targets to enable mass MMD.
8. Click **Run MMD** → progress indicator updates per iteration. Should converge
   in 30–80 iterations.
9. After MMD: overlay selector enabled. Switch material (water / lipid / iodine
   / calcium) and unit (fraction / mass) — large view recolors.
10. Right panel: 3D surface plot for the selected material renders. Slide arc-
    length slider — surface morphs across the 20 cross-sections.

**Pass criteria:**
- Lumen voxels classify as high iodine + high water, very low lipid.
- Fat ring voxels classify as high lipid (>40%), low water, ~zero iodine.
- Calcium voxels (e.g. CA scoring region if visible) — high calcium.
- **Radial symmetry expected in circular fat zones**: at a given (target_idx,
  theta), lipid_frac should be roughly flat for r > 0 mm until r approaches
  contour boundary.
- **Angular asymmetry expected in triangular extensions**: at the apex theta,
  surface should have a longer radial reach (larger r_max).

## Scenario 4 — Save / load / export

1. Click **Save** in MMD view footer → "Saved" badge appears for ~2s.
2. Quit app, restart, **Patients** → patient should now show `complete` badge
   (MMD ran + ≥1 contour finalized).
3. Re-open the patient → snake_points and finalized status restored.
4. Click **Export CSV** → file dialog → save as `57955439_mmd_surfaces.csv`.
5. Open CSV — header line should be:
   `patient_id,target_index,arc_mm,theta_deg,r_mm,lipid_frac,lipid_mass,water_frac,water_mass,iodine_frac,iodine_mass,calcium_frac,calcium_mass,total_density`
6. Spot-check rows: lipid_frac + water_frac + iodine_frac + calcium_frac ≈ 1.0
   (volume conservation).

**Pass criteria:** CSV parseable as proper CSV; row counts ≈ 20 cross-sections
× 16 thetas × ~40 radial steps minus NaN-filtered points (well-defined
contour boundary).

## Scenario 5 — Patient browser cross-cohort

1. Open **Patients** modal. Filter by `complete` status.
2. Patient 57955439 should appear with `… contours · MMD` annotation count.
3. Filter by `not_started` — should list the 27 untouched patients.
4. Filter by ID `512` — should narrow to ~5 patients with `512` prefix.

**Pass criteria:** filters work, badges accurate.

## Defect log

Record observed issues here as you go. Format:

```
[severity] scenario.step — what happened — expected — repro steps
```

| Severity | Scenario | Description | Status |
|---|---|---|---|
| | | | |

## Sign-off

- [ ] All 5 scenarios pass.
- [ ] No critical defects (P0/P1) outstanding.
- [ ] CSV from patient 57955439 archived to `/Volumes/Molloilab/Shu Nie/PVAT_Analysis/57955439_mmd_surfaces.csv` for later cohort comparison.
- [ ] Memory updated with verification result and remaining gaps for the
      66-patient batch run.
