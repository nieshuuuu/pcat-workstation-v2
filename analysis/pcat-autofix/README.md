# pcat-autofix

Post-hoc repair and re-measurement of the trees this workstation exports, plus a
reproduction of the Siemens pipeline to compare against.

Everything here reads an exported tree (`<root>/<patient>/{manifest,centerline/,lumen/,
fat/,views/}`) together with the source DICOMs, and rewrites the tree in place. Originals
are kept as `*.orig.json`; the tool always re-centres the operator's original path, so it
is idempotent and deletes nothing.

```bash
python pcat_autofix.py --all --measure-on 70 --second-energy 150
python test_pcat_syngo.py        # brute-force check of the min-cut construction
```

`viewer/index.html` is the browser viewer for the rewritten tree. Open it beside the tree
and point it at the folder.

**`NOTES.md` is the write-up** — what was wrong, what changed, the measurements, and the
decisions still open. Read that first.

## Why this exists

Two geometry defects made the exported FAI invalid rather than noisy:

1. **The operator's centerline can sit outside the vessel.** On P92 the proximal RCA path
   is 2.7–4.4 mm off the axis, in epicardial fat — the sample reads −83 HU where the
   lumen reads +450. The lumen finder cast 360 independent rays and stopped each at the
   first edge it met, so from a start point in fat most rays never met an edge and clamped
   at the 8 mm ceiling while a couple stopped at 0.7 mm. The contour that came out was a
   disc with a bite taken out of it, and the FAI ring was built on it. 28–44% of P92's
   rays were clamped and `detected` still reported success.
2. **The exported `normal` / `binormal` were not perpendicular to the vessel** (a median
   46–76° off) and did not reproduce `views/*.i16`, so the angular axis of the
   (arclength, angle) maps meant nothing.

Both are fixed here, and the FAI is now measured on the published convention — a shell
from the segmented **outer wall** outwards over a radial distance equal to the mean lumen
diameter — with the older lumen-keyed ring reported beside it so nothing changes meaning
invisibly.

## What belongs back in the workstation

The cheapest and most valuable piece: a **pre-flight gate**. Refuse any cross-section
whose centerline voxel is below roughly +150 HU on a 70 keV VMI, and report the fraction
of rays pinned at the search ceiling. Either check alone would have caught this before an
FAI was ever written out.

## Files

| | |
|---|---|
| `pcat_autofix.py` | re-centring by medialness minimum-cost path, a lumen boundary that cannot leak, outer-wall detection, the convention shell, histograms, the 70/150 keV decomposition |
| `pcat_syngo.py` | the Siemens pipeline reproduced — rising-edge-vetoed medialness, Viterbi re-centring, and the lumen as one MRF solved exactly by min-cut |
| `test_pcat_syngo.py` | brute-force verification of the min-cut against the global optimum |
| `vmi_completeness.txt` | per-patient VMI series completeness across the cohort |
| `figures/` | before/after, shell geometry, the two-patient comparison, syngo vs autofix |

## Note on data

Nothing patient-derived is committed here — no export trees, no volumes, no masks. The
`.gitignore` in this folder blocks them by pattern. `NOTES.md` refers to cases by their
study pseudonyms only.
