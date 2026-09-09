# What was wrong, what changed, and what still needs a decision

`python tools/pcat_autofix.py --all --measure-on 70 --second-energy 150`

Originals are kept as `*.orig.json` and `*_mask.orig.nii.gz`; the tool always re-centres
the operator's original path, so it is idempotent and nothing is deleted. About 10 s per
vessel once the volume is cached; the first read of a 352-slice series costs ~30 s.

---

## 1. The Pac-Man

P92's proximal RCA centerline sits **2.7–4.4 mm off the vessel axis, inside epicardial
fat**. The sample at arc 7.2 mm reads **−83 HU**; the lumen 3 mm away reads **+450 HU**.

The lumen finder cast 360 independent rays and stopped each at the first edge it met.
From a start point in fat, most rays never meet an edge and clamp at the 8 mm ceiling,
while the two rays that happen to point at the artery stop at 0.7 mm. That is the disc
with a bite out of it. It is not a display artefact — the FAI ring was built on it.

| as exported | rays clamped at the 8 mm ceiling | median lumen radius |
|---|---|---|
| P92 RCA | 28% | 1.90 mm |
| P92 LAD | 44% | 6.66 mm |
| P92 LCx | 40% | 6.23 mm |
| P106 (all three) | 0–2% | 1.6–1.9 mm |

A second, independent defect: the `normal` / `binormal` written into
`centerline/*.json` were **not perpendicular to the vessel** (a median 46–76° off) and
did not reproduce `views/*.i16`. Anything drawn from them — notably the angular axis of
the (arclength, angle) maps — meant nothing.

## 2. The fix, in three stages

**Re-centring.** A minimum-cost path through a medialness field inside a 7 mm tube
around the operator's path. Medialness is the radius of the maximal inscribed sphere in
the contrast column, so the cheapest path is by construction the one up the middle. The
operator's path stays a soft prior: this re-centres their vessel, it does not trace a new
one. Where the operator was already right (all of P106) it moves 0.15–0.25 mm.

The path is anchored to the outermost operator sample that still has a vessel near it,
not to their literal first and last click. Pinning it to the last click is what produced
P92's LAD tail: their last samples sit past the end of the opacified vessel, so the
cheapest path to that point left the artery and crossed several millimetres of nothing.

**A lumen boundary that cannot leak.** One cyclic dynamic program per cross-section over
the polar-unwrapped image, with a hard cap of 60° on the wall slope, solved exactly
around the θ wrap, then a second pass coupling the sections along the vessel. There is no
per-ray decision left to make. Climbing 0.7 → 8 mm now needs more circumference than
exists. The angular grid is 180 rays at δr = 0.05 mm, not 360 at 0.1 mm — at the coarser
radial step the smallest expressible constraint still permits a 73° wall, so the prior
had no teeth.

**An outer wall.** The same dynamic program against a cost that finds where tissue stops
being tissue and starts being fat, floored at the lumen. The adventitia has no contrast
edge — that is why the commercial packages use a different algorithm here — but for a
shell defined by an HU window, that crossing is the operational boundary anyway.
`wall_found_fraction` reports how many rays actually reached fat; the rest are a
least-bad guess and contribute nothing to the average.

Result on all six vessels: **0% of rays at the ceiling, 0% of centerline samples outside
the contrast column**, lumen diameters 2.5–3.7 mm, frames exactly orthonormal.

## 3. P92's LAD cannot be segmented past 30 mm, and that is the answer

There is a calcified spot at 26 mm, then **six millimetres with no lumen at all** and a
streak-artefact texture, and what lies beyond is a 6.4 mm-wide structure that is not the
LAD — a vein or a chamber. The operator's path ran parallel to it, 4.2–5.3 mm off, for
the whole distal third.

So the residual off-lumen fraction there was never a segmentation failure. The tool now
keeps the longest stretch of path that is genuinely inside the contrast column, bridging
gaps shorter than `lesion_gap_mm` (3 mm — a lesion belongs in the segment) and cutting at
longer ones. P92's LAD comes out as **29.8 mm usable, 10.2 mm dropped**, reported in
`qc.json → trimmed_mm`, and 0% off-lumen over what is kept. That is the honest number.

## 4. The FAI is now measured the published way

Primary shell: from the **outer vessel wall** outwards, over a radial distance equal to
the **mean lumen diameter** of the segment, gated −190…−30 HU (Antonopoulos 2017,
Oikonomou 2018, Kotanidis & Antoniades 2021). Reported alongside, so nothing changes
meaning invisibly:

- `shell_alt` — thickness scaled to the mean *outer* diameter, the other reading of the
  ambiguous phrase "a radial distance equal to the vessel diameter".
- `shell_gapped` — the same shell started 0.75 mm outside the wall. Worth looking at:
  the in-plane PSF (σ ≈ 0.375 mm) spills the ~540 HU lumen step into the first
  half-millimetre — 136 HU at 0.25 mm, 49 at 0.50, 12 at 0.75 — and because the PSF is
  fixed in mm while the shell scales with diameter, that spill is **vessel-size-coupled**
  (~19 HU of mean spill at a 2.5 mm vessel against ~10 HU at 4.5 mm). On these six
  vessels the gap moves the FAI 3–5 HU more negative.
- `legacy_ring` — the app's old lumen-keyed 1 mm gap + fixed 3 mm ring.

Every average is volume-weighted (∝ r). A polar grid is uniform in (θ, r) but the tissue
is not, and pericoronary fat has a real radial attenuation gradient.

The 20 mm depot is untouched and reported separately, as an epicardial-fat measure read
on its own terms.

## 5. Measured at 70 keV, and why that mattered

P92 was exported at **70 keV, Qr40f, QIR 2**; P106 at **65 keV, Bv40v, QIR 3**. Both
have a complete 70 keV series, so `--measure-on 70` harmonises the energy. It is not
cosmetic: P106's RCA FAI moved by **+3.3 HU** on the energy change alone.

The kernel and QIR still cannot be harmonised from what is on the share — P92 has only
Qr40f/Q2 and P106 only Bv40v/Q3. But that term turns out to be small at this kernel
sharpness: see §9, where it is measured directly rather than cited.

## 6. The reference-tissue question, answered

**Neither remote epicardial fat nor subcutaneous fat, as the primary.** The reasoning and
the evidence are in `Z:\Shu Nie\shu_pcat\tools\` history and summarised here:

A reference cancels a nuisance only if it *shares the nuisance* and *does not share the
biology*. Subcutaneous fat fails both: it sits at a different depth and path length, in a
different part of the FOV, with different beam hardening, no cardiac motion, and almost
no iodine (PCAT enhances ~22 HU, SAT under 5) — and published PCAT-minus-SAT work shows
ρ ≈ 0 and an SD that *rose* in every stratum on subtraction. It also tracks BMI, which
tracks CAC, so in an age+sex-matched but not BMI-matched design it is a mechanism for a
spurious positive. Remote epicardial fat passes the first test (ρ ≈ 0.7 with PCAT) but
the biology leak that would break it is 1.66 HU and the published estimates straddle that
value with inconsistent sign. It belongs in a **prespecified secondary as a covariate**,
never as a subtraction, and never as a ratio (HU has an arbitrary zero at water).

And the term neither depot can cancel: the reconstruction acts through the **width** of
the fat distribution rather than as an additive offset — measured here in P92, turning
QIR off moves σ by 4 HU and the FAI by only 1 HU (§9). Two windowed means carry two
different truncation biases, which subtract to a residual that varies patient to
patient.

**One correction to what I said earlier.** At 70 keV with σ ≈ 26 HU the window is *not*
truncating the fat peak — measured here, `frac_below_window` is **0.00% on all six
vessels**. The published mean is not biased at this operating point. What is fragile is
its **gain**: d(FAI)/dσ runs −0.05 to −0.35 HU per HU, so a noise difference between arms
manufactures an FAI difference from nothing. That is visible in these two patients: P92
sits at σ ≈ 22–28 (Qr40f/Q2) and P106 at σ ≈ 30–34 (Bv40v/Q3), and re-evaluating at a
fixed σ moves them in **opposite** directions by about 6 HU combined.

So what is now computed per shell, from the unwindowed histogram:

- `hist_mm3` — the whole distribution, −250…+150 HU in 2 HU bins, volume-weighted.
- `fit_mean_hu`, `fit_sigma_hu` — a Gaussian fit over [µ−3σ, µ+1.2σ], seeded from the
  modal bin. Right-clipped because a shell reaching into myocardium has a tissue shoulder
  there; **not** left-flank-only, which measures worse (kernel sensitivity +3.65 HU
  against +0.3–0.6 for the clipped-symmetric fit).
- `hu_mean_sigma30` — the published statistic re-evaluated at a fixed σ = 30, which keeps
  a number reviewers recognise while removing the noise term.
- `frac_below_window` / `frac_above_window` / `frac_in_window` — the quality gates. A
  between-arm imbalance in these invalidates any windowed-mean comparison.
- `fit_failed` — set when there is no fat peak to fit. P106's LCx shell is only 29% fat;
  its FAI is an average over whatever survived the window and should be treated as
  unusable rather than low.

## 7. The 70/150 keV pair

Siemens VMI is rank-2, so two energies exhaust the spectral information; a third adds
only correlated noise. `--second-energy 150` samples the same shell voxels at both
energies and reports `hu_low_minus_high`, `adipose_fraction` and `iodine_mg_ml` from the
gecatsim/NIST basis (adipose −104.09 / −78.85 HU, soft tissue +43.07 / +40.62, iodine
25.861 / 4.458 HU per mg/mL, with volume closure). The solve is linear, so it is applied
to the shell mean rather than voxelwise.

P106 gives adipose fraction 0.87–0.93 and iodine 0.29–0.84 mg/mL — consistent with the
+18 HU / 3.7%-of-aortic PCAT iodine measured separately.

Keep HU at 70 keV as the primary. The difference channel buys keV-invariance, not
sensitivity: contrast per unit adipose fraction is −147 HU at 70 keV but only −28 HU on
(HU₇₀ − HU₁₅₀), a 5.4× loss, and with 21 pairs sensitivity is the binding constraint.
What it *does* buy is separating "the fat has less lipid" from "the fat is holding more
contrast" from "the scan was three seconds later in the bolus" — which is the one term no
reference tissue can cancel.

## 8. Two data problems found on the share

**P92's 150 keV series is missing 229 of 352 slices** — 123 files in three blocks,
covering 12% of the RCA and 0% of the LAD and LCx. Read blind, SimpleITK takes the
spacing from the first two files and places everything after a gap at the wrong z, so the
tool now checks slice positions and suppresses the spectral result rather than returning
confident nonsense.

A sweep of all 131 patients under `Egor_true rec` (`tools/vmi_completeness.txt`):

| | |
|---|---|
| 150 keV complete | 129 |
| 150 keV short | **P92** (123/352), **P108** (258/340) |
| 150 keV missing | 0 |
| other short series | P101 (70 keV 325/383, 65 keV 304, 100 keV 334) |

So the spectral endpoint is viable cohort-wide; three patients need their recons
re-exported.

## 9. How much do the two patients actually differ?

Measured by holding the segmentation geometry fixed and re-sampling only the attenuation,
so the reconstruction is compared against itself rather than against a different contour.

**The energy, measured in P106** (six VMI energies, one geometry, one kernel/QIR):

| keV | RCA | LAD | LCx |
|---|---|---|---|
| 40 | −109.2 | −102.5 | −89.5 |
| 65 | −83.5 | −82.3 | −73.4 |
| **70** | **−80.2** | **−80.1** | **−71.7** |
| 100 | −71.4 | −73.6 | −67.7 |
| 140 | −68.4 | −71.5 | −67.1 |
| 150 | −68.1 | −71.2 | −67.0 |

Harmonising 65 → 70 keV moves P106 by **+1.7 to +3.3 HU**. The local slope at 70 keV is
0.16–0.35 HU per keV, which brackets the 0.29 the published anchors give.

**The reconstruction, measured in P92** (same patient, same 70 keV, same Qr40f; QIR 2
against FBP, i.e. QIR off):

| | FAI QIR2 | FAI FBP | ΔFAI | σ QIR2 | σ FBP | Δσ |
|---|---|---|---|---|---|---|
| RCA | −91.7 | −93.4 | −1.7 | 24.7 | 28.9 | +4.2 |
| LAD | −83.2 | −84.0 | −0.8 | 22.9 | 26.8 | +3.9 |
| LCx | −73.2 | −74.7 | −1.4 | 28.8 | 32.6 | +3.9 |

So on this kernel the reconstruction moves the **width** of the fat distribution about
three times as much as it moves the **position**. The 19 HU figure quoted earlier was for
soft kernels — the same source has Bv36 at +12 HU and Bv44 at +0.6 HU, and a 40-family
kernel sits at the small end. **The earlier warning that kernel/QIR is the largest
technical term is not supported at Qr40/Bv40.**

**The gap between the patients:**

| | as exported (70 vs 65 keV) | both at 70 keV | fitted peak, both at 70 |
|---|---|---|---|
| RCA | −8.2 | **−11.5** | −20.7 |
| LAD | −1.0 | **−3.2** | −9.8 |
| LCx | +0.2 | **−1.6** | no peak to fit |

Harmonising the energy did not shrink the gap, it **widened** it — P106 was the one at
65 keV, and moving it to 70 keV made it less negative.

**And the confound that has to be reported with it.** Across the six vessels, FAI tracks
how much of the shell is fat: r = −0.87, about **−2.9 HU per 10 percentage points**. P92's
RCA shell is 97% fat; P106's is 69%. That 28-point difference alone predicts roughly 8 HU
of the 11.5 HU RCA gap.

With two patients this cannot be separated, and the vessel ordering (RCA more negative
than LAD than LCx) is a known anatomical fact that inflates the correlation on its own.
It is a hypothesis for the real cohort, not a finding. The action is small: **record shell
fat fraction per vessel and check it is balanced between arms** — it is already in
`fat/<vessel>.json → shell.fat_fraction`. If it is imbalanced, part of any FAI difference
is how much fat surrounds the artery rather than what that fat is made of.

### The precision floor

The working grid was 0.35 mm and is now **0.25 mm**. Measured on P92: the shell mean moves
0.9 HU between 0.35 and 0.25 mm and then under 0.1 HU down to 0.15 mm, so 0.25 is where it
converges. Re-running the tool reproduces its own numbers exactly; changing only the crop
origin, which shifts the resampling phase by a sub-voxel amount, moves the FAI by
0.2–0.4 HU. **Differences under about 0.5 HU in any table here are not meaningful.**

The fitted peak barely moves across that whole grid sweep (0.4 HU from 0.45 to 0.15 mm)
where the windowed mean moves 0.9 HU — one more argument, measured here rather than
cited, for reporting the peak.

## 10. The Siemens pipeline, reproduced

`pcat_syngo.py` implements what syngo.VIA's coronary analysis actually does, so the two
can be compared on the same vessels rather than argued about. `test_pcat_syngo.py`
checks the hard part against brute force.

**Medialness with a rising-edge veto** (Gulsun & Tek). For every in-plane direction the
profile is differentiated at two scales; a boundary at radius R is credited only if the
ray reaches it *without first crossing a rising edge*, and the response is divided by the
strongest falling edge on the same ray. Nothing is thresholded anywhere, and the result
is contrast- and dose-independent by construction.

It behaves exactly as advertised. Evaluated on a grid in the operator's own plane on
P92's RCA, where the seed sits in fat:

| section | m at the operator's seed | m at the peak | distance to the peak |
|---|---|---|---|
| 0 | 0.13 | 0.98 | 2.5 mm |
| 7 | 0.27 | 1.00 | 3.1 mm |
| 14 | 0.31 | 0.98 | 4.4 mm |
| 29 | 0.40 | 0.98 | 2.7 mm |
| 50 (operator already right) | 0.93 | 0.98 | **0.2 mm** |

**One globally optimal surface by min-cut.** Every (section, ray) radius is a label in a
Markov random field with convex Huber priors coupling angular and longitudinal
neighbours, solved exactly through Ishikawa's construction. Two things make it practical:
a Huber prior is linear past its knee, so the second difference that sets the arc
capacities is identically zero beyond it and the graph is sparse; and the theta wrap is
just another pair of neighbours, where a dynamic program has to close it by hand.

264k nodes and about 7M arcs solve in **2.0 s**.

**Two silent bugs, both worth knowing about if anyone reimplements this:**

1. `scipy.sparse.csgraph.maximum_flow` accepts an int64 matrix without complaint and
   truncates it to int32 internally. A sentinel of `1 << 40` becomes a capacity of
   **zero**, every closure constraint disappears, and the solver returns a flow of 0 with
   no error at all. The symptom was every column pinned at the maximum radius.
2. The Ishikawa pairwise arcs do **not** realise `g(a-b)`. They realise
   `g(a-b) + f(a) + f(b)` for a separable `f` that must be subtracted back out of each
   column's unary cost, once per neighbour pair. Miss it and the energy is wrong by a
   term that still produces plausible-looking contours.

Neither announced itself. Both were found by brute-forcing small instances - 54
configurations across sizes, Huber knees and prior weights, all now matching the global
optimum exactly.

**What the comparison says.** On five of the six vessels the two pipelines agree to a
median axis difference of **0.17-0.27 mm** and a median lumen-diameter difference of
**0.116 mm** over 404 cross-sections - two completely different centring criteria
(distance transform of a thresholded volume against a threshold-free rising-edge veto)
and two different surface solvers landing in the same place.

The sixth is P92's LAD, and it is instructive. The operator's path there is ~6 mm off.
The Gulsun-Tek stage searches a +-3 mm grid around the *current* path and climbs the
nearest medialness ridge, so it settles on a spurious ridge (m = 0.72) and stays: 2, 4
and 6 iterations all finish 6.6-6.8 mm from the artery. It is a local optimum, not slow
convergence. Widening the grid to 7 mm recovers it only partly (3.0 mm off, 49% still
outside contrast) and costs 24 s. `pcat_autofix` traces a *global* minimum-cost path
through a 7 mm tube, so the same 6 mm error is recoverable.

That is not a defect in the Siemens design. In the product the centring stage never sees
an operator's path - it polishes the output of an automatic tree extraction that is
already close. The lesson for this tree is the narrower one: **a local ridge-climbing
re-centring is not enough when the input is a hand-drawn path that can be six
millimetres out**, which on this data it can.

## 11. The study cohort is not protocol-confounded, and P92 was never in it

`analysis/matched_pairs_RCA_primary.csv` resolves to 42 patients across the 21 frozen
pairs. Crossed against the protocol scan:

- **all 42 slots are on one protocol** — 120 kVp, Bv64v / QIR 2, **1.0 mm**
- **21 of 21 pairs have case and control on the same protocol**

So the kVp, kernel and QIR terms cancel in the pairing. The five-protocol split in the
cohort is real but it does not reach the analysis that was frozen. **P92 and P106 appear
in neither arm of any pair** — they are development cases, and the fact that they sit on
different protocols was never a study problem.

Which also settles the question of swapping P92 for one of the five Qr40v patients:
there is no point. The Qr40v group is 0.4 mm and Qr-family; the study is 1 mm and
Bv64. If a development case should look like the study, it has to come from the 42.

**And none of the 42 has any sub-millimetre reconstruction.** They carry 70 keV and
150 keV at 1 mm and 3 mm only. Every piece of segmentation work here has been done at
0.4 mm, so that gap is the one that matters.

### What 1 mm does

P106's 0.4 mm volume blurred along z to what a 1 mm slice would see, geometry untouched,
whole pipeline re-run:

| | FAI | fitted peak | σ | lumen diameter |
|---|---|---|---|---|
| RCA at 1 mm | **+0.65** | −1.05 | −2.20 | +0.004 mm |
| LAD at 1 mm | **+1.11** | −1.71 | −3.39 | −0.014 mm |
| LCx at 1 mm | **+1.03** | +4.00 | — | +0.012 mm |
| RCA at 3 mm | **+2.87** | −3.55 | −8.64 | −0.021 mm |
| LAD at 3 mm | **+5.66** | −2.09 | −10.20 | +0.015 mm |
| LCx at 3 mm | **+5.75** | +4.00 | — | +0.099 mm |

Three things worth taking from that.

**The segmentation is fine at 1 mm.** Lumen diameter moves by under 0.02 mm and the
wall-found fraction does not move at all. The pipeline can be pointed at the study data
as it stands.

**The windowed FAI moves the opposite way from the fat it is measuring.** At 1 mm the FAI
reads about 1 HU *less* negative while the fitted peak goes 1–1.7 HU *more* negative. The
whole of the FAI change is the width: σ drops 2–3 HU, and at d(FAI)/dσ ≈ −0.2 to −0.35
that is +0.6 to +1.1 HU on its own. This is the gain argument of §6, measured here rather
than cited.

**3 mm is not usable.** +2.9 to +5.8 HU, the size of the effect being hunted.

Caveat: this is a boxcar blur of a thin reconstruction, not a real 1 mm reconstruction,
which would also carry a different z-kernel. Direction and rough magnitude, not exact
numbers.

### The QIR term, measured properly this time

The five Qr40v patients each carry QIR-2 and QIR-off reconstructions of every energy at
matched 0.4 mm slices — the controlled pair P92 could not provide, because its FBP series
is 0.8 mm thick and that comparison changed two things at once. Measured on identical
voxels of deep mediastinal fat, no centerline needed, FBP minus QIR 2:

| n = 4 patients | windowed FAI | fitted peak | σ |
|---|---|---|---|
| deep fat | **−3.07 ± 0.91** | −2.28 ± 1.19 | **+15.61 ± 1.66** |
| subcutaneous | −1.59 ± 0.55 | −0.57 ± 0.44 | +13.91 ± 1.79 |

So turning QIR off moves the *width* five to seven times as much as the *position*, and
the earlier figure of ~1 HU from P92 was too small — the 0.8 mm FBP slices had suppressed
most of the noise difference. QIR 2 against QIR 3, which is what separates P92 from P106,
is one level rather than off-to-2 and is smaller again. In the frozen study it cancels
anyway, since all 42 patients are QIR 2.

Two data notes: P87's `70 kev FBP 0.4 mm` folder reports kernel `Qr40v/2` in the DICOM,
i.e. it is not FBP, and it was excluded. And there is **no clean kernel-family pair
anywhere in the cohort** — no patient has Qr and Bv at matched sharpness, QIR, energy and
thickness — so the Qr↔Bv term cannot be measured in this data and would have to come from
a re-reconstruction.

## 12. Still open

- **The ostium is not marked.** `ostium_fraction` is `null` in every manifest, so the
  RCA 10–50 mm convention cannot be applied — arclength is measured from wherever the
  operator started clicking. Detecting it from the images was tried and is not reliable:
  a geodesic walk from the coronary into the nearest large contrast pool measures the
  distance to that pool's *core*, not to the ostium, and the numbers came back 19–35 mm
  with no anatomical clustering. This is a one-click job in the workstation and it is the
  only thing left blocking full convention compliance.
- **P92's LCx is 34.8 mm** against the operator's 40.1 mm: their path wandered, so it was
  longer than the true axis between the same two stations.
- **Kernel and QIR are confounded with patient in P92 vs P106**, but not in the study:
  all 42 patients in the frozen pairs share one protocol (§11). No re-reconstruction is
  needed for the analysis as designed.
- **Do not quote −70.1 HU.** The viewer draws it for orientation and labels it as a
  120 kVp energy-integrating cut-off that does not transfer to this reconstruction.
