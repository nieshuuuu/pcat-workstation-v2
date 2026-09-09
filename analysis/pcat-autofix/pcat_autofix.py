#!/usr/bin/env python3
"""
pcat_autofix - re-centre a hand-drawn coronary centerline onto the vessel axis and
re-cut the lumen with a boundary that cannot leak, then rewrite a shu_pcat tree.

Why this exists
---------------
The exported trees were built by casting 360 independent radial rays from each
centerline sample and stopping each one at the first edge it met.  That works only
while the sample sits inside the contrast column.  On P92 the operator's proximal RCA
centerline sits 2.7-4.4 mm off the axis, in epicardial fat (-83 HU at the sample, vs
+450 HU in the lumen), so most rays never meet an edge and clamp at the 8 mm search
ceiling while the two rays that happen to point at the vessel stop at 0.7 mm.  The
contour that comes out is a disc with a bite taken out of it.

Two things are fixed here, in this order, because the second depends on the first:

  1. Re-centring.  A minimum-cost path is traced through a medialness field inside a
     tube around the operator's path.  Medialness is the radius of the maximal
     inscribed sphere in the contrast column - the same quantity VMTK's Voronoi
     centerline maximises - so the cheapest path is by construction the one that runs
     up the middle of the lumen.  The operator's path is kept as a soft prior, so this
     re-centres their vessel rather than re-tracing a new one.

  2. A boundary that cannot leak.  The lumen contour in each section is found by one
     dynamic program over the polar-unwrapped section with a hard cap on how far the
     radius may move between neighbouring angles, solved exactly around the cyclic
     wrap.  A ray that finds no edge can no longer run away on its own: it is pulled
     to whatever its neighbours agreed on.  A second pass couples the sections along
     the vessel the same way.

Everything else in the tree (gap, ring, HU gate, depot growth) keeps the app's own
definitions so the numbers stay comparable; only the geometry they are measured on
changes.

Usage
-----
    python pcat_autofix.py --root "Z:/Shu Nie/shu_pcat" --patient P92
    python pcat_autofix.py --root "Z:/Shu Nie/shu_pcat" --all

Originals are kept beside the rewritten files as <name>.orig.json.  Nothing is
deleted.  Pass --dry-run to print the QC table without writing.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import time
from dataclasses import dataclass
from math import erf, pi, sqrt

import numpy as np
import SimpleITK as sitk
from scipy import ndimage as ndi
from scipy.interpolate import splev, splprep
from scipy.ndimage import map_coordinates
from scipy.spatial import cKDTree
from skimage.graph import MCP_Geometric

VESSELS = ("RCA", "LAD", "LCx")


# ---------------------------------------------------------------------------
# tunables.  Every one of these is a length in mm or a dimensionless weight; none
# are pixel counts, so they carry over to a different scanner unchanged.
# ---------------------------------------------------------------------------


@dataclass
class Params:
    iso_mm: float = 0.25          # working grid.  Finer than the 0.83 mm in-plane
                                  # sampling on purpose: the medial axis of a 3.5 mm
                                  # vessel is not resolvable on a 0.83 mm lattice.
                                  # 0.25 is where the FAI stops moving - measured on
                                  # P92, the shell mean shifts 0.9 HU between 0.35 and
                                  # 0.25 mm and then under 0.1 HU down to 0.15 mm.  The
                                  # fitted peak is far less sensitive (0.4 HU over the
                                  # whole range), which is one more reason to prefer it.
    crop_margin_mm: float = 14.0
    tube_mm: float = 7.0          # how far the re-centred path may stray from the
                                  # operator's.  Larger than any plausible operator
                                  # error, smaller than the distance to the chambers.
    prior_mm: float = 6.0         # soft pull back towards the operator's path
    medial_ref_mm: float = 2.0    # medialness is capped here, so a cardiac chamber is
                                  # not more attractive than a normal coronary
    medial_floor_mm: float = 0.15
    snap_mm: float = 5.0          # how far an endpoint may be moved onto the axis
    snap_min_medial_mm: float = 0.6   # and how fat the thing it lands on must be, so
                                      # that a speck of noise two voxels away cannot
                                      # outbid the vessel four millimetres away
    smooth_mm: float = 0.8        # centerline smoothing sigma along arclength.  Small:
                                  # this removes voxel-scale staircase only, and a
                                  # coronary's own bend radius is 20-50 mm.
    lesion_gap_mm: float = 3.0    # a dark stretch shorter than this is a lesion and
                                  # stays in the segment; a longer one is an occlusion,
                                  # an artefact, or the trace jumping to another vessel
    min_usable_mm: float = 12.0   # below this there is no point reporting an FAI
    max_lumen_diameter_mm: float = 5.0  # sanity ceiling: anything wider is a chamber or
                                        # a vein, not a native coronary

    r_max_mm: float = 4.5         # lumen search ceiling.  No native coronary lumen has
                                  # a 4.5 mm radius; the old 8 mm ceiling is what let
                                  # a failed ray look like a plausible measurement.
    dr_mm: float = 0.05
    n_theta: int = 180            # 180 rays on a 1.75 mm lumen is one sample every
                                  # 0.06 mm of arc, still 13x finer than the 0.83 mm
                                  # pixel.  Going to 360 does not add information, it
                                  # only makes the angular prior below toothless.
    wall_slope_deg: float = 60.0  # the steepest the lumen wall may run in the polar
                                  # plane.  This is the hard cap that makes a leak
                                  # geometrically impossible rather than merely
                                  # expensive: climbing 0.7 -> 8 mm now needs more
                                  # circumference than exists.  Expressed as an angle
                                  # because the radial step allowed between adjacent
                                  # rays has to grow with radius - the arc between two
                                  # rays does.
    s_step_mm: float = 0.45       # hard cap on |dr| between adjacent sections
    edge_win_mm: float = 0.75     # half-width of the inside/outside contrast window
    w_inside: float = 0.6         # penalty for declaring lumen where it is not bright
    w_smooth_s: float = 1.2       # weight of the along-vessel coupling, 2nd pass

    calc_hu: float = 700.0        # above this is calcium, not lumen: it must not drag
                                  # the lumen estimate up, and it blocks a ray
    fat_gate: tuple = (-190.0, -30.0)

    # the outer vessel wall, and the shell the published FAI is measured in
    wall_max_mm: float = 6.0      # no coronary outer wall sits further out than this
    wall_dr_mm: float = 0.10      # coarser than the lumen's 0.05 mm on purpose: the
                                  # adventitia is thinner than one 0.83 mm pixel, so a
                                  # finer radial grid buys resolution the data has not
    wall_min_gap_mm: float = 0.15
    w_wall_out: float = 1.0       # penalty for the tissue outside the wall not being fat
    w_wall_in: float = 0.6        # penalty for the tissue inside the wall being fat
    shell_scale: str = "lumen_diameter"   # what "a radial distance equal to the vessel
                                          # diameter" is measured on.  The published
                                          # wording is ambiguous; both readings are
                                          # reported, this one is the primary.
    shell_gap_mm: float = 0.0     # the convention starts the shell at the outer wall.
    shell_gap_alt_mm: float = 0.75    # A 0.75 mm gap is reported beside it because the
                                      # in-plane PSF (sigma about 0.375 mm) spills the
                                      # 540 HU lumen step into the first half-millimetre
                                      # of the shell: 136 HU at 0.25 mm, 49 at 0.50,
                                      # 12 at 0.75.  The PSF is fixed in mm while the
                                      # shell scales with diameter, so that spill is
                                      # vessel-size-coupled - 19 HU of mean spill at a
                                      # 2.5 mm vessel against 10 HU at 4.5 mm - and a
                                      # 0.75 mm gap collapses the spread to 0.4 HU.
    ref_sigma_hu: float = 30.0    # the fixed width the standardised FAI is evaluated at
    kev_pair: tuple = (70, 150)   # the two energies the decomposition basis is for
    hist_lo: float = -250.0       # the unwindowed histogram of the shell, so a reader
    hist_hi: float = 150.0        # can see the whole fat peak and not just the part
    hist_bin: float = 2.0         # the -190/-30 window happened to keep


# ---------------------------------------------------------------------------
# volume handling
# ---------------------------------------------------------------------------


def resolve_series_dir(manifest, images_root):
    """Map the manifest's Mac-side path onto whatever this machine calls that share."""
    src = manifest.get("source_patient_dir", "")
    measured = manifest.get("measured_on", "")
    for s in manifest.get("series", []):
        if s.get("name") == measured:
            src = s.get("path", src)
            break
    tail = src.replace("\\", "/")
    low = tail.lower()
    marker = "/volumes/molloilab/"
    if low.startswith(marker):
        tail = tail[len(marker):]
    cand = os.path.join(images_root, *tail.split("/"))
    if os.path.isdir(cand):
        return cand
    raise FileNotFoundError(
        "could not find the images for '%s'\n  tried: %s\n"
        "  pass --images with the folder that holds 'Hamid/Egor_true rec'" % (src, cand))


def find_energy_series(vmi_dir, kev):
    """The non-FBP series at this keV beside the one being measured."""
    want = str(kev).lower().replace("kev", "").strip()
    best = None
    for nm in sorted(os.listdir(vmi_dir)):
        d = os.path.join(vmi_dir, nm)
        if not os.path.isdir(d) or "fbp" in nm.lower():
            continue
        m = re.search(r"(\d+)\s*ke?v", nm.lower())
        if m and m.group(1) == want:
            n = sum(1 for f in os.listdir(d) if f.lower().endswith(".dcm"))
            if best is None or n > best[0]:
                best = (n, d)
    return best[1] if best else None


def load_volume(series_dir, cache_dir):
    """Read the series once, then keep a NIfTI beside it - re-reading 352 DICOMs on
    every run is most of the wall clock and none of the interesting part."""
    os.makedirs(cache_dir, exist_ok=True)
    key = series_dir.replace("\\", "/").strip("/").replace("/", "_").replace(":", "")
    cache = os.path.join(cache_dir, key + ".nii.gz")
    if os.path.exists(cache):
        return sitk.ReadImage(cache)
    r = sitk.ImageSeriesReader()
    ids = r.GetGDCMSeriesIDs(series_dir)
    if not ids:
        raise FileNotFoundError("no DICOM series in " + series_dir)
    r.SetFileNames(r.GetGDCMSeriesFileNames(series_dir, ids[0]))
    img = r.Execute()
    sitk.WriteImage(img, cache)
    return img


def series_z(series_dir):
    """The slice positions of a series, headers only.

    Worth the extra pass: a folder can be missing slices and still read as a stack, and
    SimpleITK will then take the spacing from the first two files and place everything
    after the gap at the wrong z.  P92's 150 keV folder holds 123 of 352 slices in three
    blocks; sampled blind it would have returned confident nonsense.
    """
    r = sitk.ImageSeriesReader()
    ids = r.GetGDCMSeriesIDs(series_dir)
    if not ids:
        return np.array([]), []
    files = r.GetGDCMSeriesFileNames(series_dir, ids[0])
    rd = sitk.ImageFileReader()
    rd.LoadPrivateTagsOff()
    z = []
    for f in files:
        rd.SetFileName(f)
        rd.ReadImageInformation()
        z.append(rd.GetOrigin()[2])
    z = np.array(sorted(z))
    if len(z) < 2:
        return z, [(z[0], z[0])] if len(z) else []
    step = float(np.median(np.diff(z)))
    brk = list(np.flatnonzero(np.diff(z) > 1.9 * step))
    blocks, lo = [], 0
    for b in brk + [len(z) - 1]:
        blocks.append((float(z[lo]), float(z[b])))
        lo = b + 1
    return z, blocks


class Grid:
    """An isotropic crop around one vessel, plus the world <-> index maps."""

    def __init__(self, img, pts_mm, p):
        arr = sitk.GetArrayFromImage(img).astype(np.float32)   # [z, y, x]
        sp = np.array(img.GetSpacing(), float)                 # (x, y, z)
        org = np.array(img.GetOrigin(), float)
        d = np.array(img.GetDirection(), float).reshape(3, 3)
        if not np.allclose(d, np.eye(3), atol=1e-6):
            raise NotImplementedError("oblique acquisitions are not handled")

        idx = (pts_mm - org) / sp                              # (i, j, k)
        m = p.crop_margin_mm / sp
        lo = np.maximum(np.floor(idx.min(0) - m).astype(int), 0)
        hi = np.minimum(np.ceil(idx.max(0) + m).astype(int) + 1,
                        [arr.shape[2], arr.shape[1], arr.shape[0]])
        sub = arr[lo[2]:hi[2], lo[1]:hi[1], lo[0]:hi[0]]

        self.iso = p.iso_mm
        self.origin = lo * sp + org                            # world coord of vox 0
        self.vol = ndi.zoom(sub, (sp[2] / self.iso, sp[1] / self.iso,
                                  sp[0] / self.iso), order=1)  # [z, y, x] isotropic

    valid_z = None                 # (z0, z1) blocks the source series really holds

    def covers(self, w):
        """Fraction of world points that fall on a slice the source series actually has.

        Interpolating across a missing block returns a number, and it is meaningless.
        """
        if not self.valid_z:
            return 1.0
        z = np.atleast_2d(w)[:, 2]
        ok = np.zeros(len(z), bool)
        for a, b in self.valid_z:
            ok |= (z >= a - 0.5) & (z <= b + 0.5)
        return float(ok.mean())

    def to_vox(self, w):
        """world (x, y, z) mm -> grid index (z, y, x), float."""
        return ((np.atleast_2d(w) - self.origin) / self.iso)[:, ::-1]

    def to_world(self, v):
        """grid index (z, y, x) -> world (x, y, z) mm."""
        return np.atleast_2d(v)[:, ::-1] * self.iso + self.origin

    def sample(self, w):
        """Trilinear HU at world points of any shape (..., 3)."""
        w = np.asarray(w, float)
        v = self.to_vox(w.reshape(-1, 3)).T                    # (3, N) as z, y, x
        return map_coordinates(self.vol, v, order=1,
                               mode="nearest").reshape(w.shape[:-1])


# ---------------------------------------------------------------------------
# 1. re-centring
# ---------------------------------------------------------------------------


def estimate_hu(g, dist_to_path, p, tube_mm):
    """Lumen and background HU, read off the neighbourhood of the operator's path.

    The lumen estimate is a high percentile of everything inside the tube with calcium
    excluded, because the tube is mostly fat and the contrast column is its bright
    tail.  The background is the median of a shell that is past the wall but still
    peri-vascular, which on a coronary is fat.  The threshold is the midpoint, i.e.
    the usual full-width-half-maximum edge criterion for CT lumen sizing.
    """
    hu = g.vol[dist_to_path <= tube_mm]
    hu = hu[hu < p.calc_hu]
    lumen = float(np.percentile(hu, 99.0))
    shell = (dist_to_path > 3.0) & (dist_to_path <= 6.5)
    bg = float(np.median(g.vol[shell])) if shell.any() else -80.0
    return lumen, bg, 0.5 * (lumen + bg)


def shell_hu(g, pts, r_in, r_out):
    """Every voxel between two radii of a path, as a flat array of HU."""
    occ = np.zeros(g.vol.shape, bool)
    ii = np.clip(np.round(g.to_vox(pts)).astype(int), 0, np.array(g.vol.shape) - 1)
    occ[ii[:, 0], ii[:, 1], ii[:, 2]] = True
    d = ndi.distance_transform_edt(~occ, sampling=g.iso)
    return g.vol[(d > r_in) & (d <= r_out)]


def recentre(g, pts_mm, p, tube_mm=None, thr_override=None):
    """Trace the cheapest medialness path inside a tube around the operator's path."""
    tube_mm = p.tube_mm if tube_mm is None else tube_mm
    vox = g.to_vox(pts_mm)
    occ = np.zeros(g.vol.shape, bool)
    ii = np.clip(np.round(vox).astype(int), 0, np.array(g.vol.shape) - 1)
    occ[ii[:, 0], ii[:, 1], ii[:, 2]] = True
    dpath = ndi.distance_transform_edt(~occ, sampling=g.iso)

    lumen_hu, bg_hu, thr = estimate_hu(g, dpath, p, tube_mm)
    if thr_override is not None:
        thr = float(thr_override)

    smooth = ndi.gaussian_filter(g.vol, 0.5 / g.iso)
    lum = ndi.binary_opening(smooth > thr, np.ones((3, 3, 3)))
    medial = ndi.distance_transform_edt(lum, sampling=g.iso)   # inscribed radius, mm

    m = np.clip(medial, p.medial_floor_mm, p.medial_ref_mm)
    cost = (p.medial_ref_mm / m) ** 2 * (1.0 + (dpath / p.prior_mm) ** 2)
    cost[dpath > tube_mm] = 1e6

    def snap(v):
        k, j, i = np.round(v).astype(int)
        R = int(round(p.snap_mm / g.iso))
        sl = (slice(max(k - R, 0), k + R + 1),
              slice(max(j - R, 0), j + R + 1),
              slice(max(i - R, 0), i + R + 1))
        zz, yy, xx = np.mgrid[sl]
        dd = np.sqrt(((zz - k) * g.iso) ** 2 + ((yy - j) * g.iso) ** 2
                     + ((xx - i) * g.iso) ** 2)
        # The nearest thing that is actually a vessel, rather than the strongest thing
        # in range: a cardiac chamber must not outbid the artery the operator pointed
        # at, and neither must a two-voxel speck of noise that happens to be closer.
        med = medial[sl]
        ok = (dd <= p.snap_mm) & (med >= p.snap_min_medial_mm)
        if ok.any():
            w = np.unravel_index(int(np.argmax(np.where(ok, -(dd - 1.5 * med), -9e9))),
                                 dd.shape)
            return (int(zz[w]), int(yy[w]), int(xx[w])), float(dd[w]), float(med[w]), True
        w = np.unravel_index(int(np.argmax(np.where(dd <= p.snap_mm,
                                                    med - 0.35 * dd, -9e9))), dd.shape)
        return (int(zz[w]), int(yy[w]), int(xx[w])), float(dd[w]), float(med[w]), False

    def anchor(order):
        """Pin the path to the outermost operator sample that still has a vessel near it.

        Pinning it to the literal first and last click is what produced P92's LAD tail:
        the operator's last few samples sit past the end of the opacified vessel, so the
        cheapest path to that point leaves the artery and crosses several millimetres of
        nothing to reach it.  Walking inwards until the anchor is a real vessel keeps
        the whole path on the artery, and gives up only the samples that had nothing
        under them in the first place.
        """
        limit = max(1, len(order) // 3)
        for n, s in enumerate(order[:limit]):
            v, moved, med_v, ok = snap(vox[s])
            if ok:
                return v, moved, med_v, n
        v, moved, med_v, _ = snap(vox[order[0]])
        return v, moved, med_v, 0

    a, moved_a, med_a, skip_a = anchor(range(len(vox)))
    b, moved_b, med_b, skip_b = anchor(range(len(vox) - 1, -1, -1))

    mcp = MCP_Geometric(cost, sampling=(g.iso,) * 3)
    mcp.find_costs([a], [b])
    world = g.to_world(np.array(mcp.traceback(b), float))

    info = dict(lumen_hu=lumen_hu, bg_hu=bg_hu, threshold_hu=thr, tube_mm=tube_mm,
                start_moved_mm=moved_a, end_moved_mm=moved_b,
                start_medial_mm=med_a, end_medial_mm=med_b,
                anchor_skipped=[skip_a, skip_b],
                anchors=(pts_mm[skip_a], pts_mm[len(pts_mm) - 1 - skip_b]))
    return world, medial, thr, info


def resample_centerline(world, n, p, ends=None):
    """Smooth the voxel-stepped path and lay n samples on it at equal arclength.

    `ends` is the operator's first and last sample.  The re-centred path is trimmed to
    the feet of those two points, so the segment covers the same stretch of vessel the
    operator chose - otherwise snapping an endpoint sideways onto the axis can also
    slide it a millimetre or two along, and the reported segment length drifts.
    """
    step = np.linalg.norm(np.diff(world, axis=0), axis=1)
    w = world[np.r_[True, step > 1e-6]]
    if len(w) < 4:
        raise ValueError("path too short to smooth")

    s = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(w, axis=0), axis=1))]
    if p.smooth_mm > 0:
        w = ndi.gaussian_filter1d(w, p.smooth_mm / max(np.median(np.diff(s)), 1e-6),
                                  axis=0, mode="nearest")

    # arclength re-parameterisation, twice, because smoothing shortens the path
    for _ in range(2):
        s = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(w, axis=0), axis=1))]
        tck, _ = splprep(w.T, u=s / s[-1], s=0, k=3)
        w = np.array(splev(np.linspace(0, 1, max(n * 8, 400)), tck)).T

    s = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(w, axis=0), axis=1))]
    a, b = 0.0, s[-1]
    if ends is not None:
        a = s[int(np.argmin(np.linalg.norm(w - ends[0], axis=1)))]
        b = s[int(np.argmin(np.linalg.norm(w - ends[1], axis=1)))]
        if b < a:
            a, b = b, a
    target = np.linspace(a, b, n)
    out = np.column_stack([np.interp(target, s, w[:, d]) for d in range(3)])
    return out, target - target[0]


def trim_to_lumen(g, pts, arc, thr, p):
    """Keep the longest stretch of the path that is actually inside the contrast column.

    A short dark patch is a lesion and belongs in the result, so gaps up to
    `lesion_gap_mm` are bridged rather than cut.  A long one is something else: P92's
    LAD has a calcified spot at 26 mm followed by six millimetres with no lumen at all
    and a streak-artefact texture, and what lies beyond that is a 6 mm-wide structure
    that is not the LAD.  Measuring pericoronary fat around an axis nobody can verify is
    the error this whole tool exists to remove, so the far side of a gap that long is
    dropped and the loss is reported rather than papered over.
    """
    n = len(pts)
    inside = g.sample(pts) > thr
    none = dict(proximal_mm=0.0, distal_mm=0.0, kept_mm=float(arc[-1]),
                interior_gap_mm=0.0)
    if inside.all() or not inside.any():
        return pts, arc, none

    ds = arc[-1] / max(n - 1, 1)
    # a closing with a structuring element `k` samples long bridges gaps of up to k-1
    # samples, so k is the gap we are willing to call a lesion plus one
    bridge = max(int(round(p.lesion_gap_mm / ds)) + 1, 2)
    filled = ndi.binary_closing(inside, np.ones(bridge))
    lab, nlab = ndi.label(filled)
    if nlab == 0:
        return pts, arc, none
    sizes = ndi.sum(filled, lab, range(1, nlab + 1))
    keep = int(np.argmax(sizes)) + 1
    idx = np.nonzero(lab == keep)[0]
    a, b = int(idx[0]), int(idx[-1])
    if arc[b] - arc[a] < p.min_usable_mm:
        return pts, arc, dict(none, note="less than %.0f mm of the path was in contrast"
                              % p.min_usable_mm)

    cut = dict(proximal_mm=float(arc[a]), distal_mm=float(arc[-1] - arc[b]),
               kept_mm=float(arc[b] - arc[a]),
               interior_gap_mm=float(ds * (~inside[a:b + 1]).sum()))
    s = np.linspace(arc[a], arc[b], n)
    out = np.column_stack([np.interp(s, arc, pts[:, d]) for d in range(3)])
    return out, s - s[0], cut


def frames(pts):
    """A rotation-minimising orthonormal frame (double reflection, Wang et al. 2008).

    A Frenet frame spins wherever the curve is straight and flips through inflections,
    which would make the angular axis of any (arclength, angle) map meaningless.  This
    one carries the reference direction along the vessel with no twist, so angle zero
    means the same side of the vessel from one end to the other.
    """
    t = np.gradient(pts, axis=0)
    t /= np.linalg.norm(t, axis=1, keepdims=True)

    seed = np.array([0.0, 1.0, 0.0])
    if abs(seed @ t[0]) > 0.9:
        seed = np.array([0.0, 0.0, 1.0])
    u = np.cross(t[0], seed)

    U = np.zeros_like(pts)
    U[0] = u / np.linalg.norm(u)
    for i in range(len(pts) - 1):
        v1 = pts[i + 1] - pts[i]
        c1 = float(v1 @ v1)
        if c1 < 1e-12:
            U[i + 1] = U[i]
            continue
        uL = U[i] - (2.0 / c1) * (v1 @ U[i]) * v1
        tL = t[i] - (2.0 / c1) * (v1 @ t[i]) * v1
        v2 = t[i + 1] - tL
        c2 = float(v2 @ v2)
        u2 = uL if c2 < 1e-12 else uL - (2.0 / c2) * (v2 @ uL) * v2
        u2 -= (u2 @ t[i + 1]) * t[i + 1]
        U[i + 1] = u2 / max(np.linalg.norm(u2), 1e-12)

    V = np.cross(t, U)
    return t, U, V / np.linalg.norm(V, axis=1, keepdims=True)


# ---------------------------------------------------------------------------
# 2. a lumen boundary that cannot leak
# ---------------------------------------------------------------------------


def polar_stack(g, pts, U, V, p, r_max, dr):
    """I(section, angle, radius), sampled in each cross-section plane."""
    th = np.arange(p.n_theta) * 2 * np.pi / p.n_theta
    rr = np.arange(0.0, r_max + 1e-9, dr)
    c, s = np.cos(th)[:, None], np.sin(th)[:, None]
    off = rr[None, :, None] * (c[:, :, None] * U[:, None, None, :]
                               + s[:, :, None] * V[:, None, None, :])
    return g.sample(pts[:, None, None, :] + off), th, rr


def radial_step_limit(rr, n_theta, slope_deg):
    """How many radial bins the contour may move between two adjacent rays.

    The arc between neighbouring rays is 2*pi*r/A, so a fixed bin limit would mean a
    near-vertical wall near the axis and a gentle one far out.  Converting a wall
    slope into bins keeps the constraint the same shape everywhere.
    """
    dr = rr[1] - rr[0]
    arc = 2 * np.pi * np.maximum(rr, dr) / n_theta
    return np.maximum(1, np.ceil(np.tan(np.radians(slope_deg)) * arc / dr)).astype(int)


def _windowed_min(acc, step_limit, sizes):
    """min over r' with |r - r'| <= step_limit[r], along the last axis.

    A variable-width window, done as one fixed-width minimum filter per distinct width
    and then a select.  There are only a handful of distinct widths, and each filter is
    a single C pass, which is what keeps the whole dynamic program in the tens of
    milliseconds instead of the tens of seconds a Python gather would cost.
    """
    out = np.empty_like(acc)
    for d in sizes:
        m = ndi.minimum_filter1d(acc, size=2 * int(d) + 1, axis=-1, mode="nearest")
        sel = step_limit == d
        out[..., sel] = m[..., sel]
    return out


def cyclic_dp(cost, step_limit):
    """Exact minimum-cost closed contours through cost[section, angle, radius].

    The contour has to close, and a dynamic program needs an acyclic chain, so the open
    chain is run from every possible starting radius and only the runs that come back to
    within one step of their own start are allowed to compete.  That is what makes the
    contour join up cleanly instead of showing a seam at angle zero.

    Two passes.  The first carries all R starts at once and only needs the running cost,
    so it can use the filter above; it says where each section's contour should begin.
    The second re-runs those single starts keeping backpointers, which is cheap because
    there is now one start per section rather than R.  Every section in the vessel goes
    through both passes together, as one array.

    `step_limit` is per radius, so the wall-slope cap holds at every radius rather than
    only near the axis.
    """
    Z, A, R = cost.shape
    INF = np.float32(1e12)
    step_limit = np.ascontiguousarray(np.broadcast_to(np.asarray(step_limit), (R,)))
    sizes = np.unique(step_limit)

    # pass 1: cost only, all starts
    acc = np.full((Z, R, R), INF, np.float32)               # acc[z, start, radius]
    d = np.arange(R)
    acc[:, d, d] = cost[:, 0, :]
    for a in range(1, A):
        acc = _windowed_min(acc, step_limit, sizes) + cost[:, a, :][:, None, :]
    close_ok = np.abs(d[:, None] - d[None, :]) <= step_limit[:, None]
    flat = np.where(close_ok[None], acc, INF).reshape(Z, -1)
    best = np.argmin(flat, axis=1)
    s0 = (best // R).astype(np.int32)

    # pass 2: the chosen start per section, with backpointers
    offs = np.arange(-int(step_limit.max()), int(step_limit.max()) + 1)
    idx = np.clip(d[:, None] + offs[None, :], 0, R - 1)
    valid = ((d[:, None] + offs[None, :] >= 0) & (d[:, None] + offs[None, :] < R)
             & (np.abs(offs)[None, :] <= step_limit[:, None]))
    acc2 = np.full((Z, R), INF, np.float32)
    acc2[np.arange(Z), s0] = cost[np.arange(Z), 0, s0]
    back = np.zeros((A, Z, R), np.int16)
    for a in range(1, A):
        cand = np.where(valid[None], acc2[:, idx], INF)      # (z, r, offset)
        j = np.argmin(cand, axis=2)
        acc2 = np.take_along_axis(cand, j[..., None], 2)[..., 0] + cost[:, a, :]
        back[a] = np.clip(d[None, :] + offs[j], 0, R - 1)

    close2 = np.where(np.abs(s0[:, None] - d[None, :]) <= step_limit[s0][:, None],
                      acc2, INF)
    out = np.empty((Z, A), np.int32)
    out[:, A - 1] = np.argmin(close2, axis=1)
    zz = np.arange(Z)
    for a in range(A - 1, 0, -1):
        out[:, a - 1] = back[a, zz, out[:, a]]
    return out


def lumen_cost(I, rr, thr, p):
    """Cost of putting the boundary at each (angle, radius) of one section.

    Three terms, all in HU so the weights mean something:
      - the step: bright just inside minus dark just outside, negated so a strong
        edge is cheap.  This is the only term that locates the wall.
      - inside evidence: the mean HU shortfall of the interior.  Without it a contour
        could sit happily at radius zero, where there is no edge to disagree with.
      - calcium: a ray may not pass through it, because the lumen ends there.
    """
    w = max(int(round(p.edge_win_mm / p.dr_mm)), 1)
    R = I.shape[-1]
    box = ndi.uniform_filter1d(I, w, axis=-1, mode="nearest")
    half = w // 2 + 1
    inner = np.roll(box, half, axis=-1)
    outer = np.roll(box, -half, axis=-1)
    inner[..., :half] = box[..., :1]
    outer[..., -half:] = box[..., -1:]

    n = np.arange(1, R + 1)
    deficit = np.cumsum(np.maximum(thr - I, 0.0), axis=-1) / n

    cost = -(inner - outer) + p.w_inside * deficit
    cost[..., 0] += 1e3                                    # never collapse to a point
    cost[np.cumsum(I > p.calc_hu, axis=-1) > 0] += 1e3     # stop at calcium
    return np.ascontiguousarray(cost, np.float32)


def segment_lumen(g, pts, U, V, thr, p):
    I, th, rr = polar_stack(g, pts, U, V, p, p.r_max_mm, p.dr_mm)
    n_s = len(pts)
    base = lumen_cost(I, rr, thr, p)

    step_bins = radial_step_limit(rr, p.n_theta, p.wall_slope_deg)
    r_idx = cyclic_dp(base, step_bins)

    # second pass: pull each section towards its neighbours along the vessel, which is
    # what stops one bad section from disagreeing with the two beside it.
    prior = ndi.gaussian_filter1d(r_idx.astype(np.float32), 1.5, axis=0, mode="nearest")
    lim = p.s_step_mm / p.dr_mm
    grid = np.arange(base.shape[2])[None, None, :]
    pen = np.clip(np.abs(grid - prior[:, :, None]) - lim, 0, None) ** 2
    r_idx = cyclic_dp(base + p.w_smooth_s * pen, step_bins)

    r_mm = rr[r_idx]
    # honest per-section confidence, replacing the `detected` flag that used to report
    # success on sections where four rays in five had clamped at the ceiling.
    hu_axis = g.sample(pts)
    conf = dict(axis_hu=hu_axis.tolist(),
                axis_in_lumen=[bool(x) for x in hu_axis > thr],
                frac_at_ceiling=(r_mm >= p.r_max_mm - 1e-6).mean(1).tolist())
    return r_mm, th, conf


# ---------------------------------------------------------------------------
# 3. the outer vessel wall
# ---------------------------------------------------------------------------


def wall_cost(I, rr, r_lumen, p):
    """Cost of putting the outer wall at each (angle, radius).

    The adventitia has no contrast edge to find - that is why the commercial packages
    use a different algorithm here than for the lumen.  What is findable is where the
    tissue stops being tissue and starts being fat, and for a shell defined by an HU
    window that crossing is the operational definition of the boundary anyway.  Three
    terms, all in HU:

      - be at the crossing: |I - fat_hi| is small exactly where tissue becomes fat.
      - what is outside should be fat, so brighter-than-fat outside is penalised.
      - what is inside should not be, so fat inside is penalised too - otherwise the
        boundary slides outward through the depot to wherever the profile is flattest.
    """
    lo, hi = p.fat_gate
    w = max(int(round(p.edge_win_mm / p.wall_dr_mm)), 1)
    box = ndi.uniform_filter1d(I, w, axis=-1, mode="nearest")
    half = w // 2 + 1
    inner = np.roll(box, half, axis=-1)
    outer = np.roll(box, -half, axis=-1)
    inner[..., :half] = box[..., :1]
    outer[..., -half:] = box[..., -1:]

    cost = (np.abs(I - hi)
            + p.w_wall_out * np.maximum(outer - hi, 0.0)
            + p.w_wall_in * np.maximum(hi - inner, 0.0))
    # the wall cannot be inside the lumen, and calcium is wall, not fat
    floor = r_lumen + p.wall_min_gap_mm
    cost = cost + np.where(rr[None, None, :] < floor[:, :, None], 1e3, 0.0)
    return np.ascontiguousarray(cost, np.float32)


def segment_outer_wall(g, pts, U, V, r_lumen, p):
    I, th, rr = polar_stack(g, pts, U, V, p, p.wall_max_mm, p.wall_dr_mm)
    base = wall_cost(I, rr, r_lumen, p)
    step = radial_step_limit(rr, p.n_theta, p.wall_slope_deg)
    idx = cyclic_dp(base, step)
    prior = ndi.gaussian_filter1d(idx.astype(np.float32), 1.5, axis=0, mode="nearest")
    lim = p.s_step_mm / p.wall_dr_mm
    grid = np.arange(base.shape[2])[None, None, :]
    pen = np.clip(np.abs(grid - prior[:, :, None]) - lim, 0, None) ** 2
    idx = cyclic_dp(base + p.w_smooth_s * pen, step)
    r_wall = rr[idx]

    # honest per-ray confidence: a ray that never sees fat has no wall to find, and the
    # boundary it returns is the least-bad guess rather than a measurement
    lo, hi = p.fat_gate
    reach = np.zeros(r_wall.shape, bool)
    ai = np.arange(r_wall.shape[1])
    for s in range(r_wall.shape[0]):
        j = np.clip((r_wall[s] / p.wall_dr_mm).astype(int), 0, I.shape[2] - 2)
        ahead = np.minimum(j + int(round(1.0 / p.wall_dr_mm)), I.shape[2] - 1)
        reach[s] = I[s, ai, ahead] <= hi
    return r_wall, reach


# ---------------------------------------------------------------------------
# 4. fat
# ---------------------------------------------------------------------------


def fat_peak(hist, centres):
    """Where the fat peak actually sits, and how wide it is, from the whole histogram.

    The windowed mean is a mean over [-190, -30] of a distribution that does not end at
    those numbers, so a shift of the peak changes both which voxels are kept and the
    average of the ones that are.  The mode and a moment fit around it do not care where
    the window is, which is what makes them comparable across keV and kernel.
    """
    sm = ndi.gaussian_filter1d(hist.astype(float), 2.0)
    band = (centres >= -200) & (centres <= -20)
    nan = dict(mode_hu=float("nan"), fit_mean_hu=float("nan"),
               fit_sigma_hu=float("nan"))
    if not band.any() or sm[band].max() <= 0:
        return nan
    i = int(np.flatnonzero(band)[np.argmax(sm[band])])
    mu, sg = float(centres[i]), 28.0

    # A Gaussian is a parabola in log space, so the fit is one weighted quadratic per
    # iteration.  The window is [mu - 3s, mu + 1.2s]: clipped on the right because a
    # shell that reaches into myocardium has a tissue shoulder there, but NOT left-only.
    # Fitting the left flank alone looks safer and measures worse - it turns the width
    # into the only thing setting the peak position, so it inherits the kernel's effect
    # on noise, which is exactly the nuisance this statistic exists to escape.
    for _ in range(4):
        sel = ((centres > mu - 3.0 * sg) & (centres < mu + 1.2 * sg)
               & (hist > 0.05 * hist.max()))
        if sel.sum() < 5:
            break
        x, y = centres[sel], hist[sel].astype(float)
        try:
            a, b, _c = np.polyfit(x, np.log(y), 2, w=np.sqrt(y))
        except (np.linalg.LinAlgError, ValueError):
            break
        if a >= 0:
            break
        mu, sg = float(-b / (2 * a)), float(np.sqrt(-1.0 / (2 * a)))
        # A shell that is mostly myocardium has a tissue shoulder taller than the fat
        # peak, and the fit will walk into it given the chance.  It is not allowed to
        # leave the peak it was seeded on.
        if (not (np.isfinite(mu) and np.isfinite(sg))
                or abs(mu - centres[i]) > 30.0 or mu > -35.0):
            return dict(mode_hu=float(centres[i]), fit_mean_hu=float(centres[i]),
                        fit_sigma_hu=float("nan"), fit_failed=True)
        sg = float(np.clip(sg, 5.0, 80.0))
    return dict(mode_hu=float(centres[i]), fit_mean_hu=mu, fit_sigma_hu=sg)


def truncated_mean(mu, sigma, lo, hi):
    """Mean of a Gaussian(mu, sigma) restricted to [lo, hi].

    This is what the published FAI measures, so evaluating it at a FIXED sigma turns the
    familiar number into one that no longer moves when the reconstruction noise does.
    d(FAI)/d(sigma) runs -0.05 to -0.35 HU per HU of width, so a ten-HU noise difference
    between two arms manufactures a two-to-four HU FAI difference out of nothing - the
    same size as the effect being looked for.
    """
    if not (np.isfinite(mu) and np.isfinite(sigma)) or sigma <= 0:
        return float("nan")
    a, b = (lo - mu) / sigma, (hi - mu) / sigma
    pdf = lambda x: np.exp(-0.5 * x * x) / np.sqrt(2 * np.pi)
    cdf = lambda x: 0.5 * (1 + np.math.erf(x / np.sqrt(2))) if hasattr(np, "math") \
        else 0.5 * (1 + erf(x / np.sqrt(2)))
    z = cdf(b) - cdf(a)
    if z <= 1e-12:
        return float("nan")
    return float(mu + sigma * (pdf(a) - pdf(b)) / z)


def shell_stats(I, rr, r_in, thickness, p, dvol, with_maps=False):
    """Volume-weighted attenuation statistics over one shell, plus its histogram.

    A polar grid is uniform in (theta, r) but the tissue is not - an annulus at 5 mm
    carries twice the volume of one at 2.5 mm - so every average here is weighted by r.
    Counting samples instead tilts the result towards the inner shell, and pericoronary
    fat has a real radial attenuation gradient.
    """
    lo, hi = p.fat_gate
    sel = ((rr[None, None, :] >= r_in[:, :, None])
           & (rr[None, None, :] <= (r_in + thickness)[:, :, None]))
    hu = I[sel]
    wt = np.broadcast_to(rr[None, None, :], I.shape)[sel]
    fat = (hu >= lo) & (hu <= hi)

    edges = np.arange(p.hist_lo, p.hist_hi + p.hist_bin, p.hist_bin)
    hist, _ = np.histogram(hu, bins=edges, weights=wt * dvol)
    centres = 0.5 * (edges[:-1] + edges[1:])

    out = dict(shell_volume_mm3=float(wt.sum() * dvol),
               fat_volume_mm3=float(wt[fat].sum() * dvol),
               fat_fraction=float(wt[fat].sum() / wt.sum()) if wt.sum() else float("nan"),
               hist_lo=p.hist_lo, hist_bin=p.hist_bin,
               hist_mm3=[round(float(x), 4) for x in hist])
    out.update(fat_peak(hist, centres))
    tot = hist.sum()
    if tot > 0:
        out.update(frac_below_window=float(hist[centres < lo].sum() / tot),
                   frac_above_window=float(hist[centres > hi].sum() / tot),
                   frac_in_window=float(hist[(centres >= lo) & (centres <= hi)].sum()
                                        / tot))
    # the published statistic, re-evaluated at a fixed width, so it stops tracking the
    # reconstruction noise.  Reported next to the raw one, never instead of it.
    out["hu_mean_sigma%d" % int(p.ref_sigma_hu)] = (
        float("nan") if out.get("fit_failed")
        else truncated_mean(out["fit_mean_hu"], p.ref_sigma_hu, lo, hi))

    if fat.any():
        fw, fh = wt[fat], hu[fat]
        mean = float(np.average(fh, weights=fw))
        order = np.argsort(fh)
        cw = np.cumsum(fw[order])
        out.update(hu_mean=mean,
                   hu_std=float(np.sqrt(np.average((fh - mean) ** 2, weights=fw))),
                   hu_median=float(fh[order][np.searchsorted(cw, 0.5 * cw[-1])]))
    else:
        out.update(hu_mean=float("nan"), hu_std=float("nan"),
                   hu_median=float("nan"))

    if with_maps:
        band = np.where(sel & (I >= lo) & (I <= hi), I, 0.0)
        bw = np.where(sel & (I >= lo) & (I <= hi), rr[None, None, :], 0.0)
        with np.errstate(invalid="ignore", divide="ignore"):
            out["fai_map"] = (band * bw).sum(2) / np.where(bw.sum(2) > 0, bw.sum(2),
                                                           np.nan)
            out["fai_per_section"] = ((band * bw).sum((1, 2))
                                      / np.where(bw.sum((1, 2)) > 0, bw.sum((1, 2)),
                                                 np.nan))
    return out


def decompose(hu70, hu150, kev=(70, 150)):
    """Split a shell mean into how much lipid is there and how much iodine is in it.

    Siemens VMI is rank-2, so two energies exhaust the spectral information and a third
    adds only correlated noise.  The solve is linear, so applying it to the shell mean
    and averaging the voxels first are algebraically the same thing - and doing it on
    the mean avoids amplifying per-voxel noise through the inverse.

    Basis (HU relative to water, from the gecatsim/NIST tables): adipose -104.09 at
    70 keV and -78.85 at 150 keV; soft tissue +43.07 and +40.62; iodine +25.861 and
    +4.458 HU per mg/mL.  With volume closure f_adipose + f_soft = 1.
    """
    if tuple(kev) != (70, 150):
        return dict(adipose_fraction=float("nan"), iodine_mg_ml=float("nan"),
                    basis_note="the basis is tabulated at 70 and 150 keV; this pair is "
                               "%s/%s keV, so only the raw difference is reported" % kev)
    if not (np.isfinite(hu70) and np.isfinite(hu150)):
        return dict(adipose_fraction=float("nan"), iodine_mg_ml=float("nan"))
    a, b = hu70 - 43.07, hu150 - 40.62
    return dict(adipose_fraction=float(0.0018320 * a - 0.0106269 * b),
                iodine_mg_ml=float(0.0490921 * a - 0.0604690 * b))


def measure_fat(g, pts, U, V, r_lumen, r_wall, wall_ok, th, p,
                gap_mm=1.0, ring_mm=3.0, max_radius_mm=20.0, break_mm=1.0, g2=None):
    """Every shell worth reporting, measured on one polar sampling of the vessel.

    The primary is the published convention: adipose tissue from the OUTER VESSEL WALL
    outwards, over a radial distance equal to the vessel diameter, gated to -190..-30 HU
    (Antonopoulos 2017, Oikonomou 2018, Kotanidis & Antoniades 2021).  The app's own
    lumen-keyed 1 mm gap + fixed 3 mm ring is kept beside it so nothing already measured
    changes meaning without being visible, and the ambiguity in "the vessel diameter" -
    lumen or outer wall - is resolved by reporting both.
    """
    lo, hi = p.fat_gate
    rr = np.arange(0.0, max_radius_mm + 1e-9, p.dr_mm)
    c, s = np.cos(th)[:, None], np.sin(th)[:, None]
    off = rr[None, :, None] * (c[:, :, None] * U[:, None, None, :]
                               + s[:, :, None] * V[:, None, None, :])
    I = g.sample(pts[:, None, None, :] + off)              # (S, A, R)

    ds = float(np.mean(np.diff(
        np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(pts, axis=0), axis=1))])))
    dvol = p.dr_mm * (2 * np.pi / len(th)) * ds

    d_lumen = float(2 * np.mean(np.sqrt((r_lumen ** 2).mean(1))))
    d_outer = float(2 * np.mean(np.sqrt((r_wall ** 2).mean(1))))
    thickness = d_lumen if p.shell_scale == "lumen_diameter" else d_outer

    primary = shell_stats(I, rr, r_wall + p.shell_gap_mm, thickness, p, dvol,
                          with_maps=True)
    alt = shell_stats(I, rr, r_wall + p.shell_gap_mm,
                      d_outer if p.shell_scale == "lumen_diameter" else d_lumen,
                      p, dvol)
    gapped = shell_stats(I, rr, r_wall + p.shell_gap_alt_mm, thickness, p, dvol)
    legacy = shell_stats(I, rr, r_lumen + gap_mm, ring_mm, p, dvol)

    # depot: walk out from the wall, tolerate a break_mm non-fat gap, stop at the first
    # longer one.  Unchanged, and reported separately - this is an epicardial-fat
    # measure, not the FAI shell, and the two are read independently.
    fat = (I >= lo) & (I <= hi)
    S, A, R = fat.shape
    gap_bins = int(round(break_mm / p.dr_mm))
    i0 = np.clip(np.round(r_wall / p.dr_mm).astype(int), 0, R - 1)
    r_out = np.zeros((S, A))
    open_at_max = np.zeros((S, A), bool)
    for a_ in range(S):
        row = fat[a_]
        for b_ in range(A):
            i, run, last = i0[a_, b_], 0, i0[a_, b_]
            while i < R:
                if row[b_, i]:
                    run, last = 0, i
                else:
                    run += 1
                    if run > gap_bins:
                        break
                i += 1
            r_out[a_, b_] = rr[min(last, R - 1)]
            open_at_max[a_, b_] = i >= R

    spectral = None
    if g2 is not None:
        # the same voxels at both energies: select on the 70 keV fat gate, then average
        # each energy over that one set, so the two means describe the same tissue
        sel = ((rr[None, None, :] >= (r_wall + p.shell_gap_mm)[:, :, None])
               & (rr[None, None, :] <= (r_wall + p.shell_gap_mm + thickness)[:, :, None]))
        pos = pts[:, None, None, :] + off
        J = g2.sample(pos)
        w = np.broadcast_to(rr[None, None, :], I.shape)[sel]
        a70, a150 = I[sel], J[sel]
        gate = (a70 >= lo) & (a70 <= hi)
        if gate.any():
            m70 = float(np.average(a70[gate], weights=w[gate]))
            m150 = float(np.average(a150[gate], weights=w[gate]))
            spectral = dict(kev=list(p.kev_pair), hu_low=m70, hu_high=m150,
                            hu_low_minus_high=m70 - m150,
                            **decompose(m70, m150, p.kev_pair))
            spectral["coverage"] = g2.covers(pts)

    return dict(
        gate_hu=[lo, hi], spectral=spectral,
        shell=dict(primary, scale=p.shell_scale, thickness_mm=thickness,
                   starts_at="outer wall"),
        shell_alt=dict(alt, scale=("outer_diameter" if p.shell_scale == "lumen_diameter"
                                   else "lumen_diameter"),
                       thickness_mm=(d_outer if p.shell_scale == "lumen_diameter"
                                     else d_lumen), starts_at="outer wall"),
        shell_gapped=dict(gapped, scale=p.shell_scale, thickness_mm=thickness,
                          gap_mm=p.shell_gap_alt_mm,
                          starts_at="%.2f mm outside the outer wall"
                                    % p.shell_gap_alt_mm),
        legacy_ring=dict(legacy, gap_mm=gap_mm, ring_mm=ring_mm, starts_at="lumen"),
        mean_lumen_diameter_mm=d_lumen, mean_outer_diameter_mm=d_outer,
        wall_found_fraction=float(np.mean(wall_ok)),
        r_wall_mm=r_wall,
        r_out_mm=r_out, thickness_mm=np.maximum(r_out - r_wall, 0.0),
        open_at_max=open_at_max, max_radius_mm=max_radius_mm, break_mm=break_mm)


# ---------------------------------------------------------------------------
# view stack
# ---------------------------------------------------------------------------


def write_mask(g, img, pts, U, V, r_lumen, th, fat, p, out_path):
    """Rasterise the result back into the image grid, so it can be dropped on the
    original series in any viewer.

    1 = the fat inside the FAI shell, which is what the reported number averages
    2 = the lumen
    3 = fat in the depot beyond the shell, which is read separately
    """
    r_out = fat["r_out_mm"]
    r_wall = fat["r_wall_mm"]
    start = r_wall + p.shell_gap_mm
    stop = start + fat["shell"]["thickness_mm"]
    reach = float(r_out.max()) + 1.0

    tree = cKDTree(pts)
    # only voxels that could possibly be in reach, in grid index space
    lo = g.to_vox(pts.min(0) - reach)[0]
    hi = g.to_vox(pts.max(0) + reach)[0]
    sl = tuple(slice(max(int(np.floor(min(a, b))), 0),
                     min(int(np.ceil(max(a, b))) + 1, n))
               for a, b, n in zip(lo, hi, g.vol.shape))
    zz, yy, xx = np.mgrid[sl]
    world = g.to_world(np.column_stack([zz.ravel(), yy.ravel(), xx.ravel()]))
    d, s_i = tree.query(world, distance_upper_bound=reach)
    ok = np.isfinite(d)

    out = np.zeros(g.vol.shape, np.uint8)
    if ok.any():
        w = world[ok]
        si = s_i[ok]
        vec = w - pts[si]
        # Keep only what lies inside the slab of its own section, so the end caps are
        # flat instead of hemispherical - a Voronoi assignment alone would sweep in
        # everything past the last sample and inflate the volume.
        tang = np.gradient(pts, axis=0)
        tang /= np.linalg.norm(tang, axis=1, keepdims=True)
        ds = float(np.median(np.linalg.norm(np.diff(pts, axis=0), axis=1)))
        axial = np.abs((vec * tang[si]).sum(1))
        rad = np.linalg.norm(vec, axis=1)
        ang = np.arctan2((vec * V[si]).sum(1), (vec * U[si]).sum(1)) % (2 * np.pi)
        ai = np.clip((ang / (2 * np.pi) * len(th)).astype(int), 0, len(th) - 1)
        hu = g.vol[zz.ravel()[ok], yy.ravel()[ok], xx.ravel()[ok]]
        gate = (hu >= fat["gate_hu"][0]) & (hu <= fat["gate_hu"][1])
        inslab = axial <= 0.5 * ds + g.iso
        in_shell = (rad >= start[si, ai]) & (rad <= stop[si, ai]) & gate & inslab
        in_depot = (rad > stop[si, ai]) & (rad <= r_out[si, ai]) & gate & inslab
        in_lumen = (rad <= r_lumen[si, ai]) & inslab
        lab = np.where(in_shell, 1, np.where(in_lumen, 2,
                                             np.where(in_depot, 3, 0))).astype(np.uint8)
        out[zz.ravel()[ok], yy.ravel()[ok], xx.ravel()[ok]] = lab

    m = sitk.GetImageFromArray(out)
    m.SetSpacing((g.iso,) * 3)
    m.SetOrigin(tuple(g.origin))
    m.SetDirection(img.GetDirection())
    sitk.WriteImage(m, out_path, True)
    return dict(shell=int((out == 1).sum()), lumen=int((out == 2).sum()),
                depot=int((out == 3).sum()))


def write_cpr(g, pts, arc, U, V, r_lumen, out_path, n_ang=8, half_mm=9.0,
              step_mm=0.2):
    """A stretched curved-planar reformat: the vessel laid out flat, one image per
    rotation angle.  This is the view a reader can actually judge - a mis-placed axis
    or a leaking boundary shows up here at a glance, where twelve separate discs do
    not - and it is the picture the commercial packages put on screen first."""
    cols = np.linspace(0, arc[-1], max(int(round(arc[-1] / step_mm)), 2))
    across = np.arange(-half_mm, half_mm + 1e-9, step_mm)
    cen = np.column_stack([np.interp(cols, arc, pts[:, d]) for d in range(3)])
    stack = np.empty((n_ang, len(across), len(cols)), np.int16)
    edges = np.empty((n_ang, 2, len(cols)), np.float32)
    A = r_lumen.shape[1]
    for k in range(n_ang):
        phi = k * np.pi / n_ang           # a plane, so half a turn covers every view
        e = (np.cos(phi) * np.array([np.interp(cols, arc, U[:, d]) for d in range(3)])
             + np.sin(phi) * np.array([np.interp(cols, arc, V[:, d]) for d in range(3)])).T
        e /= np.linalg.norm(e, axis=1, keepdims=True)
        pos = cen[None, :, :] + across[:, None, None] * e[None, :, :]
        stack[k] = np.clip(np.round(g.sample(pos)), -32768, 32767).astype(np.int16)
        a0 = int(round(phi / (2 * np.pi) * A)) % A
        a1 = (a0 + A // 2) % A
        edges[k, 0] = np.interp(cols, arc, r_lumen[:, a0])
        edges[k, 1] = np.interp(cols, arc, r_lumen[:, a1])
    stack.tofile(out_path)
    return dict(angles=n_ang, half_width_mm=half_mm, step_mm=step_mm,
                width=len(cols), height=len(across), length_mm=float(arc[-1]),
                lumen_edge_mm=edges.tolist())


def write_views(g, pts, U, V, out_path, n_sections=12, pixels=128, half_mm=10.0):
    idx = np.unique(np.linspace(0, len(pts) - 1, n_sections).round().astype(int))
    u = np.linspace(-half_mm, half_mm, pixels)
    stack = np.empty((len(idx), pixels, pixels), np.int16)
    for n, s in enumerate(idx):
        gxy = (pts[s][None, None, :] + u[None, :, None] * U[s][None, None, :]
               + u[:, None, None] * V[s][None, None, :])
        stack[n] = np.clip(np.round(g.sample(gxy)), -32768, 32767).astype(np.int16)
    stack.tofile(out_path)
    return idx, pixels, half_mm


# ---------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------


def _curvature(pts):
    out = {}
    for sig in (0, 2, 5, 10):
        q = ndi.gaussian_filter1d(pts, sig, axis=0, mode="nearest") if sig else pts
        d1 = np.gradient(q, axis=0)
        d2 = np.gradient(d1, axis=0)
        den = np.linalg.norm(d1, axis=1) ** 3
        num = np.linalg.norm(np.cross(d1, d2), axis=1)
        out[str(sig)] = np.where(den > 1e-9, num / np.maximum(den, 1e-9), 0.0).tolist()
    return out


def process_vessel(root, patient, vessel, img, p, write=True, img2=None,
                   valid_z2=None, series=None, series2=None):
    base = os.path.join(root, patient)
    cpath = os.path.join(base, "centerline", vessel + ".json")
    if not os.path.exists(cpath):
        return None
    # Always re-centre the operator's original path, never a path this tool already
    # produced: re-running would otherwise smooth and trim the same vessel again and
    # again, and the segment would creep shorter with every pass.
    src = cpath.replace(".json", ".orig.json")
    with open(src if os.path.exists(src) else cpath) as fh:
        C = json.load(fh)
    pts0 = np.array(C["samples_mm"], float)
    n = len(pts0)

    t = time.time()
    g = Grid(img, pts0, p)
    g2 = None
    if img2 is not None:
        g2 = Grid(img2, pts0, p)
        g2.valid_z = valid_z2
    t_grid = time.time() - t

    t = time.time()

    def attempt(thr_override):
        path, _m, thr, info = recentre(g, pts0, p, thr_override=thr_override)
        pts, arc = resample_centerline(path, n, p, ends=info["anchors"])
        off = float(np.mean(g.sample(pts) < thr))
        return off, path, thr, info, pts, arc

    _off, path, thr, info, pts, arc = attempt(None)

    # Second pass on a threshold read off the vessel rather than off its neighbourhood.
    # The first pass has to guess the lumen level from a high percentile of everything
    # inside the tube, which on P92's LAD included part of a heart chamber and pushed the
    # threshold to 292 HU - high enough to erode the distal vessel out of the medialness
    # field entirely.  Once there is an axis, the lumen level is simply the attenuation
    # along it, and the half-maximum against the local peri-vascular tissue follows.
    axis_hu = g.sample(pts)
    inside = axis_hu > thr
    if inside.sum() >= 8:
        lumen2 = float(np.median(axis_hu[inside]))
        shell = shell_hu(g, pts, 3.0, 6.5)
        bg2 = float(np.median(shell)) if shell.size else info["bg_hu"]
        thr2 = 0.5 * (lumen2 + bg2)
        if abs(thr2 - thr) > 15.0:
            alt = attempt(thr2)
            # accept only if it does not make things worse and the result is still a
            # coronary: a lower threshold fuses the artery to whatever it touches, and
            # the medial axis of a fused blob is not the axis of the vessel
            wide = ndi.distance_transform_edt(
                ndi.binary_opening(ndi.gaussian_filter(g.vol, 0.5 / g.iso) > alt[2],
                                   np.ones((3, 3, 3))), sampling=g.iso)
            med_on_path = wide[tuple(np.clip(np.round(g.to_vox(alt[4])).astype(int), 0,
                                             np.array(g.vol.shape) - 1).T)]
            if (alt[0] <= _off
                    and np.median(med_on_path) * 2 <= p.max_lumen_diameter_mm):
                _off, path, thr, info, pts, arc = alt
                info["threshold_refined_to"] = thr2

    # If the operator's segment runs off the end of the opacified vessel there is
    # nothing to re-centre onto out there, and dragging the path into fat to reach
    # their last click would quietly poison the FAI.  Give up those millimetres
    # instead, and say how many were given up.
    pts, arc, trimmed = trim_to_lumen(g, pts, arc, thr, p)
    ds0 = float(np.median(np.linalg.norm(np.diff(pts0, axis=0), axis=1)))
    trimmed["anchor_proximal_mm"] = info["anchor_skipped"][0] * ds0
    trimmed["anchor_distal_mm"] = info["anchor_skipped"][1] * ds0
    tangent, U, V = frames(pts)
    t_centre = time.time() - t

    t = time.time()
    r_lumen, th, conf = segment_lumen(g, pts, U, V, thr, p)
    t_lumen = time.time() - t

    t = time.time()
    r_wall, wall_ok = segment_outer_wall(g, pts, U, V, r_lumen, p)
    t_wall = time.time() - t

    t = time.time()
    fat = measure_fat(g, pts, U, V, r_lumen, r_wall, wall_ok, th, p, g2=g2)
    t_fat = time.time() - t

    d_moved = cKDTree(pts).query(pts0)[0]
    hu_before, hu_after = g.sample(pts0), g.sample(pts)
    r_eq = np.sqrt((r_lumen ** 2).mean(1))

    qc = dict(
        vessel=vessel, n_sections=int(n), length_mm=float(arc[-1]),
        measured_on=series, second_energy_series=series2,
        operator_length_mm=float(np.linalg.norm(np.diff(pts0, axis=0), axis=1).sum()),
        trimmed_mm=trimmed, tube_mm=info["tube_mm"],
        threshold_hu=thr, lumen_hu=info["lumen_hu"], background_hu=info["bg_hu"],
        recentre_mm=dict(median=float(np.median(d_moved)),
                         p90=float(np.percentile(d_moved, 90)),
                         max=float(d_moved.max())),
        axis_hu_before=dict(median=float(np.median(hu_before)),
                            min=float(hu_before.min()),
                            frac_outside_lumen=float(np.mean(hu_before < thr))),
        axis_hu_after=dict(median=float(np.median(hu_after)),
                           min=float(hu_after.min()),
                           frac_outside_lumen=float(np.mean(hu_after < thr))),
        lumen_diameter_mm=dict(median=float(2 * np.median(r_eq)),
                               min=float(2 * r_eq.min()), max=float(2 * r_eq.max())),
        frac_rays_at_ceiling=float(np.mean(r_lumen >= p.r_max_mm - 1e-6)),
        outer_diameter_mm=fat["mean_outer_diameter_mm"],
        wall_found_fraction=fat["wall_found_fraction"],
        shell_thickness_mm=fat["shell"]["thickness_mm"],
        fai_hu=fat["shell"]["hu_mean"],
        fai_legacy_hu=fat["legacy_ring"]["hu_mean"],
        fat_peak_hu=fat["shell"]["fit_mean_hu"],
        fat_peak_sigma_hu=fat["shell"]["fit_sigma_hu"],
        spectral=fat["spectral"],
        seconds=dict(grid=t_grid, recentre=t_centre, lumen=t_lumen, wall=t_wall,
                     fat=t_fat))

    if not write:
        return qc

    def backup(pth):
        b = pth.replace(".json", ".orig.json")
        if os.path.exists(pth) and not os.path.exists(b):
            os.replace(pth, b)

    backup(cpath)
    with open(cpath, "w") as fh:
        json.dump(dict(
            vessel=vessel, seeds_mm=C.get("seeds_mm", []),
            ostium_fraction=C.get("ostium_fraction"),
            samples_mm=pts.tolist(), arclength_mm=arc.tolist(),
            tangent=tangent.tolist(), normal=U.tolist(), binormal=V.tolist(),
            curvature_per_mm=_curvature(pts),
            trusted=[bool(x) for x in hu_after > thr],
            recentred_from_mm=pts0.tolist(), recentre_offset_mm=d_moved.tolist(),
            trimmed_mm=trimmed,
            method="medialness minimum-cost path in a %.0f mm tube around the "
                   "operator's path; rotation-minimising frame" % p.tube_mm), fh)

    lpath = os.path.join(base, "lumen", vessel + ".json")
    backup(lpath)
    with open(lpath, "w") as fh:
        json.dump(dict(
            vessel=vessel, theta_rad=th.tolist(), r_theta_mm=r_lumen.tolist(),
            r_eq_mm=r_eq.tolist(), detected=conf["axis_in_lumen"],
            axis_hu=conf["axis_hu"], frac_at_ceiling=conf["frac_at_ceiling"],
            max_search_radius_mm=p.r_max_mm,
            method="cyclic dynamic programming over the polar section, wall slope "
                   "capped at %.0f deg, coupled along the vessel"
                   % p.wall_slope_deg), fh)

    fpath = os.path.join(base, "fat", vessel + ".json")
    backup(fpath)
    shell = {k: v for k, v in fat["shell"].items()
             if k not in ("fai_map", "fai_per_section")}
    shell["fai_per_section_hu"] = [None if not np.isfinite(x) else float(x)
                                   for x in fat["shell"]["fai_per_section"]]
    out = dict(
        vessel=vessel, segment_achieved_mm=float(arc[-1]), gate_hu=fat["gate_hu"],
        measured_on=series, second_energy_series=series2,
        mean_lumen_diameter_mm=fat["mean_lumen_diameter_mm"],
        mean_outer_diameter_mm=fat["mean_outer_diameter_mm"],
        wall_found_fraction=fat["wall_found_fraction"],
        r_wall_mm=fat["r_wall_mm"].tolist(),
        shell=shell, shell_alt=fat["shell_alt"], shell_gapped=fat["shell_gapped"],
        legacy_ring=fat["legacy_ring"], spectral=fat["spectral"],
        # the published cut-off, kept only so the viewer can draw where it would fall.
        # It is a 120 kVp energy-integrating number and does not transfer to a virtual
        # monoenergetic reconstruction on a photon-counting scanner.
        reference_cutoff_hu=-70.1,
        reference_cutoff_note="Oikonomou 2018 / CRISP-CT, 120 kVp EID. Not validated "
                              "at this keV, kernel or detector; do not classify on it.",
        method="shell from the outer wall outwards over a radial distance equal to the "
               "mean %s diameter (Antonopoulos 2017, Oikonomou 2018, Kotanidis 2021)"
               % fat["shell"]["scale"].replace("_", " "))
    ds = arc[-1] / max(n - 1, 1)
    dth = 2 * np.pi / len(th)
    vol = float((fat["thickness_mm"] * (fat["r_out_mm"] + fat["r_wall_mm"]) / 2
                 * dth).sum() * ds)
    # kept separate on purpose: this is an epicardial-fat depot measure, read on its own
    # terms, not the FAI shell above
    out["depot"] = dict(
        r_out_mm=fat["r_out_mm"].tolist(), thickness_mm=fat["thickness_mm"].tolist(),
        mean_thickness_mm=fat["thickness_mm"].mean(1).tolist(),
        open_at_max=fat["open_at_max"].tolist(), gate_hu=fat["gate_hu"],
        starts_at="outer wall", max_radius_mm=fat["max_radius_mm"],
        break_mm=fat["break_mm"], volume_mm3=vol)

    mpath = os.path.join(base, "fat", vessel + "_mask.nii.gz")
    if os.path.exists(mpath) and not os.path.exists(mpath.replace(".nii.gz",
                                                                 ".orig.nii.gz")):
        os.replace(mpath, mpath.replace(".nii.gz", ".orig.nii.gz"))
    out["mask"] = dict(
        file="fat/" + vessel + "_mask.nii.gz", spacing_mm=p.iso_mm,
        labels={"1": "fat inside the FAI shell", "2": "lumen",
                "3": "fat in the depot beyond the shell"},
        n_voxels=write_mask(g, img, pts, U, V, r_lumen, th, fat, p, mpath))
    out["depot"]["mask_file"] = "fat/" + vessel + "_mask.nii.gz"
    out["depot"]["mask_spacing_mm"] = p.iso_mm
    with open(fpath, "w") as fh:
        json.dump(out, fh)

    vdir = os.path.join(base, "views")
    os.makedirs(vdir, exist_ok=True)
    idx, px, hm = write_views(g, pts, U, V,
                              os.path.join(vdir, vessel + "_sections.i16"))
    cpr = write_cpr(g, pts, arc, U, V, r_lumen,
                    os.path.join(vdir, vessel + "_cpr.i16"))

    # the (arclength, angle) FAI map, as int16 HU with a sentinel where the ring held
    # no fat.  It is 58 kB this way and 1.5 MB as JSON.
    fm = fat["shell"]["fai_map"]
    np.where(np.isfinite(fm), np.clip(fm, -1000, 1000), -32768
             ).astype(np.int16).tofile(os.path.join(vdir, vessel + "_faimap.i16"))

    vpath = os.path.join(vdir, vessel + ".json")
    backup(vpath)
    with open(vpath, "w") as fh:
        json.dump(dict(vessel=vessel, pixels=px, half_width_mm=hm,
                       section_index=idx.tolist(),
                       arc_mm=[float(arc[i]) for i in idx],
                       arclength_mm=arc.tolist(),
                       data_file="views/" + vessel + "_sections.i16",
                       cpr_file="views/" + vessel + "_cpr.i16", cpr=cpr,
                       faimap_file="views/" + vessel + "_faimap.i16",
                       faimap=dict(sections=int(fm.shape[0]), angles=int(fm.shape[1]),
                                   nodata=-32768),
                       axis_note="column -> +normal, row -> +binormal; centre pixel "
                                 "is the centerline"), fh)
    return qc


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=r"Z:\Shu Nie\shu_pcat")
    ap.add_argument("--images", default="Z:\\",
                    help="the folder that stands in for /Volumes/Molloilab")
    ap.add_argument("--patient", action="append", default=[])
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--vessels", nargs="*", default=list(VESSELS))
    ap.add_argument("--cache", default=None, help="where the NIfTI cache lives")
    ap.add_argument("--measure-on", default=None, metavar="keV",
                    help="measure every patient on this VMI energy instead of whatever "
                         "the manifest recorded. P92 was exported at 70 keV and P106 at "
                         "65 keV, a systematic offset of about 2 HU before any biology.")
    ap.add_argument("--second-energy", default=None, metavar="keV",
                    help="a second VMI energy of the same acquisition, e.g. 150. Its "
                         "shell mean is decomposed with the primary into an adipose "
                         "fraction and an iodine concentration.")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args(argv)

    pats = a.patient
    if a.all or not pats:
        pats = sorted(d for d in os.listdir(a.root)
                      if os.path.exists(os.path.join(a.root, d, "manifest.json")))
    cache = a.cache or os.path.join(a.root, ".cache")
    p = Params()
    if a.measure_on and a.second_energy:
        p.kev_pair = (int(a.measure_on), int(a.second_energy))

    allqc = {}
    for pat in pats:
        with open(os.path.join(a.root, pat, "manifest.json")) as fh:
            man = json.load(fh)
        t = time.time()
        sdir = resolve_series_dir(man, a.images)
        if a.measure_on:
            alt = find_energy_series(os.path.dirname(sdir), a.measure_on)
            if alt is None:
                print("  %s has no %s keV series - measuring on the manifest's choice"
                      % (pat, a.measure_on))
            else:
                sdir = alt
        img = load_volume(sdir, cache)
        print("\n%s  %s  %s  loaded in %.1fs"
              % (pat, os.path.basename(sdir), img.GetSize(), time.time() - t))

        img2, valid_z2, d2 = None, None, None
        if a.second_energy:
            d2 = find_energy_series(os.path.dirname(sdir), a.second_energy)
            if d2 is None:
                print("  no %s keV series beside this one - spectral pair skipped"
                      % a.second_energy)
            else:
                z, blocks = series_z(d2)
                if len(blocks) > 1:
                    print("  %s keV holds %d slices in %d blocks - that folder is "
                          "missing slices, so coverage is checked per vessel"
                          % (a.second_energy, len(z), len(blocks)))
                img2, valid_z2 = load_volume(d2, cache), blocks

        allqc[pat] = {}
        for v in a.vessels:
            qc = process_vessel(a.root, pat, v, img, p, write=not a.dry_run,
                                img2=img2, valid_z2=valid_z2,
                                series=os.path.basename(sdir),
                                series2=os.path.basename(d2) if img2 is not None
                                else None)
            if qc is None:
                continue
            allqc[pat][v] = qc
            s = qc["seconds"]
            print("  %-4s len %5.1fmm  moved med %4.2f max %4.2fmm | axis HU %6.0f -> "
                  "%6.0f (%3.0f%% -> %3.0f%% off-lumen) | lumen d %4.2fmm  ceiling "
                  "%4.1f%%  FAI %6.1f (peak %6.1f) | %4.1fs"
                  % (v, qc["length_mm"], qc["recentre_mm"]["median"],
                     qc["recentre_mm"]["max"], qc["axis_hu_before"]["median"],
                     qc["axis_hu_after"]["median"],
                     100 * qc["axis_hu_before"]["frac_outside_lumen"],
                     100 * qc["axis_hu_after"]["frac_outside_lumen"],
                     qc["lumen_diameter_mm"]["median"],
                     100 * qc["frac_rays_at_ceiling"], qc["fai_hu"],
                     qc["fat_peak_hu"], sum(s.values())))
        if not a.dry_run:
            with open(os.path.join(a.root, pat, "qc.json"), "w") as fh:
                json.dump(allqc[pat], fh, indent=1)
    return allqc


if __name__ == "__main__":
    main()
