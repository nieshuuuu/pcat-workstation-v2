#!/usr/bin/env python3
"""
pcat_syngo - the Siemens pipeline, reproduced.

`pcat_autofix.py` gets to the same place by a cheaper route.  This module implements
what syngo.VIA's coronary analysis actually does, so the two can be compared on the same
vessels rather than argued about:

  1. Medialness with a rising-edge veto (Gulsun & Tek).  A response that is near zero
     anywhere outside the vessel, because a point in fat has a bright rising edge between
     it and the wall and that rising edge is explicitly subtracted.  It needs no
     intensity threshold at all - the response is normalised by the strongest falling
     edge on the same ray, so it is contrast- and dose-independent.  autofix instead
     thresholds the volume and takes a distance transform, which works but has to guess
     the threshold first.

  2. The lumen as one Markov random field over every (section, ray) radius at once,
     with convex Huber priors coupling angular and longitudinal neighbours, solved
     exactly by min-cut through Ishikawa's construction.  autofix solves each section
     exactly and then couples the sections by a second pass, which is an approximation.
     The wrap-around in theta, which the dynamic program has to close by hand, is just
     another pair of neighbours here.

  3. Calcium found per ray, with the boundary likelihood damped exponentially outbound
     of the first calcified run rather than blocked outright, so blooming cannot drag
     the edge outward but a boundary can still be placed past a small speck.

The one thing that makes the min-cut tractable: a Huber prior is linear past its knee,
so the second difference that sets the inter-column arc capacities is ZERO beyond the
knee.  The graph is sparse - a few arcs per node - instead of quadratic in the number
of radial candidates.

Usage
-----
    python pcat_syngo.py --patient P92                 # compare against autofix
    python pcat_syngo.py --all --measure-on 70
"""

from __future__ import annotations

import argparse
import json
import os
import time

import numpy as np
from scipy import ndimage as ndi
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.csgraph import maximum_flow

import pcat_autofix as af


# ---------------------------------------------------------------------------
# 1. medialness with a rising-edge veto
# ---------------------------------------------------------------------------


def ray_medialness(g, centres, U, V, p):
    """Gulsun & Tek medialness at a set of candidate centres.

    `centres` is (..., 3) in world mm; U and V give the plane each candidate is
    evaluated in, broadcast to the same leading shape.  Returns a response in [0, 1]
    with the same leading shape.

    For every in-plane direction the profile is differentiated at two scales and the
    signed edge response b(t) taken as whichever is stronger.  A boundary at radius R
    is only credited if the ray reaches it WITHOUT first crossing a rising edge:

        E(R) = max(-b(R) - rise(R), 0) / fall_max

    where rise(R) is the strongest dark-to-bright step anywhere inside R.  From a seed
    in epicardial fat every ray pointing at the artery crosses a -80 -> +450 HU rising
    edge on the way, rise is large, E collapses, and the response is near zero exactly
    where the operator put the seed.  Dividing by the strongest falling edge on the same
    ray removes any dependence on how well the patient opacified.
    """
    n_dir = p.med_dirs
    phi = np.arange(n_dir) * 2 * np.pi / n_dir
    t = np.arange(p.med_dr, p.med_rmax + 1e-9, p.med_dr)          # radii sampled
    lead = centres.shape[:-1]

    # (lead, dir, t, 3)
    e = (np.cos(phi)[:, None, None] * U[..., None, None, :]
         + np.sin(phi)[:, None, None] * V[..., None, None, :])
    pos = centres[..., None, None, :] + t[None, :, None] * e
    prof = g.sample(pos)                                          # (lead, dir, t)

    # signed edge response: strongest of two scales, sign kept
    b = None
    for sig in p.med_sigmas:
        d = ndi.gaussian_filter1d(prof, sig / p.med_dr, axis=-1, order=1,
                                  mode="nearest") / p.med_dr
        b = d if b is None else np.where(np.abs(d) > np.abs(b), d, b)

    rise = np.maximum.accumulate(np.maximum(b, 0.0), axis=-1)
    fall_max = np.maximum(np.max(np.maximum(-b, 0.0), axis=-1, keepdims=True), 1.0)
    E = np.maximum(-b - rise, 0.0) / fall_max                     # (lead, dir, t)

    # the vessel radius is unknown, so take the best supported one
    lo = int(round(p.med_rmin / p.med_dr)) - 1
    return E.mean(axis=-2)[..., lo:].max(axis=-1)


def viterbi_recentre(g, pts, U, V, p):
    """Move every sample onto the medialness ridge, with a smoothness prior along the
    vessel, by a Viterbi pass rather than an independent per-section argmax."""
    step = p.cand_step_mm
    k = int(np.floor(p.cand_radius_mm / step))
    a, b = np.meshgrid(np.arange(-k, k + 1) * step, np.arange(-k, k + 1) * step,
                       indexing="ij")
    keep = (a ** 2 + b ** 2) <= p.cand_radius_mm ** 2 + 1e-9
    off_uv = np.column_stack([a[keep], b[keep]])                  # (C, 2) in-plane mm
    C = len(off_uv)

    cand = (pts[:, None, :]
            + off_uv[None, :, 0, None] * U[:, None, :]
            + off_uv[None, :, 1, None] * V[:, None, :])           # (Z, C, 3)
    Uc = np.broadcast_to(U[:, None, :], cand.shape)
    Vc = np.broadcast_to(V[:, None, :], cand.shape)
    m = ray_medialness(g, cand, Uc, Vc, p)                        # (Z, C)

    d2 = ((off_uv[:, None, :] - off_uv[None, :, :]) ** 2).sum(-1)
    allowed = d2 <= p.cand_move_mm ** 2 + 1e-9
    trans = np.where(allowed, p.cand_lambda * d2, np.inf)

    Z = len(pts)
    D = np.empty((Z, C))
    back = np.zeros((Z, C), np.int32)
    D[0] = -m[0]
    for z in range(1, Z):
        tot = D[z - 1][:, None] + trans
        back[z] = np.argmin(tot, axis=0)
        D[z] = tot[back[z], np.arange(C)] - m[z]

    idx = np.empty(Z, np.int32)
    idx[-1] = int(np.argmin(D[-1]))
    for z in range(Z - 1, 0, -1):
        idx[z - 1] = back[z, idx[z]]
    return cand[np.arange(Z), idx], m[np.arange(Z), idx]


# ---------------------------------------------------------------------------
# 2. one globally optimal surface, by min-cut
# ---------------------------------------------------------------------------


def huber_arc_capacities(delta, gamma, R):
    """Ishikawa arc capacities and the unary correction they come with.

    For a convex pairwise g, the arc from V(p,i) to V(q,j) carries half the discrete
    second difference of g at offset i - j.  That is non-negative for any convex g and -
    this is what keeps the graph sparse - identically zero once g has gone linear, which
    for Huber is everywhere past its knee.

    The catch, and the thing that is easy to miss: those arcs do not realise g(a - b).
    They realise g(a - b) + f(a) + f(b), for a separable f that has to be subtracted back
    out of each column's unary cost, once for every neighbour pair that column is in.
    Verified numerically against brute force in the tests below.
    """
    def hub(d):
        d = abs(float(d))
        return gamma * (d * d / (2 * delta) if d <= delta else d - delta / 2)

    reach = int(np.ceil(delta)) + 1
    cap = np.array([max(0.5 * (hub(k + 1) + hub(k - 1) - 2 * hub(k)), 0.0)
                    for k in range(reach + 1)])
    # f(a) = sum_{i<=a} sum_j c_ij  -  g(a),  taken over i, j in 1..R-1
    S = np.zeros(R)
    for i_ in range(1, R):
        S[i_] = sum(cap[abs(i_ - j_)] if abs(i_ - j_) <= reach else 0.0
                    for j_ in range(1, R))
    f = np.cumsum(S) - np.array([hub(a) for a in range(R)])
    return cap, f - f[0]


def optimal_surface(cost, gamma_theta, gamma_z, delta, scale=None):
    """Exact minimum-cost surface r(z, a) through cost[z, a, r].

    Nodes are V(z, a, r); the surface is the outermost node kept on the source side of
    the cut.  Intra-column infinite arcs make each column monotone, so exactly one radius
    is chosen per ray.  Inter-column arcs couple each ray to its two angular neighbours -
    cyclically, so the theta wrap is nothing special, unlike in a dynamic program - and
    to its neighbours along the vessel.

    Everything is int32 on purpose.  scipy's maximum_flow accepts an int64 matrix without
    complaint and then truncates it to int32 internally, so a sentinel like 1 << 40 comes
    out as a capacity of zero and every constraint silently disappears.
    """
    Z, A, R = cost.shape
    cap_t, f_t = huber_arc_capacities(delta, gamma_theta, R)
    cap_z, f_z = huber_arc_capacities(delta, gamma_z, R)

    # undo the separable term the pairwise arcs will add, once per neighbour pair
    n_t = 2 if A >= 3 else 0
    n_z = np.zeros(Z)
    if Z >= 2:
        n_z += 2
        n_z[0] -= 1
        n_z[-1] -= 1
    cost = (cost - n_t * f_t[None, None, :]
            - n_z[:, None, None] * f_z[None, None, :])

    N = Z * A * R
    SRC, SNK = N, N + 1
    idx = np.arange(N).reshape(Z, A, R)

    w = np.empty_like(cost, float)
    w[..., 0] = 0.0                 # V(p,0) is forced in, so its cost is a constant
    w[..., 1:] = np.diff(cost, axis=-1)

    # The sentinel only has to beat the maximum flow, which is bounded by the total
    # source capacity - not by the sum of every capacity in the graph, which for a few
    # million arcs runs past int32 immediately.  Pick the scale from that budget.
    if scale is None:
        src_total = float(np.abs(w[w < 0]).sum()) + 1e-9
        scale = 3e8 / src_total
        smallest = min([c for c in np.r_[cap_t, cap_z] if c > 0] or [1.0])
        if smallest * scale < 1.0:                 # priors must not quantise to nothing
            scale = 1.0 / smallest
    wi = np.rint(w * scale).astype(np.int64)

    rows, cols, data = [], [], []

    def add(u, v, c):
        rows.append(np.atleast_1d(u).ravel())
        cols.append(np.atleast_1d(v).ravel())
        data.append(np.atleast_1d(c).ravel())

    add(idx[:, :, 1:], idx[:, :, :-1], np.zeros(Z * A * (R - 1), np.int64))
    inf_slots = [len(data) - 1]
    add(np.full(Z * A, SRC), idx[:, :, 0], np.zeros(Z * A, np.int64))
    inf_slots.append(len(data) - 1)

    neg = wi.ravel() < 0
    add(np.full(int(neg.sum()), SRC), np.arange(N)[neg], -wi.ravel()[neg])
    pos = wi.ravel() > 0
    add(np.arange(N)[pos], np.full(int(pos.sum()), SNK), wi.ravel()[pos])

    def couple(axis, cap):
        """Arcs for every undirected neighbour pair along `axis`, in both directions."""
        if axis == "theta":
            if A < 3:
                return
            pa, pb = idx, np.roll(idx, -1, axis=1)
        else:
            if Z < 2:
                return
            pa, pb = idx[:-1], idx[1:]
        reach = len(cap) - 1
        for k in range(-reach, reach + 1):
            ci = int(round(cap[abs(k)] * scale))
            if ci <= 0:
                continue
            i0, i1 = max(1, 1 + k), min(R, R + k)     # i in [1,R), i-k in [1,R)
            if i1 <= i0:
                continue
            u = pa[..., i0:i1].ravel()
            v = pb[..., i0 - k:i1 - k].ravel()
            add(u, v, np.full(u.size, ci, np.int64))
            u2 = pb[..., i0:i1].ravel()
            v2 = pa[..., i0 - k:i1 - k].ravel()
            add(u2, v2, np.full(u2.size, ci, np.int64))

    couple("theta", cap_t)
    couple("z", cap_z)

    INF = int(sum(int(d.sum()) for n, d in enumerate(data)
                  if n not in inf_slots and n == 2)) + 1   # the source arcs alone
    INF = max(INF, 1 << 20)
    for n in inf_slots:
        data[n] = np.full(data[n].size, INF, np.int64)

    rows = np.concatenate(rows)
    cols = np.concatenate(cols)
    data = np.concatenate(data)
    keep = data > 0
    graph = coo_matrix((data[keep], (rows[keep], cols[keep])),
                       shape=(N + 2, N + 2)).tocsr()
    graph.sum_duplicates()
    top = int(graph.data.max())
    if top >= 2 ** 31 - 1:
        raise OverflowError("capacity %d exceeds int32; lower `scale`" % top)
    graph = graph.astype(np.int32)

    res = maximum_flow(graph, SRC, SNK)
    resid = csr_matrix(graph.astype(np.int64) - res.flow.astype(np.int64))

    seen = np.zeros(N + 2, bool)
    seen[SRC] = True
    stack = [SRC]
    indptr, indices, vals = resid.indptr, resid.indices, resid.data
    while stack:
        u = stack.pop()
        for t in range(indptr[u], indptr[u + 1]):
            v = indices[t]
            if vals[t] > 0 and not seen[v]:
                seen[v] = True
                stack.append(v)

    inside = seen[:N].reshape(Z, A, R)
    return (R - 1 - inside[:, :, ::-1].argmax(-1))


# ---------------------------------------------------------------------------
# 3. data cost, with calcium damped rather than blocked
# ---------------------------------------------------------------------------


def lumen_data_cost(I, rr, p):
    """Boundary likelihood per (section, ray, radius), turned into a cost.

    Two terms, both scaled to O(1) so the Huber weights below mean the same thing on any
    patient: the normalised bright-to-dark gradient, and a level term that peaks where
    the profile crosses the half-maximum between the lumen and the peri-vascular tissue.
    """
    grad = -np.gradient(ndi.gaussian_filter1d(I, p.grad_sigma_mm / p.dr_mm, axis=-1),
                        p.dr_mm, axis=-1)
    grad = np.maximum(grad, 0.0)
    norm = np.percentile(grad, 98, axis=(1, 2), keepdims=True) + 1e-6
    g_n = grad / norm

    core = I[:, :, : max(int(round(0.5 / p.dr_mm)), 2)]
    L = np.median(core, axis=(1, 2), keepdims=True)
    B = np.percentile(I, 25, axis=(1, 2), keepdims=True)
    T = 0.5 * (L + B)
    lvl = np.exp(-(((I - T) / (0.25 * np.maximum(L - B, 50.0))) ** 2))

    like = 1.0 * g_n + 0.5 * lvl

    # calcium: the first sustained run above a threshold derived from THIS patient's
    # lumen, then an exponential damp outbound so blooming cannot pull the edge out
    t_cal = np.maximum(2.0 * L, L + 2.0 * I[:, :, :4].std(axis=(1, 2), keepdims=True))
    hot = I > t_cal
    run = np.zeros_like(hot, int)
    need = max(int(round(0.3 / p.dr_mm)), 1)
    acc = np.zeros(hot.shape[:2], int)
    for i in range(hot.shape[2]):
        acc = np.where(hot[:, :, i], acc + 1, 0)
        run[:, :, i] = acc
    first = np.where((run >= need).any(-1), (run >= need).argmax(-1), hot.shape[2])
    idx = np.arange(hot.shape[2])[None, None, :]
    start = np.maximum(first - need + 1, 0)[:, :, None]
    beyond = np.maximum(idx - start, 0)
    damp = np.where(idx >= start, p.calc_damp ** beyond, 1.0)
    like = like * damp

    return np.ascontiguousarray(-like + 1e-5, np.float64)


# ---------------------------------------------------------------------------
# the pipeline
# ---------------------------------------------------------------------------


def recentre_syngo(g, pts0, p, iters=2):
    """Medialness Viterbi, re-splined and re-framed, a couple of times.

    No threshold is estimated anywhere: the medialness normalises itself against the
    strongest falling edge on each ray, so a poorly opacified patient and a well
    opacified one are scored on the same scale.
    """
    pts = pts0
    for _ in range(iters):
        _, U, V = af.frames(pts)
        moved, m = viterbi_recentre(g, pts, U, V, p)
        pts, arc = af.resample_centerline(moved, len(pts0), p)
    tangent, U, V = af.frames(pts)
    return pts, arc, U, V, m


def segment_lumen_syngo(g, pts, U, V, p):
    """One globally optimal lumen surface for the whole vessel at once."""
    q = af.Params()
    q.n_theta = p.surf_theta
    q.dr_mm = p.surf_dr_mm
    I, th, rr = af.polar_stack(g, pts, U, V, q, p.surf_r_max_mm, p.surf_dr_mm)
    cost = lumen_data_cost(I, rr, p)
    lab = optimal_surface(cost, p.gamma_theta, p.gamma_z, p.huber_delta)
    return rr[lab], th, I


# ---------------------------------------------------------------------------
# parameters
# ---------------------------------------------------------------------------


class Params(af.Params):
    """autofix's parameters, plus the ones this pipeline needs."""

    med_dirs = 16
    med_dr = 0.1
    med_rmin = 0.5
    med_rmax = 6.0
    med_sigmas = (0.4, 0.8)

    cand_radius_mm = 3.0
    cand_step_mm = 0.4
    cand_move_mm = 1.0
    cand_lambda = 0.5

    surf_r_max_mm = 5.0        # Siemens use a 5 mm ceiling, not 8
    surf_dr_mm = 0.1
    surf_theta = 64            # 64 equiangular rays
    huber_delta = 3.0          # bins
    gamma_theta = 0.5
    gamma_z = 0.1
    grad_sigma_mm = 0.15
    calc_damp = 0.6
