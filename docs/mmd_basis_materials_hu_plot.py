"""
Hounsfield number (HU) of candidate MMD basis materials at 70 keV vs 150 keV,
computed from NIST XCOM mass attenuation coefficients.

HU(E) = 1000 * (mu_mat(E) - mu_water(E)) / mu_water(E)
      = 1000 * (rho_mat * (mu/rho)_mat(E) / (rho_water * (mu/rho)_water(E)) - 1)

Water is (0, 0) by definition. Pure iodine is off-scale for clinical imaging
(~+155000 HU at 70 keV), so an additional point is shown for iodine at
10 mg/mL in water (representative coronary-lumen opacification in CCTA).

A Compton-only reference line has slope 1 (same relative HU at both energies);
deviation toward the x-axis indicates excess photoelectric absorption at
70 keV (high-Z signature).
"""
from __future__ import annotations

import itertools

import matplotlib.pyplot as plt
import numpy as np

# Mass attenuation coefficients (cm^2/g) from NIST XCOM.
MU_RHO = {
    "Air":            {"mu70": 0.1759, "mu150": 0.1356, "rho": 0.001205},
    "Water":          {"mu70": 0.1929, "mu150": 0.1505, "rho": 1.000},
    "Lipid":          {"mu70": 0.1861, "mu150": 0.1460, "rho": 0.950},
    "Collagen":       {"mu70": 0.1875, "mu150": 0.1483, "rho": 1.350},
    "Hydroxyapatite": {"mu70": 0.4481, "mu150": 0.1818, "rho": 3.156},
    "Iodine (pure)":  {"mu70": 6.097,  "mu150": 0.7961, "rho": 4.933},
}

# Linear attenuation of water (cm^-1) — the HU reference.
MU_W_70 = MU_RHO["Water"]["mu70"] * MU_RHO["Water"]["rho"]
MU_W_150 = MU_RHO["Water"]["mu150"] * MU_RHO["Water"]["rho"]

def hu_at(mu_mat_70: float, mu_mat_150: float) -> tuple[float, float]:
    """Convert linear attenuation (cm^-1) at 70/150 keV to HU at each energy."""
    return (
        1000.0 * (mu_mat_70 - MU_W_70) / MU_W_70,
        1000.0 * (mu_mat_150 - MU_W_150) / MU_W_150,
    )

# Build points for each pure material.
points = {}
for name, m in MU_RHO.items():
    mu70 = m["mu70"] * m["rho"]
    mu150 = m["mu150"] * m["rho"]
    points[name] = hu_at(mu70, mu150)

# Clinical iodine concentration: 10 mg/mL in water (typical CCTA lumen).
c_iodine_g_per_cm3 = 0.010
mu70_10mgmL  = c_iodine_g_per_cm3 * MU_RHO["Iodine (pure)"]["mu70"]  + (1.0 - c_iodine_g_per_cm3) * MU_W_70
mu150_10mgmL = c_iodine_g_per_cm3 * MU_RHO["Iodine (pure)"]["mu150"] + (1.0 - c_iodine_g_per_cm3) * MU_W_150
points["Iodine (10 mg/mL)"] = hu_at(mu70_10mgmL, mu150_10mgmL)

# Style per material.
STYLE = {
    "Air":                {"color": "#17becf", "marker": "o"},
    "Water":              {"color": "#1f77b4", "marker": "o"},
    "Lipid":              {"color": "#ff7f0e", "marker": "o"},
    "Collagen":           {"color": "#2ca02c", "marker": "o"},
    "Hydroxyapatite":     {"color": "#d62728", "marker": "o"},
    "Iodine (pure)":      {"color": "#9467bd", "marker": "o"},
    "Iodine (10 mg/mL)":  {"color": "#c49cf0", "marker": "o"},
}

# Connecting polyline order.
ORDER = ["Air", "Water", "Lipid", "Collagen", "Hydroxyapatite", "Iodine (10 mg/mL)", "Iodine (pure)"]

# ---------------------------------------------------------------------------
# Two-panel plot: full range (symlog) + clinical zoom
# ---------------------------------------------------------------------------

fig, (axL, axR) = plt.subplots(1, 2, figsize=(14.5, 6.8), dpi=140)

def draw(ax, ylim_lo=None, ylim_hi=None, xlim_lo=None, xlim_hi=None, symlog=True):
    xs = [points[n][0] for n in ORDER]
    ys = [points[n][1] for n in ORDER]
    # Every pairwise edge. Triangles formed by any 3 points become visible;
    # large/open triangles = good MMD basis, near-collinear = degenerate.
    for a, b in itertools.combinations(ORDER, 2):
        xa, ya = points[a]
        xb, yb = points[b]
        ax.plot([xa, xb], [ya, yb], color="#999999", lw=0.9, ls="--",
                alpha=0.7, zorder=1)
    for name in ORDER:
        x, y = points[name]
        s = STYLE[name]
        ax.scatter(x, y, s=180,
                   marker=s["marker"], color=s["color"],
                   edgecolor="black", linewidth=0.8, zorder=3)
        # Offset label so it doesn't overlap the marker.
        dx = 0.12 if x >= 0 else -0.12
        ax.annotate(
            name, xy=(x, y),
            xytext=(8, 8), textcoords="offset points",
            fontsize=10, color=s["color"], fontweight="bold",
        )

    # Compton-only slope = 1 reference line through origin.
    lim_hi = xlim_hi if xlim_hi is not None else max(xs) * 1.1
    lim_lo = xlim_lo if xlim_lo is not None else min(xs) * 1.1
    ref = np.linspace(lim_lo, lim_hi, 500)
    ax.plot(ref, ref, color="#bbbbbb", lw=0.8, zorder=0, label="Compton-only (slope = 1)")

    ax.axhline(0, color="#dddddd", lw=0.6, zorder=0)
    ax.axvline(0, color="#dddddd", lw=0.6, zorder=0)

    if symlog:
        ax.set_xscale("symlog", linthresh=100)
        ax.set_yscale("symlog", linthresh=100)
    ax.set_xlabel("HU at 70 keV", fontsize=11)
    ax.set_ylabel("HU at 150 keV", fontsize=11)
    ax.grid(which="both", alpha=0.3)
    if ylim_lo is not None:
        ax.set_ylim(ylim_lo, ylim_hi)
    if xlim_lo is not None:
        ax.set_xlim(xlim_lo, xlim_hi)
    ax.legend(loc="lower right", fontsize=9, framealpha=0.9)

# Left: full symlog range (shows pure iodine and air).
draw(axL, symlog=True, xlim_lo=-1500, xlim_hi=300_000, ylim_lo=-1500, ylim_hi=300_000)
axL.set_title("HU(70) vs HU(150) — full range (symlog)", fontsize=12)

# Right: linear zoom covering air (-1000 HU) to heavily calcified (-> HAp).
draw(axR, symlog=False, xlim_lo=-1200, xlim_hi=7500, ylim_lo=-1200, ylim_hi=3500)
axR.set_title("HU(70) vs HU(150) — clinical range (linear)", fontsize=12)

fig.suptitle("Hounsfield response of MMD basis candidates\n(NIST XCOM, water reference)",
             fontsize=13, y=1.02)
plt.tight_layout()

out_path = "/Users/shunie/Developer/PCAT/pcat-workstation-v2/docs/mmd_basis_materials_hu.png"
plt.savefig(out_path, dpi=180, bbox_inches="tight")
print(f"wrote {out_path}")

# Print table.
print()
print(f"{'Material':<22}  {'HU @ 70 keV':>14}  {'HU @ 150 keV':>14}  {'slope 150/70':>14}")
print("-" * 72)
for name in ORDER:
    h70, h150 = points[name]
    slope = (h150 / h70) if abs(h70) > 1e-6 else float("nan")
    print(f"{name:<22}  {h70:>14.1f}  {h150:>14.1f}  {slope:>14.3f}")
