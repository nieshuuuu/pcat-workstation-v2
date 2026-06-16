"""
Plot mass attenuation coefficients (mu/rho) at 70 keV vs 150 keV for candidate
MMD basis materials. Values from NIST XCOM (photon total, cm^2/g).

Compound values (water, lipid, collagen, hydroxyapatite) are computed as
mass-weighted sums of elemental mu/rho, using log-log interpolation of the
NIST-tabulated elemental values at bracketing energies (60/80 keV for 70 keV
and 100/200 keV for 150 keV).

Iodine is elemental.

Lipid composition: ICRP Adipose Tissue (H 11.4%, C 59.8%, N 0.7%, O 27.8%, rho 0.95).
Collagen composition: H 6.5%, C 33.0%, N 11.0%, O 50.0% (rho 1.35).
Hydroxyapatite: Ca10(PO4)6(OH)2 -> Ca 39.89%, P 18.50%, O 41.41%, H 0.20% (rho 3.156).

Axis: mu/rho (cm^2/g). Low-Z materials (water, lipid, collagen) cluster near
the diagonal; Compton dominates so they are spectrally indistinguishable.
HAp lifts off (Ca + P photoelectric at 70 keV). Iodine lifts off dramatically
(Z=53, K-edge at 33 keV). The off-diagonal displacement is MMD's separating
signal.
"""
from __future__ import annotations

import itertools

import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------------------
# NIST XCOM mass attenuation coefficients, cm^2/g
# (photon total, incoherent + coherent + photoelectric)
# ---------------------------------------------------------------------------

# Material -> (mu/rho at 70 keV, mu/rho at 150 keV)  in cm^2/g
#   Water         : NIST "Water, Liquid"          rho = 1.000 g/cm^3
#   Lipid (adip.) : NIST "Tissue, Adipose (ICRP)" rho = 0.950 g/cm^3
#   Collagen      : weighted from elemental NIST  rho = 1.350 g/cm^3
#   HAp           : weighted from elemental NIST  rho = 3.156 g/cm^3
#   Iodine        : NIST elemental Z=53            rho = 4.933 g/cm^3
MATERIALS = {
    "Air":             {"mu70": 0.1759, "mu150": 0.1356, "rho": 0.001205, "marker": "o", "color": "#17becf"},
    "Water":           {"mu70": 0.1929, "mu150": 0.1505, "rho": 1.000, "marker": "o", "color": "#1f77b4"},
    "Lipid":           {"mu70": 0.1861, "mu150": 0.1460, "rho": 0.950, "marker": "o", "color": "#ff7f0e"},
    "Collagen":        {"mu70": 0.1875, "mu150": 0.1483, "rho": 1.350, "marker": "o", "color": "#2ca02c"},
    "Hydroxyapatite":  {"mu70": 0.4481, "mu150": 0.1818, "rho": 3.156, "marker": "o", "color": "#d62728"},
    "Iodine":          {"mu70": 6.097,  "mu150": 0.7961, "rho": 4.933, "marker": "o", "color": "#9467bd"},
}

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

fig, ax = plt.subplots(figsize=(8.5, 7.0), dpi=140)

xs = [m["mu70"] for m in MATERIALS.values()]
ys = [m["mu150"] for m in MATERIALS.values()]

# Draw every pairwise edge so basis-material triangles become visible.
# Thick edges (water cluster <-> high-Z) are long; thin ones (within the
# water/lipid/collagen cluster) collapse to near-zero length -> the cluster
# is the degenerate triangle that cannot serve as an MMD basis.
names = list(MATERIALS.keys())
for a, b in itertools.combinations(names, 2):
    xa, ya = MATERIALS[a]["mu70"], MATERIALS[a]["mu150"]
    xb, yb = MATERIALS[b]["mu70"], MATERIALS[b]["mu150"]
    ax.plot([xa, xb], [ya, yb], color="#999999", lw=0.9, ls="--",
            alpha=0.7, zorder=1)

# Plot each point.
for name, m in MATERIALS.items():
    ax.scatter(
        m["mu70"], m["mu150"],
        s=180,
        marker=m["marker"],
        color=m["color"],
        edgecolor="black",
        linewidth=0.8,
        zorder=3,
        label=f"{name} (rho={m['rho']:.2f} g/cm^3)",
    )
    # Annotate.
    dx = 1.10
    dy = 1.00
    ax.annotate(
        name,
        xy=(m["mu70"], m["mu150"]),
        xytext=(m["mu70"] * dx, m["mu150"] * dy),
        fontsize=10,
        color=m["color"],
        fontweight="bold",
    )

# Compton-only diagonal reference (mu/rho roughly constant ratio for low-Z).
ref_x = np.array([0.15, 10.0])
ratio = MATERIALS["Water"]["mu150"] / MATERIALS["Water"]["mu70"]
ax.plot(ref_x, ref_x * ratio, color="#cccccc", lw=0.8, zorder=0,
        label=f"Compton-only slope (ratio = {ratio:.3f})")

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel(r"$\mu/\rho$ at 70 keV  (cm$^2$/g)", fontsize=11)
ax.set_ylabel(r"$\mu/\rho$ at 150 keV  (cm$^2$/g)", fontsize=11)
ax.set_title("NIST mass attenuation: MMD basis candidates\n70 keV vs 150 keV", fontsize=12)
ax.grid(which="both", alpha=0.3)
ax.legend(loc="lower right", fontsize=9, framealpha=0.9)

# Margins.
ax.set_xlim(0.12, 12.0)
ax.set_ylim(0.12, 1.2)

plt.tight_layout()

out_path = "/Users/shunie/Developer/PCAT/pcat-workstation-v2/docs/mmd_basis_materials.png"
plt.savefig(out_path, dpi=180, bbox_inches="tight")
print(f"wrote {out_path}")

# Also print the numeric table for the record.
print()
print(f"{'Material':<18}  {'mu/rho @ 70 keV':>16}  {'mu/rho @ 150 keV':>18}  {'ratio (150/70)':>16}")
print("-" * 78)
for name, m in MATERIALS.items():
    ratio_m = m["mu150"] / m["mu70"]
    print(f"{name:<18}  {m['mu70']:>16.4f}  {m['mu150']:>18.4f}  {ratio_m:>16.3f}")
