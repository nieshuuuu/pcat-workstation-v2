"""
How much volumetric composition change in pericoronary fat is detectable on CT?

Model PCAT as a 3-component volumetric mix of triglyceride (lipid), water,
and collagen-like ECM. Plot HU vs water volume fraction on a mono-energetic
70 keV reconstruction, with typical healthy / inflamed PCAT ranges marked
and noise bands for per-voxel and ring-averaged detection.

Physics: at 70 keV, soft-tissue mu/rho is nearly constant, so HU is linear
in volume fractions -> slope ~1.1 HU per percentage-point water replacing
lipid. The entire clinical FAI signal (healthy ~ -95 HU, inflamed ~ -70 HU)
is this single density axis.
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------------------
# NIST XCOM mu/rho (cm^2/g) at 70 keV, plus density (g/cm^3).
# Triglyceride is represented by triolein (C57H104O6, rho 0.915).
# Collagen is the same elemental weighting used in the basis plots.
# ---------------------------------------------------------------------------

COMP = {
    # mu/rho at 70 keV, mu/rho at 150 keV, rho.
    "Triglyceride": dict(mu70=0.1818, mu150=0.1424, rho=0.915),
    "Water":        dict(mu70=0.1929, mu150=0.1505, rho=1.000),
    "Collagen":     dict(mu70=0.1875, mu150=0.1483, rho=1.350),
}

def mu_lin(m, energy):
    return m[f"mu{energy}"] * m["rho"]

def hu(mu_mat, mu_water):
    return 1000.0 * (mu_mat - mu_water) / mu_water

MU_W_70 = mu_lin(COMP["Water"], 70)
MU_W_150 = mu_lin(COMP["Water"], 150)

# ---------------------------------------------------------------------------
# PCAT mixture model
#
# f_L + f_W + f_C = 1 by volume (ignoring microvasculature ~2%, cells ~5%).
# Sweep f_W from 0 to 40% (covers pure fat to heavily inflamed). For each
# f_W, split the remainder between lipid and collagen using a fixed ECM
# fraction typical of adipose tissue (4% baseline, rising with inflammation).
# ---------------------------------------------------------------------------

def ecm_fraction(fw):
    # Linear rise from 4% at 0% water to 8% at 40% water (histology).
    return 0.04 + 0.10 * fw

def pcat_hu_at_water_fraction(fw, energy):
    fc = ecm_fraction(fw)
    fl = max(0.0, 1.0 - fw - fc)
    mu = (
        fl * mu_lin(COMP["Triglyceride"], energy)
        + fw * mu_lin(COMP["Water"], energy)
        + fc * mu_lin(COMP["Collagen"], energy)
    )
    mu_w = MU_W_70 if energy == 70 else MU_W_150
    return hu(mu, mu_w)

fws = np.linspace(0.0, 0.40, 401)
hu_70 = np.array([pcat_hu_at_water_fraction(f, 70) for f in fws])
hu_150 = np.array([pcat_hu_at_water_fraction(f, 150) for f in fws])

# Sensitivity (HU per percentage-point water) — finite difference.
d_hu_per_pp = (hu_70[-1] - hu_70[0]) / (100 * (fws[-1] - fws[0]))

# ---------------------------------------------------------------------------
# Clinical anchors
# ---------------------------------------------------------------------------

# Antonopoulos 2017: mean FAI healthy -95, inflamed -70 (observational).
# Mapping to water fraction using this exact model.
def water_frac_for_hu(target_hu, energy=70):
    return float(np.interp(target_hu, hu_70 if energy == 70 else hu_150, fws))

anchors = {
    "Pure fat (no water)":       (0.00, "#1f6fb4"),
    "Healthy PCAT (-95 HU)":     (water_frac_for_hu(-95), "#2ca02c"),
    "Inflamed PCAT (-70 HU)":    (water_frac_for_hu(-70), "#ff7f0e"),
    "Severe inflammation (-50)": (water_frac_for_hu(-50), "#d62728"),
}

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

fig, ax = plt.subplots(figsize=(10.5, 6.2), dpi=140)

# Main curve.
ax.plot(100 * fws, hu_70, color="#1f77b4", lw=2.0,
        label=f"70 keV VMI  (slope {d_hu_per_pp:.2f} HU / %pp water)")
ax.plot(100 * fws, hu_150, color="#9467bd", lw=1.5, ls="--",
        label="150 keV VMI (flatter, less Compton sensitivity)")

# Noise bands.
# Per-voxel CT noise: ~10 HU std in 0.5x0.5x0.5 mm voxel, typical CCTA.
# Ring-averaged FAI: ~2 HU std (thousands of voxels in the 3mm ring).
# Shade bands at +/- 1 std around each anchor HU to show detectability.
for name, (fw, color) in anchors.items():
    h = pcat_hu_at_water_fraction(fw, 70)
    ax.scatter(100 * fw, h, s=90, color=color, edgecolor="black",
               linewidth=0.8, zorder=5)
    # Label.
    ax.annotate(
        name,
        xy=(100 * fw, h),
        xytext=(6, 8),
        textcoords="offset points",
        fontsize=10, color=color, fontweight="bold",
    )
    # Ring-FAI precision band around anchor.
    ax.fill_betweenx([h - 2, h + 2], 100 * fw - 1.8, 100 * fw + 1.8,
                     color=color, alpha=0.15, zorder=2)

# Per-voxel noise band around the whole curve.
ax.fill_between(100 * fws, hu_70 - 10, hu_70 + 10,
                color="#cccccc", alpha=0.35, zorder=1,
                label="Per-voxel noise band (~10 HU std)")
# Ring-averaged precision band.
ax.fill_between(100 * fws, hu_70 - 2, hu_70 + 2,
                color="#888888", alpha=0.35, zorder=1,
                label="Ring-averaged precision (~2 HU std)")

# Readable reference lines.
ax.axhline(0, color="#dddddd", lw=0.6)
ax.axhline(-100, color="#dddddd", lw=0.6)
ax.set_xlabel("Water volume fraction in PCAT  (%)", fontsize=11)
ax.set_ylabel("HU at 70 keV VMI", fontsize=11)
ax.set_title("PCAT: volumetric composition -> CT signal\n"
             "(triglyceride + water + collagen mixture model)",
             fontsize=12)
ax.grid(which="both", alpha=0.3)
ax.legend(loc="lower right", fontsize=9, framealpha=0.9)
ax.set_xlim(0, 40)
ax.set_ylim(-115, 10)

plt.tight_layout()
out_path = "/Users/shunie/Developer/PCAT/pcat-workstation-v2/docs/pcat_detectability.png"
plt.savefig(out_path, dpi=180, bbox_inches="tight")
print(f"wrote {out_path}")
print()

# Summary table.
print(f"{'State':<28} {'f_water':>10} {'HU@70':>8} {'HU@150':>8}")
print("-" * 58)
for name, (fw, _) in anchors.items():
    print(f"{name:<28} {fw*100:>9.1f}% {pcat_hu_at_water_fraction(fw,70):>8.1f} "
          f"{pcat_hu_at_water_fraction(fw,150):>8.1f}")

print()
print(f"70 keV sensitivity: {d_hu_per_pp:.2f} HU per percentage-point water")
print(f"Per-voxel detectable delta-fw (SNR=3, sigma=10 HU):   "
      f"{30 / d_hu_per_pp:.1f} pp")
print(f"Ring-averaged detectable delta-fw (SNR=3, sigma=2 HU):"
      f" {6 / d_hu_per_pp:.1f} pp")
