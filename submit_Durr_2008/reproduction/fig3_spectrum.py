"""
Figure 3: The light hadron spectrum of QCD.

Reproduces Figure 3 from Durr et al., Science 322, 1224 (2008).

Plots lattice QCD results (Xi-set) for hadron masses compared to
experimental values. Horizontal lines and bands show experimental
values with decay widths for resonances.

Output:
    ../plots/fig3_spectrum.png
    ../data/fig3_spectrum.csv
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, ".")
from lattice_qcd import (
    EXPERIMENTAL_MASSES, LATTICE_XI_SET, DECAY_WIDTHS,
    INPUT_PARTICLES, SCALE_PARTICLES_XI, PREDICTED_PARTICLES_XI,
)

# Particle ordering for the plot (matching Fig 3)
PLOT_ORDER = ["pi", "K", "rho", "K*", "N", "Lambda", "Sigma", "Xi",
              "Delta", "Sigma*", "Xi*", "Omega"]

print("=" * 60)
print("Figure 3: Light Hadron Spectrum of QCD")
print("=" * 60)

# ---- Collect data ----
x_pos = np.arange(len(PLOT_ORDER))
exp_masses = [EXPERIMENTAL_MASSES[p] for p in PLOT_ORDER]

lat_central = []
lat_err = []
is_input = []
is_scale = []

for p in PLOT_ORDER:
    if p in LATTICE_XI_SET:
        c, s, sy = LATTICE_XI_SET[p]
        total_err = np.sqrt(s**2 + sy**2)
        lat_central.append(c)
        lat_err.append(total_err)
    else:
        lat_central.append(EXPERIMENTAL_MASSES[p])
        lat_err.append(0.0)

    is_input.append(p in INPUT_PARTICLES)
    is_scale.append(p in SCALE_PARTICLES_XI)

lat_central = np.array(lat_central)
lat_err = np.array(lat_err)

# ---- Plot ----
fig, ax = plt.subplots(1, 1, figsize=(10, 6))

# Experimental values: lines + decay width bands
for i, p in enumerate(PLOT_ORDER):
    m = EXPERIMENTAL_MASSES[p]
    ax.hlines(m, i - 0.35, i + 0.35, colors="black", linewidth=1.0)
    if p in DECAY_WIDTHS:
        w = DECAY_WIDTHS[p]
        ax.fill_between([i - 0.35, i + 0.35], m - w / 2, m + w / 2,
                        color="gold", alpha=0.4, zorder=1)

# Lattice results
for i, p in enumerate(PLOT_ORDER):
    if p in INPUT_PARTICLES or p in SCALE_PARTICLES_XI:
        # Input particles: no error bar
        ax.plot(i, lat_central[i], "o", color="royalblue", markersize=7, zorder=3)
    else:
        ax.errorbar(i, lat_central[i], yerr=lat_err[i],
                     fmt="o", color="royalblue", markersize=7,
                     capsize=3, capthick=1.5, elinewidth=1.5, zorder=3)

# Labels
ax.set_xticks(x_pos)
ax.set_xticklabels([r"$\pi$", r"$K$", r"$\rho$", r"$K^*$",
                     r"$N$", r"$\Lambda$", r"$\Sigma$", r"$\Xi$",
                     r"$\Delta$", r"$\Sigma^*$", r"$\Xi^*$", r"$\Omega$"],
                    fontsize=12)

ax.set_ylabel("Mass [GeV]", fontsize=13)
ax.set_title("Light Hadron Spectrum of QCD (Xi-set)", fontsize=14)

# Add separator between mesons and baryons
ax.axvline(3.5, color="gray", linestyle="--", alpha=0.5)
ax.axvline(7.5, color="gray", linestyle="--", alpha=0.5)
ax.text(1.5, 0.05, "Mesons", ha="center", fontsize=10, color="gray")
ax.text(5.5, 0.05, "Octet", ha="center", fontsize=10, color="gray")
ax.text(9.5, 0.05, "Decuplet", ha="center", fontsize=10, color="gray")

ax.set_ylim(0, 2.0)
ax.set_xlim(-0.5, len(PLOT_ORDER) - 0.5)
ax.grid(axis="y", alpha=0.3)
plt.tight_layout()
plt.savefig("../plots/fig3_spectrum.png", dpi=300)
print("Saved: ../plots/fig3_spectrum.png")

# ---- Export CSV ----
with open("../data/fig3_spectrum.csv", "w") as f:
    f.write("# Figure 3: Light hadron spectrum - lattice QCD vs experiment\n")
    f.write("# Columns: particle, exp_mass_GeV, lat_mass_GeV, lat_stat_err, lat_sys_err, lat_total_err\n")
    f.write("particle,exp_mass_GeV,lat_mass_GeV,lat_stat_err,lat_sys_err,lat_total_err\n")
    for p in PLOT_ORDER:
        exp_m = EXPERIMENTAL_MASSES[p]
        if p in LATTICE_XI_SET:
            c, s, sy = LATTICE_XI_SET[p]
            total = np.sqrt(s**2 + sy**2)
            f.write(f"{p},{exp_m:.8f},{c:.8f},{s:.8f},{sy:.8f},{total:.8f}\n")
        else:
            f.write(f"{p},{exp_m:.8f},{exp_m:.8f},0.00000000,0.00000000,0.00000000\n")

print("Saved: ../data/fig3_spectrum.csv")

# Print comparison table
print("\n" + "-" * 70)
print(f"{'Particle':>10s} {'Exp (GeV)':>10s} {'Lat (GeV)':>12s} {'Diff (MeV)':>12s} {'sigma':>8s}")
print("-" * 70)
for p in PREDICTED_PARTICLES_XI:
    exp_m = EXPERIMENTAL_MASSES[p]
    c, s, sy = LATTICE_XI_SET[p]
    total = np.sqrt(s**2 + sy**2)
    diff = (c - exp_m) * 1000
    sigma = abs(diff) / (total * 1000) if total > 0 else 0
    print(f"{p:>10s} {exp_m:10.3f} {c:8.3f}({s*1000:.0f})({sy*1000:.0f}) {diff:+10.1f} {sigma:8.1f}")

print("-" * 70)
print("Done.")
