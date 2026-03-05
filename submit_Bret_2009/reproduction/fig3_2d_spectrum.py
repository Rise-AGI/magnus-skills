"""
Figure 3: 2D unstable spectrum in (Zx, Zz) plane.
Three panels for R=0, R=1/100, R=1/1836.

Parameters: alpha=0.1, gamma_b=2, OmegaB=3.
Reproduces Fig. 3 of Bret (2009).
"""
import sys
import os
import numpy as np

sys.path.insert(0, os.path.dirname(__file__) if "__file__" in dir() else ".")
from dielectric_tensor import compute_2d_spectrum, plasma_params

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Parameters (Fig. 3 caption)
alpha = 0.1
gamma_b = 2.0
OmegaB = 3.0

# Grid: Zx from 0 to ~10, Zz from 0 to ~10 (Fig. 3 shows range ~0-6)
nZx = 50
nZz = 60
Zx_arr = np.linspace(0.0, 8.0, nZx)
Zz_arr = np.linspace(0.01, 8.0, nZz)

os.makedirs("../data", exist_ok=True)
os.makedirs("../plots", exist_ok=True)

R_values = {"R0": 0.0, "R100": 1.0/100, "R1836": 1.0/1836}
panels = {}

for label, R in R_values.items():
    print(f"Computing 2D spectrum for {label} (R={R:.6f})...")
    if R == 0:
        # R=0 means no proton motion; set R to very small value
        R_eff = 1e-10
    else:
        R_eff = R
    delta = compute_2d_spectrum(Zx_arr, Zz_arr, alpha, gamma_b, OmegaB, R_eff,
                                verbose=True)
    panels[label] = delta

    # Save data
    np.savetxt(f"../data/fig3_spectrum_{label}.csv", delta, delimiter=",",
               header=f"2D growth rate map. Rows=Zx(0..8,{nZx}pts), Cols=Zz(0.01..8,{nZz}pts). alpha={alpha},gamma_b={gamma_b},OmegaB={OmegaB},R={R}",
               fmt="%.6e")

# Plot
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
titles = [r"(A) $R=0$", r"(B) $R=1/100$", r"(C) $R=1/1836$"]
for ax, (label, delta), title in zip(axes, panels.items(), titles):
    vmax = delta.max() if delta.max() > 0 else 0.1
    im = ax.contourf(Zz_arr, Zx_arr, delta, levels=20, cmap="hot_r")
    ax.set_xlabel(r"$Z_z$")
    ax.set_ylabel(r"$Z_x$")
    ax.set_title(title)
    plt.colorbar(im, ax=ax, label=r"$\delta$")

fig.suptitle(r"Fig 3: $\alpha=0.1,\ \gamma_b=2,\ \Omega_B=3$", fontsize=14)
fig.tight_layout()
fig.savefig("../plots/fig3_2d_spectrum.png", dpi=150)
plt.close()

print("Done. Saved data/fig3_spectrum_*.csv and plots/fig3_2d_spectrum.png")
