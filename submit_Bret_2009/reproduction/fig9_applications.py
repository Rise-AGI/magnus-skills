"""
Figure 9: 2D unstable spectra for astrophysical applications.
Three panels: Solar Flares, Intergalactic Streams, Relativistic Shocks.

Parameters from Table 3:
  Solar Flares: alpha=1/8, gamma_b=1.34, OmegaB=0 (and 1), R=1/1836
  Intergalactic: alpha=1, gamma_b=1.005 (beta=1/10), OmegaB=0, R=1/1836
  Rel. Shocks: alpha=1, gamma_b=10, OmegaB=0, R=0 (pair plasma)

Reproduces Fig. 9 of Bret (2009).
"""
import sys
import os
import numpy as np

sys.path.insert(0, os.path.dirname(__file__) if "__file__" in dir() else ".")
from dielectric_tensor import compute_2d_spectrum

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

os.makedirs("../data", exist_ok=True)
os.makedirs("../plots", exist_ok=True)

# Application parameters (Table 3)
cases = {
    "solar_flares": {"alpha": 1.0/8, "gamma_b": 1.34, "OmegaB": 1.0, "R": 1.0/1836},
    "intergalactic": {"alpha": 1.0, "gamma_b": 1.005, "OmegaB": 0.0, "R": 1.0/1836},
    "rel_shocks": {"alpha": 1.0, "gamma_b": 10.0, "OmegaB": 0.0, "R": 1e-10},  # R=0 pair plasma
}

nZx = 40
nZz = 50

fig, axes = plt.subplots(1, 3, figsize=(15, 5))
titles = ["(A) Solar Flares", "(B) Intergalactic Streams", "(C) Relativistic Shocks"]

for idx, (name, params) in enumerate(cases.items()):
    alpha = params["alpha"]
    gamma_b = params["gamma_b"]
    OmegaB = params["OmegaB"]
    R = params["R"]

    # Adaptive grid range based on physics
    if name == "solar_flares":
        Zx_max, Zz_max = 5.0, 5.0
    elif name == "intergalactic":
        Zx_max, Zz_max = 3.0, 3.0
    else:
        Zx_max, Zz_max = 5.0, 5.0

    Zx_arr = np.linspace(0.0, Zx_max, nZx)
    Zz_arr = np.linspace(0.01, Zz_max, nZz)

    print(f"Computing {name}: alpha={alpha}, gamma_b={gamma_b}, OmegaB={OmegaB}, R={R:.2e}")
    delta = compute_2d_spectrum(Zx_arr, Zz_arr, alpha, gamma_b, OmegaB, R, verbose=True)

    # Save
    np.savetxt(f"../data/fig9_{name}.csv", delta, delimiter=",",
               header=f"Rows=Zx(0..{Zx_max},{nZx}pts), Cols=Zz(0.01..{Zz_max},{nZz}pts). {params}",
               fmt="%.6e")

    # Plot
    ax = axes[idx]
    vmax = delta.max() if delta.max() > 0 else 0.1
    im = ax.contourf(Zz_arr, Zx_arr, delta, levels=20, cmap="hot_r")
    ax.set_xlabel(r"$Z_z$")
    ax.set_ylabel(r"$Z_x$")
    ax.set_title(titles[idx])
    plt.colorbar(im, ax=ax, label=r"$\delta$")

fig.suptitle("Fig 9: 2D spectra for astrophysical applications", fontsize=14)
fig.tight_layout()
fig.savefig("../plots/fig9_applications.png", dpi=150)
plt.close()

print("Done. Saved fig9 data and plots.")
