"""
Figure 1: Growth rate delta (omega_p units) of four unstable modes
for flow-aligned wave vectors (Zx=0) vs Zz = k_par * vb / omega_p.

Parameters: alpha=0.1, gamma_b=20, OmegaB=1, R=1/1836.
Reproduces Fig. 1 of Bret (2009).
"""
import sys
import os
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(".")))
sys.path.insert(0, os.path.dirname(__file__) if "__file__" in dir() else ".")
from dielectric_tensor import (
    scan_parallel_modes, plasma_params, newton_complex,
    Tzz_par, Txx_pm_iTxy_par
)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Parameters (from Fig. 1 caption)
alpha = 0.1
gamma_b = 20.0
OmegaB = 1.0
R = 1.0 / 1836.0
beta, gamma_p = plasma_params(alpha, gamma_b, R)

# Zz range: the spectrum extends over ~3 orders of magnitude (Fig. 1)
# Two-Stream peaks near Zz ~ 1, Buneman near Zz ~ 1/alpha = 10
Zz_fine = np.logspace(-1, np.log10(30), 300)

print("Computing parallel growth rates for Fig 1...")
print(f"  alpha={alpha}, gamma_b={gamma_b}, OmegaB={OmegaB}, R={R:.6f}")

delta_es, delta_em_m, delta_em_p = scan_parallel_modes(
    Zz_fine, alpha, gamma_b, OmegaB, R)

# Save data
os.makedirs("../data", exist_ok=True)
header = ("Zz,delta_electrostatic,delta_em_minus,delta_em_plus\n"
          "# Fig 1: alpha=0.1, gamma_b=20, OmegaB=1, R=1/1836\n"
          "# Electrostatic = Two-Stream + Buneman merged\n"
          "# EM minus = Txx - iTxy (Bell-like)\n"
          "# EM plus = Txx + iTxy\n")
data = np.column_stack([Zz_fine, delta_es, delta_em_m, delta_em_p])
np.savetxt("../data/fig1_parallel_growth.csv", data, delimiter=",",
           header="Zz,delta_electrostatic,delta_em_minus,delta_em_plus",
           fmt="%.8e")

# Plot
os.makedirs("../plots", exist_ok=True)
fig, ax = plt.subplots(figsize=(8, 5))
mask_es = delta_es > 1e-6
mask_em_m = delta_em_m > 1e-6
mask_em_p = delta_em_p > 1e-6

if mask_es.any():
    ax.semilogy(Zz_fine[mask_es], delta_es[mask_es], "k-", lw=2, label="Electrostatic (TS+Bun)")
if mask_em_m.any():
    ax.semilogy(Zz_fine[mask_em_m], delta_em_m[mask_em_m], "b-", lw=1, label=r"$T_{xx}-iT_{xy}$ (Bell-like)")
if mask_em_p.any():
    ax.semilogy(Zz_fine[mask_em_p], delta_em_p[mask_em_p], "b--", lw=1, label=r"$T_{xx}+iT_{xy}$")

# Analytical predictions
delta_ts = np.sqrt(3) / 2**(4./3) * alpha**(1./3) / gamma_b
delta_bun = np.sqrt(3) / 2**(4./3) * R**(1./3)
ax.axhline(delta_ts, color="r", ls=":", alpha=0.5, label=f"TS analytical ({delta_ts:.4f})")
ax.axhline(delta_bun, color="g", ls=":", alpha=0.5, label=f"Bun analytical ({delta_bun:.4f})")

ax.set_xlabel(r"$Z_z = k_\parallel v_b / \omega_p$")
ax.set_ylabel(r"$\delta$ ($\omega_p$ units)")
ax.set_title(r"Fig 1: $\alpha=0.1,\ \gamma_b=20,\ \Omega_B=1,\ R=1/1836$")
ax.legend(fontsize=8)
ax.set_xlim(0.1, 30)
ax.set_ylim(1e-4, 0.1)
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig("../plots/fig1_parallel_growth.png", dpi=150)
plt.close()

print(f"Done. Max electrostatic delta = {delta_es.max():.6f}")
print(f"  Expected Two-Stream: {delta_ts:.6f}")
print(f"  Expected Buneman: {delta_bun:.6f}")
print(f"Saved data/fig1_parallel_growth.csv and plots/fig1_parallel_growth.png")
