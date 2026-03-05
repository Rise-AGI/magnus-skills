"""
Figure 2: Minimum classical action S_Gmin^N vs topological charge |N|^2.
Reproduces Fig. 3 from Fukaya & Onogi (2003).

For Luscher's action, the minimum action in each topological sector
is computed from the classical configuration with uniform field strength.
S_Gmin^N scales quadratically with |N|.
"""
import sys
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(__file__))
from schwinger_luscher import luscher_action_min

os.makedirs("../data", exist_ok=True)
os.makedirs("../plots", exist_ok=True)

L = 16
beta = 0.5
epsilon = 1.0

N_values = np.arange(0, 11)
S_min = np.array([luscher_action_min(L, N, beta, epsilon) for N in N_values])
N_sq = N_values ** 2

# Fit to quadratic S = a * N^2
mask = S_min < np.inf
coeffs = np.polyfit(N_sq[mask], S_min[mask], 1)
S_fit = np.polyval(coeffs, N_sq)

print("S_Gmin^N values:")
for N, S in zip(N_values, S_min):
    print(f"  N={N}: S_Gmin = {S:.8f}")
print(f"Quadratic fit: S = {coeffs[0]:.8f} * N^2 + {coeffs[1]:.8f}")

# --- Save data ---
data = np.column_stack([N_values, N_sq, S_min, S_fit])
header = ("# Minimum classical action S_Gmin^N vs topological charge\n"
          "# Luscher action, L=16, beta=0.5, epsilon=1.0\n"
          "# Columns: N, N_squared, S_Gmin, S_fit_quadratic\n"
          "N,N_squared,S_Gmin,S_fit_quadratic")
np.savetxt("../data/fig2_classical_action.csv", data, delimiter=",",
           header=header, comments="", fmt="%.8f")
print("Saved data/fig2_classical_action.csv")

# --- Plot ---
fig, ax = plt.subplots(figsize=(8, 6))
ax.plot(N_sq[mask], S_min[mask], 'ko', markersize=8, label="Data")
ax.plot(N_sq, S_fit, 'r--', linewidth=1.5,
        label=f"Fit: $S = {coeffs[0]:.4f} N^2 + {coeffs[1]:.4f}$")
ax.set_xlabel(r"$|N|^2$", fontsize=14)
ax.set_ylabel(r"$S_{G,\mathrm{min}}^{N}$", fontsize=14)
ax.set_title(r"Minimum action vs $|N|^2$ (L\"uscher action, $\beta=0.5$, $\epsilon=1.0$)")
ax.legend(fontsize=12)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("../plots/fig2_classical_action.png", dpi=150, bbox_inches="tight")
print("Saved plots/fig2_classical_action.png")
