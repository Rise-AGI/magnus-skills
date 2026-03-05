"""
Figure 3: Fermion determinant ratio Det^N vs |N|.
Reproduces Fig. 5 from Fukaya & Onogi (2003).

Det^N = integral(det(D_DW)^2/det(D_AP)^2) for sector N
        / integral(det(D_DW)^2/det(D_AP)^2) for sector 0

Computed on classical gauge configurations with moduli integration.
Shows results for fermion masses m=0.1 and m=0.2.
"""
import sys
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(__file__))
from schwinger_luscher import compute_DetN

os.makedirs("../data", exist_ok=True)
os.makedirs("../plots", exist_ok=True)

# Parameters from the paper
L = 16
L3 = 6
M = 0.9
n_nu = 3  # Use 3x3 moduli points (paper uses 5x5 but that's slow)

N_values = [0, 1, 2, 3, 4]
masses = [0.1, 0.2]

results = {}
for m_f in masses:
    print(f"\nComputing Det^N for m={m_f}...")
    det_vals = []
    for N in N_values:
        print(f"  N={N}...", end=" ", flush=True)
        det_N = compute_DetN(L, L3, M, m_f, N, n_nu=n_nu)
        det_vals.append(det_N)
        print(f"Det^{N} = {det_N:.8e}")
    results[m_f] = det_vals

# --- Save data ---
data = np.column_stack([
    N_values,
    results[0.1],
    results[0.2]
])
header = ("# Fermion determinant ratio Det^N vs |N|\n"
          "# L=16, L3=6, M=0.9, n_nu=3x3 moduli integration\n"
          "# Columns: N, DetN_m0.1, DetN_m0.2\n"
          "N,DetN_m0.1,DetN_m0.2")
np.savetxt("../data/fig3_fermion_determinant.csv", data, delimiter=",",
           header=header, comments="", fmt="%.8e")
print("Saved data/fig3_fermion_determinant.csv")

# --- Plot ---
fig, ax = plt.subplots(figsize=(8, 6))
for m_f in masses:
    ax.semilogy(N_values, results[m_f], 'o-', markersize=8,
                label=f"$m = {m_f}$")

ax.set_xlabel(r"$|N|$", fontsize=14)
ax.set_ylabel(r"$Det^N$", fontsize=14)
ax.set_title(r"Fermion determinant ratio $Det^N$ vs $|N|$")
ax.legend(fontsize=12)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("../plots/fig3_fermion_determinant.png", dpi=150, bbox_inches="tight")
print("Saved plots/fig3_fermion_determinant.png")
