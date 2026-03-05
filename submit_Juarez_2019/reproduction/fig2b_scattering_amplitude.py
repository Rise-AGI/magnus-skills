"""
Fig 2b: Scattering amplitude |b_0|^2 for thin cylinders.

Parameters: R = 0.1a, n = 10, 511x511 array, l_max = 1, k ~ q
Plots |b_0|^2 = |sum_l b_0l|^2 vs nqR.
"""

import sys
sys.path.insert(0, '.')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mie_scattering import compute_scattered_coefficients

# Parameters from the paper
R_over_a = 0.1
n_refr = 10
N_array = 511
l_max = 1
k_factor = 1.001  # k approx q

# Frequency range: nqR from 0.5 to 8
nqR_array = np.linspace(0.5, 8.0, 400)
damping_per_a = 0.005

print("Computing scattering amplitude for 511x511 array (k ~ q)...")

results = compute_scattered_coefficients(
    nqR_array, R_over_a, n_refr, N_array, l_max,
    k_factor=k_factor, damping_per_a=damping_per_a
)

# |b_0|^2 = |sum_l b_0l|^2
b0_total = np.zeros(len(nqR_array), dtype=complex)
for l in range(-l_max, l_max + 1):
    b0_total += results[l]
b0_sq = np.abs(b0_total)**2

fig, ax = plt.subplots(figsize=(7, 5))
ax.semilogy(nqR_array, b0_sq, 'b-', linewidth=1)
ax.set_xlabel(r'$nqR$', fontsize=13)
ax.set_ylabel(r'$|b_0|^2$', fontsize=13)
ax.set_title(r'Fig. 2(b): Scattering amplitude ($R=0.1a$, $n=10$, $511\times511$)', fontsize=11)
ax.set_xlim(0.5, 8.0)
ax.grid(True, alpha=0.3, linestyle='--')
plt.tight_layout()
plt.savefig('../plots/fig2b_scattering_amplitude.png', dpi=300, bbox_inches='tight')
print("Saved ../plots/fig2b_scattering_amplitude.png")
plt.close()

# Export CSV
with open('../data/fig2b_scattering_amplitude.csv', 'w') as f:
    f.write("# Figure 2b: |b_0|^2 vs nqR\n")
    f.write("# Parameters: R/a=0.1, n=10, 511x511 array, l_max=1, k~q\n")
    f.write("nqR,b0_squared\n")
    for i in range(len(nqR_array)):
        f.write(f"{nqR_array[i]:.8f},{b0_sq[i]:.8f}\n")
print("Exported ../data/fig2b_scattering_amplitude.csv")
