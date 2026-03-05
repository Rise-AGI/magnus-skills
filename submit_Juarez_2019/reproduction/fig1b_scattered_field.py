"""
Fig 1b: Scattered field coefficients for thin cylinders.

Parameters: R = 0.1a, n = 10, 511x511 array, l_max = 1, k = 1.01q
Plots sum of Im(b_0l) vs nqR.
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
k_factor = 1.01

# Frequency range: nqR from 0.5 to 8
nqR_array = np.linspace(0.5, 8.0, 400)

# Damping to mimic infinite array (eliminate edge reflections)
damping_per_a = 0.005

print("Computing scattering coefficients for 511x511 array...")
print(f"Parameters: R/a={R_over_a}, n={n_refr}, l_max={l_max}, k={k_factor}q")
print(f"nqR range: {nqR_array[0]:.1f} to {nqR_array[-1]:.1f}, {len(nqR_array)} points")

results = compute_scattered_coefficients(
    nqR_array, R_over_a, n_refr, N_array, l_max,
    k_factor=k_factor, damping_per_a=damping_per_a
)

# Sum of Im(b_0l) over all l
sum_imag = np.zeros(len(nqR_array))
for l in range(-l_max, l_max + 1):
    sum_imag += results[l].imag

fig, ax = plt.subplots(figsize=(7, 5))
ax.plot(nqR_array, sum_imag, 'b-', linewidth=1)
ax.set_xlabel(r'$nqR$', fontsize=13)
ax.set_ylabel(r"$\sum_l b''_{0l}$", fontsize=13)
ax.set_title(r'Fig. 1(b): Scattered field coefficients ($R=0.1a$, $n=10$, $511\times511$)', fontsize=11)
ax.set_xlim(0.5, 8.0)
ax.grid(True, alpha=0.3, linestyle='--')
plt.tight_layout()
plt.savefig('../plots/fig1b_scattered_field.png', dpi=300, bbox_inches='tight')
print("Saved ../plots/fig1b_scattered_field.png")
plt.close()

# Export CSV
with open('../data/fig1b_scattered_field.csv', 'w') as f:
    f.write("# Figure 1b: Sum of Im(b_0l) vs nqR\n")
    f.write("# Parameters: R/a=0.1, n=10, 511x511 array, l_max=1, k=1.01q\n")
    f.write("nqR,sum_Im_b0l\n")
    for i in range(len(nqR_array)):
        f.write(f"{nqR_array[i]:.8f},{sum_imag[i]:.8f}\n")
print("Exported ../data/fig1b_scattered_field.csv")
