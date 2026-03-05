"""
Fig 4: Individual scattering coefficients |b_0l| for l=0,1,2.

(a) Single cylinder
(b) 11x11 array
(c) 201x201 array

Parameters: R = 0.35a, n = 4, k = q
"""

import sys
sys.path.insert(0, '.')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mie_scattering import (
    single_cylinder_coefficients,
    compute_scattered_coefficients,
)

R_over_a = 0.35
n_refr = 4
l_max = 2
k_factor = 1.001  # k ~ q

# nqR range: 0 to 8 (covers first few Bessel zeros)
nqR_array = np.linspace(0.3, 8.0, 300)

a = 1.0
R = R_over_a * a

# --- (a) Single cylinder ---
print("Computing single cylinder scattering coefficients...")
b_single = {l: np.zeros(len(nqR_array)) for l in range(l_max + 1)}
for idx, nqR in enumerate(nqR_array):
    q = nqR / (n_refr * R)
    coeffs = single_cylinder_coefficients(q, R, n_refr, l_max)
    for l in range(l_max + 1):
        # For l > 0, sum |b_l| + |b_{-l}|
        if l == 0:
            b_single[l][idx] = np.abs(coeffs[0])
        else:
            b_single[l][idx] = np.abs(coeffs[l]) + np.abs(coeffs[-l])

# --- (b) 11x11 array ---
print("Computing 11x11 array scattering coefficients...")
results_11 = compute_scattered_coefficients(
    nqR_array, R_over_a, n_refr, 11, l_max,
    k_factor=k_factor, damping_per_a=0.01
)
b_11 = {l: np.zeros(len(nqR_array)) for l in range(l_max + 1)}
for idx in range(len(nqR_array)):
    for l in range(l_max + 1):
        if l == 0:
            b_11[l][idx] = np.abs(results_11[0][idx])
        else:
            b_11[l][idx] = np.abs(results_11[l][idx]) + np.abs(results_11[-l][idx])

# --- (c) 201x201 array ---
print("Computing 201x201 array scattering coefficients...")
results_201 = compute_scattered_coefficients(
    nqR_array, R_over_a, n_refr, 201, l_max,
    k_factor=k_factor, damping_per_a=0.003
)
b_201 = {l: np.zeros(len(nqR_array)) for l in range(l_max + 1)}
for idx in range(len(nqR_array)):
    for l in range(l_max + 1):
        if l == 0:
            b_201[l][idx] = np.abs(results_201[0][idx])
        else:
            b_201[l][idx] = np.abs(results_201[l][idx]) + np.abs(results_201[-l][idx])

# --- Plotting ---
fig, axes = plt.subplots(3, 1, figsize=(8, 12), sharex=True)
colors = ['b', 'r', 'g']
labels_l = [r'$l=0$', r'$l=1$', r'$l=2$']
titles = [
    r'(a) Single cylinder',
    r'(b) $11\times11$ array',
    r'(c) $201\times201$ array'
]
data_sets = [b_single, b_11, b_201]

for ax, bdata, title in zip(axes, data_sets, titles):
    for l in range(l_max + 1):
        ax.plot(nqR_array, bdata[l], colors[l], linewidth=1, label=labels_l[l])
    ax.set_ylabel(r'$|b_{0l}|$', fontsize=12)
    ax.set_title(title, fontsize=12)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3, linestyle='--')

axes[-1].set_xlabel(r'$nqR$', fontsize=13)
axes[-1].set_xlim(0.3, 8.0)

fig.suptitle(r'Fig. 4: Scattering coefficients ($R=0.35a$, $n=4$)', fontsize=13, y=1.01)
plt.tight_layout()
plt.savefig('../plots/fig4_scattering_coefficients.png', dpi=300, bbox_inches='tight')
print("Saved ../plots/fig4_scattering_coefficients.png")
plt.close()

# Export CSV
with open('../data/fig4_single_cylinder.csv', 'w') as f:
    f.write("# Figure 4a: Single cylinder |b_0l| vs nqR\n")
    f.write("# Parameters: R/a=0.35, n=4\n")
    f.write("nqR,abs_b0_l0,abs_b0_l1,abs_b0_l2\n")
    for i in range(len(nqR_array)):
        f.write(f"{nqR_array[i]:.8f},{b_single[0][i]:.8f},{b_single[1][i]:.8f},{b_single[2][i]:.8f}\n")

with open('../data/fig4_array_11x11.csv', 'w') as f:
    f.write("# Figure 4b: 11x11 array |b_0l| vs nqR\n")
    f.write("# Parameters: R/a=0.35, n=4, k~q\n")
    f.write("nqR,abs_b0_l0,abs_b0_l1,abs_b0_l2\n")
    for i in range(len(nqR_array)):
        f.write(f"{nqR_array[i]:.8f},{b_11[0][i]:.8f},{b_11[1][i]:.8f},{b_11[2][i]:.8f}\n")

with open('../data/fig4_array_201x201.csv', 'w') as f:
    f.write("# Figure 4c: 201x201 array |b_0l| vs nqR\n")
    f.write("# Parameters: R/a=0.35, n=4, k~q\n")
    f.write("nqR,abs_b0_l0,abs_b0_l1,abs_b0_l2\n")
    for i in range(len(nqR_array)):
        f.write(f"{nqR_array[i]:.8f},{b_201[0][i]:.8f},{b_201[1][i]:.8f},{b_201[2][i]:.8f}\n")

print("Exported CSVs for fig4")
