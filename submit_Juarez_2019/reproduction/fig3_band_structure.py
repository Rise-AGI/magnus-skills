"""
Fig 3 (lower panel): Band structure from scattering amplitude.

Computes |b_0|^2 as a function of both nqR and ka, displayed as a 2D
color map (logarithmic scale) showing the photonic band structure.

Parameters: R = 0.35a, n = 4, 201x201 array, l_max = 2.
Smoothed as eta/(eta^2 + 1/|b_0|^2) with eta = 0.1.
"""

import sys
sys.path.insert(0, '.')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mie_scattering import solve_bloch_scattering

R_over_a = 0.35
n_refr = 4
l_max = 2
N_array = 201
eta = 0.1

a = 1.0
R = R_over_a * a

# nqR range: 0.5 to 7.0  (vertical axis: frequency)
# ka range: 0.05 to pi    (horizontal axis: wavevector)
n_freq = 120
n_kvec = 80

nqR_vals = np.linspace(0.5, 7.0, n_freq)
ka_vals = np.linspace(0.05, np.pi, n_kvec)

# Output array: smoothed scattering amplitude
intensity = np.zeros((n_freq, n_kvec))

print(f"Computing band structure: {n_freq} x {n_kvec} = {n_freq * n_kvec} points")
print(f"Parameters: R/a={R_over_a}, n={n_refr}, l_max={l_max}, N={N_array}")

damping_per_a = 0.003

for j, ka in enumerate(ka_vals):
    if j % 10 == 0:
        print(f"  ka column {j+1}/{n_kvec} (ka={ka:.3f})")
    k_x = ka / a
    for i, nqR in enumerate(nqR_vals):
        q = nqR / (n_refr * R)
        try:
            beta = solve_bloch_scattering(q, k_x, R, n_refr, a, N_array, l_max, damping_per_a / a)
            b0 = sum(beta[l] for l in range(-l_max, l_max + 1))
            b0_sq = np.abs(b0)**2
            # Smoothing: eta / (eta^2 + 1/|b_0|^2)
            if b0_sq > 1e-30:
                intensity[i, j] = eta / (eta**2 + 1.0 / b0_sq)
            else:
                intensity[i, j] = 0.0
        except (np.linalg.LinAlgError, ValueError, FloatingPointError):
            intensity[i, j] = 0.0

# Replace zeros with small value for log scale
intensity[intensity <= 0] = 1e-15

fig, ax = plt.subplots(figsize=(8, 6))
im = ax.pcolormesh(
    ka_vals, nqR_vals, np.log10(intensity),
    shading='auto', cmap='hot'
)
ax.set_xlabel(r'$ka$', fontsize=13)
ax.set_ylabel(r'$nqR$', fontsize=13)
ax.set_title(r'Fig. 3 (lower): Scattering amplitude band structure ($R=0.35a$, $n=4$)', fontsize=11)
plt.colorbar(im, label=r'$\log_{10}[\eta/(\eta^2 + 1/|b_0|^2)]$')
ax.set_xlim(ka_vals[0], ka_vals[-1])
ax.set_ylim(nqR_vals[0], nqR_vals[-1])
plt.tight_layout()
plt.savefig('../plots/fig3_band_structure.png', dpi=300, bbox_inches='tight')
print("Saved ../plots/fig3_band_structure.png")
plt.close()

# Export CSV
with open('../data/fig3_band_structure.csv', 'w') as f:
    f.write("# Figure 3 (lower panel): Band structure from scattering amplitude\n")
    f.write("# Parameters: R/a=0.35, n=4, 201x201 array, l_max=2, eta=0.1\n")
    f.write("# Rows: nqR values, Columns: ka values\n")
    f.write("# First row: ka values\n")
    f.write("nqR_ka," + ",".join(f"{ka:.8f}" for ka in ka_vals) + "\n")
    for i, nqR in enumerate(nqR_vals):
        row = ",".join(f"{intensity[i,j]:.8e}" for j in range(n_kvec))
        f.write(f"{nqR:.8f},{row}\n")
print("Exported ../data/fig3_band_structure.csv")
