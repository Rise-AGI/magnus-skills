"""
Figure 1: Band structure of 1D photonic crystal (infinite).

Computes band structure using transfer matrix method for a 1D periodic
stack with eps=13 (GaAs) and air, l/P=0.2.

This establishes which frequencies fall in the band gap, essential for
understanding the reflection effect in Fig 6 of the paper.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(sys.argv[0])) if sys.argv[0] else '.')
from photonic_crystal import band_structure_trace, find_bandgaps

# Parameters from the paper (Fig 5)
eps_rod = 13.0
eps_air = 1.0
l_frac = 0.2
d_air = 1.0 - l_frac

# Compute band structure
n_omega = 3000
omega_max = 0.6
omegas = np.linspace(0.001, omega_max, n_omega)

band_data = []
for omega in omegas:
    ht = band_structure_trace(omega, eps_rod, eps_air, l_frac, d_air)
    if abs(ht) <= 1.0:
        ky = np.arccos(np.clip(ht, -1.0, 1.0))
        band_data.append((ky / np.pi, omega))
        if ky > 0.01:
            band_data.append((1.0 - ky / np.pi, omega))

band_data = np.array(band_data) if band_data else np.zeros((0, 2))

# Find band gaps
gap_regions = find_bandgaps(eps_rod, eps_air, l_frac, d_air, n_omega=3000, omega_max=omega_max)

# Save data
script_dir = os.path.dirname(os.path.abspath(sys.argv[0])) if sys.argv[0] else '.'
data_dir = os.path.join(script_dir, '..', 'data')
plots_dir = os.path.join(script_dir, '..', 'plots')
os.makedirs(data_dir, exist_ok=True)
os.makedirs(plots_dir, exist_ok=True)

np.savetxt(os.path.join(data_dir, 'fig1_bandstructure.csv'),
           band_data, delimiter=',',
           header='ky_over_pi,omega_normalized',
           comments='# Band structure of 1D photonic crystal (eps=13/air, l/P=0.2)\n# ')

with open(os.path.join(data_dir, 'fig1_bandgaps.csv'), 'w') as f:
    f.write('# Band gaps of 1D photonic crystal (eps=13/air, l/P=0.2)\n')
    f.write('gap_low,gap_high\n')
    for gl, gh in gap_regions:
        f.write(f'{gl:.8f},{gh:.8f}\n')

# Plot
fig, ax = plt.subplots(1, 1, figsize=(6, 8))
if len(band_data) > 0:
    ax.scatter(band_data[:, 0], band_data[:, 1], s=0.3, c='blue', alpha=0.5)
for i, (gl, gh) in enumerate(gap_regions):
    ax.axhspan(gl, gh, alpha=0.15, color='red',
               label='Band gap' if i == 0 else '')
ax.axhline(y=0.3, color='green', linestyle='--', alpha=0.5, label=r'$\omega P/(2\pi c)=0.3$')
ax.set_xlabel(r'$k_y P / \pi$')
ax.set_ylabel(r'$\omega P / (2\pi c)$')
ax.set_title('1D Photonic Crystal Band Structure\n'
             r'$\varepsilon_{rod}=13$ (GaAs), $\varepsilon_{air}=1$, $l/P=0.2$')
ax.set_xlim(0, 1)
ax.set_ylim(0, omega_max)
ax.legend()
ax.grid(True, alpha=0.3)
fig.savefig(os.path.join(plots_dir, 'fig1_bandstructure.png'), dpi=150, bbox_inches='tight')
plt.close(fig)

print(f"Band gaps found: {len(gap_regions)}")
for i, (gl, gh) in enumerate(gap_regions):
    print(f"  Gap {i+1}: {gl:.4f} - {gh:.4f} (width: {gh-gl:.4f})")
print(f"omega=0.3 is in band gap: {not any(abs(band_structure_trace(0.3, eps_rod, eps_air, l_frac, d_air)) <= 1.0 for _ in [1])}")
print("Done: fig1_bandstructure.png")
