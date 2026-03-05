"""
Figure 2: Upper bound transmittance Tu vs bend radius R.

Reproduces Figure 2 from Hasegawa et al., Appl. Phys. Lett. 84, 1835 (2004).
Silver-air interface, theta = 90 degrees, wavelengths 500, 600, 700 nm.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(".")))
sys.path.insert(0, os.path.dirname(__file__) if '__file__' in dir() else '.')

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from spp_core import compute_Tu_vs_R, silver_epsilon, c

theta = np.pi / 2  # 90 degrees
wavelengths = [500, 600, 700]
labels = ['500 nm', '600 nm', '700 nm']
styles = ['-.', '--', '-']

# Main curves: 120 points per wavelength
R_array = np.logspace(np.log10(1e-6), np.log10(200e-6), 120)
R_um = R_array * 1e6

data_dir = os.path.join(os.path.dirname(__file__) if '__file__' in dir() else '.', '..', 'data')
plots_dir = os.path.join(os.path.dirname(__file__) if '__file__' in dir() else '.', '..', 'plots')
os.makedirs(data_dir, exist_ok=True)
os.makedirs(plots_dir, exist_ok=True)

all_Tu = {}
for wl in wavelengths:
    eps_i = silver_epsilon(wl)
    print(f"Computing Tu for lambda={wl}nm (eps_i={eps_i:.3f})...")
    Tu, m_arr = compute_Tu_vs_R(wl, R_array, theta)
    all_Tu[wl] = Tu
    idx_peak = np.nanargmax(Tu)
    print(f"  Peak Tu={Tu[idx_peak]:.4f} at R={R_um[idx_peak]:.1f} um")

# Save CSV
csv_path = os.path.join(data_dir, 'fig2_transmittance_bound.csv')
header = 'R_um,' + ','.join(f'Tu_{wl}nm' for wl in wavelengths)
data_matrix = np.column_stack([R_um] + [all_Tu[wl] for wl in wavelengths])
np.savetxt(csv_path, data_matrix, delimiter=',', header=header, comments='', fmt='%.8e')
print(f"Data saved to {csv_path}")

# Plot main figure
fig, ax = plt.subplots(1, 1, figsize=(8, 6))
for wl, style, label in zip(wavelengths, styles, labels):
    ax.plot(R_um, all_Tu[wl], style, label=r'$\lambda = $' + label, linewidth=1.5)

ax.set_xlabel(r'Bend Radius R ($\mu$m)', fontsize=12)
ax.set_ylabel(r'$T_u$', fontsize=14)
ax.set_title(r'Upper Bound Transmittance $T_u$ (Silver-Air, $\theta = 90°$)')
ax.legend(fontsize=11)
ax.set_xlim([0, 200])
ax.set_ylim([0, 1])
ax.grid(True, alpha=0.3)

# Inset: Tu map (coarser grid for speed)
print("Computing inset Tu map...")
wl_fine = np.linspace(400, 800, 25)
R_fine = np.logspace(np.log10(1e-6), np.log10(200e-6), 40)
Tu_map = np.zeros((len(wl_fine), len(R_fine)))
for i, wl in enumerate(wl_fine):
    Tu_row, _ = compute_Tu_vs_R(wl, R_fine, theta)
    Tu_map[i, :] = Tu_row
    print(f"  Inset: {i+1}/{len(wl_fine)}")

R_fine_um = R_fine * 1e6
ax_inset = fig.add_axes([0.55, 0.45, 0.35, 0.35])
ax_inset.pcolormesh(R_fine_um, wl_fine, Tu_map, cmap='gray_r', shading='auto', vmin=0, vmax=1)
ax_inset.set_xlabel(r'R ($\mu$m)', fontsize=9)
ax_inset.set_ylabel(r'$\lambda$ (nm)', fontsize=9)
ax_inset.set_title(r'$T_u$', fontsize=10)
ax_inset.tick_params(labelsize=8)

plot_path = os.path.join(plots_dir, 'fig2_transmittance_bound.png')
plt.savefig(plot_path, dpi=150, bbox_inches='tight')
print(f"Plot saved to {plot_path}")

# Save inset data
inset_path = os.path.join(data_dir, 'fig2_inset_Tu_map.csv')
with open(inset_path, 'w') as f:
    f.write(f"# Tu map: rows=wavelength(nm), cols=R(um)\n")
    f.write('# wavelengths_nm: ' + ','.join(f'{w:.1f}' for w in wl_fine) + '\n')
    f.write('# R_um: ' + ','.join(f'{r:.4f}' for r in R_fine_um) + '\n')
    for i in range(len(wl_fine)):
        f.write(','.join(f'{Tu_map[i,j]:.8e}' for j in range(len(R_fine))) + '\n')
print(f"Inset data saved to {inset_path}")
