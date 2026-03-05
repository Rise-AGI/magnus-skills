"""
Figure 1: Phonon band structure and density of states of FCC-Al.

Reproduces Fig. 1 from Togo & Tanaka, Scripta Materialia 108, 1-5 (2015).
"""

import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(os.path.realpath(sys.argv[0]))))
from phonon_model import (get_al_model, compute_bandstructure, compute_dos,
                           AL_A0, AL_MASS, generate_fcc_shells, fit_al_force_constants)

print("Fitting force constants...")
fc, err = fit_al_force_constants()
print(f"Fit error: {err:.4f}")
for i, (a, b) in enumerate(fc):
    print(f"  Shell {i+1}: alpha={a:.4f}, beta={b:.4f} N/m")

a = AL_A0
mass = AL_MASS
shells = generate_fcc_shells(a, n_shells=4)

# Compute band structure: L -> Gamma -> X -> W (same as paper)
print("Computing band structure...")
path_labels = ['L', 'G', 'X', 'W']
n_per_seg = 80
dist, freq, ticks, labels = compute_bandstructure(
    path_labels, n_per_seg, shells, fc, mass, a)

# Compute DOS
print("Computing DOS...")
freq_dos, dos = compute_dos(shells, fc, mass, a, n_grid=35, n_bins=300, sigma_THz=0.12)

# Save data
print("Saving data...")
outdir = os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(os.path.realpath(sys.argv[0])))), 'data')
os.makedirs(outdir, exist_ok=True)

# Band structure CSV
header = "distance,branch_1_THz,branch_2_THz,branch_3_THz"
data = np.column_stack([dist, freq])
np.savetxt(os.path.join(outdir, 'fig1_bandstructure.csv'), data,
           delimiter=',', header=header, fmt='%.8f')

# DOS CSV
header = "frequency_THz,dos_states_per_THz_per_atom"
data_dos = np.column_stack([freq_dos, dos])
np.savetxt(os.path.join(outdir, 'fig1_dos.csv'), data_dos,
           delimiter=',', header=header, fmt='%.8f')

# Plot
print("Plotting...")
plotdir = os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(os.path.realpath(sys.argv[0])))), 'plots')
os.makedirs(plotdir, exist_ok=True)

fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(8, 6),
                                gridspec_kw={'width_ratios': [3, 1], 'wspace': 0.05})

# Band structure
for branch in range(3):
    ax1.plot(dist, freq[:, branch], 'r-', linewidth=1.5)

for t in ticks:
    ax1.axvline(t, color='gray', linewidth=0.5, linestyle='-')

tick_labels_display = [r'L', r'$\Gamma$', r'X', r'W']
ax1.set_xticks(ticks)
ax1.set_xticklabels(tick_labels_display, fontsize=12)
ax1.set_ylabel('Frequency (THz)', fontsize=13)
ax1.set_xlabel('Wave vector', fontsize=13)
ax1.set_ylim(0, None)
ax1.set_xlim(dist[0], dist[-1])

# DOS
ax2.plot(dos, freq_dos, 'b-', linewidth=1.5)
ax2.set_xlabel('DOS\n(states/THz' + r'$\cdot$' + 'atom)', fontsize=11)
ax2.set_xlim(0, None)
ax2.set_ylim(0, None)
ax2.tick_params(axis='y', labelleft=False)

fig.suptitle('Phonon band structure and DOS of Al', fontsize=14, y=0.98)
plt.tight_layout()
plt.savefig(os.path.join(plotdir, 'fig1_bandstructure_dos.png'), dpi=150, bbox_inches='tight')
print("Done: fig1_bandstructure_dos.png")
