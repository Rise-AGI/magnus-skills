"""
Figures 10-11: Kagome fiber attenuation spectra.

Reproduces the real part of n_eff and the imaginary part (attenuation)
as a function of wavelength for 19-cell and 1-cell kagome-structured fibers.

Reference: Fig. 10-11 of Pomplun et al., phys. stat. sol. (a) 204, 3822 (2007)
Parameters from Table 2:
  Layout A: 19-cell, pitch=10.9um, strut=0.51um
  Layout B: 1-cell, pitch=11.8um, strut=0.67um
"""
import sys
import os
import numpy as np

sys.path.insert(0, ".")
from fem_waveguide import kagome_attenuation_spectrum

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ---- Parameters from Table 2 ----
layouts = {
    'A': {'core_type': '19-cell', 'pitch': 10.9, 'strut_width': 0.51},  # um
    'B': {'core_type': '1-cell', 'pitch': 11.8, 'strut_width': 0.67},   # um
}

# Wavelength range: 500-1500 nm
wavelengths = np.linspace(500, 1500, 500)

# ---- Compute spectra ----
results = {}
for name, params in layouts.items():
    re_neff, im_neff, conf = kagome_attenuation_spectrum(
        wavelengths, params['pitch'], params['strut_width'], params['core_type']
    )
    results[name] = (re_neff, im_neff, conf)

# ---- Export CSV ----
os.makedirs('../data', exist_ok=True)

with open('../data/fig10_re_neff.csv', 'w') as f:
    f.write('# Real part of n_eff vs wavelength for kagome fibers\n')
    f.write('wavelength_nm,re_neff_19cell,re_neff_1cell\n')
    re_A, _, _ = results['A']
    re_B, _, _ = results['B']
    for i in range(len(wavelengths)):
        f.write(f'{wavelengths[i]:.4f},{re_A[i]:.8e},{re_B[i]:.8e}\n')

with open('../data/fig11_im_neff.csv', 'w') as f:
    f.write('# Imaginary part of n_eff vs wavelength for kagome fibers (attenuation)\n')
    f.write('wavelength_nm,im_neff_19cell,im_neff_1cell\n')
    _, im_A, _ = results['A']
    _, im_B, _ = results['B']
    for i in range(len(wavelengths)):
        f.write(f'{wavelengths[i]:.4f},{im_A[i]:.8e},{im_B[i]:.8e}\n')

with open('../data/fig13_confinement.csv', 'w') as f:
    f.write('# Core confinement ratio vs wavelength for kagome fibers\n')
    f.write('wavelength_nm,confinement_19cell,confinement_1cell\n')
    _, _, c_A = results['A']
    _, _, c_B = results['B']
    for i in range(len(wavelengths)):
        f.write(f'{wavelengths[i]:.4f},{c_A[i]:.8e},{c_B[i]:.8e}\n')

# ---- Plot Fig 10: Re(n_eff) ----
fig, ax = plt.subplots(1, 1, figsize=(10, 5))
ax.plot(wavelengths, results['A'][0], 'b-', label='19-cell (Layout A)')
ax.plot(wavelengths, results['B'][0], 'r-', label='1-cell (Layout B)')
ax.set_xlabel('Wavelength [nm]')
ax.set_ylabel('Re(n_eff)')
ax.set_title('Real part of effective refractive index')
ax.legend()
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('../plots/fig10_re_neff.png', dpi=150, bbox_inches='tight')
plt.close()

# ---- Plot Fig 11: Im(n_eff) (attenuation) ----
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

for ax, (name, label) in zip(axes, [('A', '19-cell'), ('B', '1-cell')]):
    _, im_neff, _ = results[name]
    ax.semilogy(wavelengths, im_neff, 'b-', linewidth=0.5)
    ax.set_xlabel('Wavelength [nm]')
    ax.set_ylabel('Im(n_eff)')
    ax.set_title(f'Attenuation spectrum ({label})')
    ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('../plots/fig11_attenuation.png', dpi=150, bbox_inches='tight')
plt.close()

# ---- Plot Fig 13: Confinement ----
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

for ax, (name, label) in zip(axes, [('A', '19-cell'), ('B', '1-cell')]):
    _, _, conf = results[name]
    ax.plot(wavelengths, conf, 'b-', linewidth=0.5)
    ax.set_xlabel('Wavelength [nm]')
    ax.set_ylabel('E_core / E_total')
    ax.set_title(f'Core confinement ({label})')
    ax.set_ylim([0, 1])
    ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('../plots/fig13_confinement.png', dpi=150, bbox_inches='tight')
plt.close()

print("Figures 10-11-13 complete.")
