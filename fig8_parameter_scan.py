"""
Figure 8: Parameter sensitivity of Im(n_eff) for HCPCF design optimization.

Reproduces the imaginary part of the effective refractive index as a function
of geometric parameters: (a) number of cladding rings, (b) pitch Lambda,
(c) core surround thickness t, (d) strut thickness w, (e) hole edge radius r.

Reference: Fig. 8 of Pomplun et al., phys. stat. sol. (a) 204, 3822 (2007)
Default parameters: Lambda=1550nm, r=300nm, w=50nm, t=170nm, 6 rings, lambda=589nm
"""
import sys
import os
import numpy as np

sys.path.insert(0, ".")
from fem_waveguide import pcf_loss_model, pcf_loss_vs_parameter

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ---- Default parameters ----
defaults = {
    'wavelength': 589.0,    # nm
    'pitch': 1550.0,        # nm
    'hole_radius': 300.0,   # nm
    'strut_width': 50.0,    # nm
    'n_rings': 6,
    'core_surround_thickness': 170.0  # nm
}

# ---- Parameter scans ----
# (a) Number of cladding rings
n_rings_vals = np.arange(1, 11)
im_neff_rings = pcf_loss_vs_parameter('n_rings', n_rings_vals, defaults)

# (b) Pitch Lambda
pitch_vals = np.linspace(1200, 1800, 200)
im_neff_pitch = pcf_loss_vs_parameter('pitch', pitch_vals, defaults)

# (c) Core surround thickness t
t_vals = np.linspace(50, 300, 200)
im_neff_t = pcf_loss_vs_parameter('core_surround', t_vals, defaults)

# (d) Strut thickness w
w_vals = np.linspace(10, 100, 200)
im_neff_w = pcf_loss_vs_parameter('strut_width', w_vals, defaults)

# (e) Hole edge radius r
r_vals = np.linspace(100, 500, 200)
im_neff_r = pcf_loss_vs_parameter('hole_radius', r_vals, defaults)

# ---- Export CSV ----
os.makedirs('../data', exist_ok=True)

with open('../data/fig8a_rings.csv', 'w') as f:
    f.write('# Im(n_eff) vs number of cladding rings\n')
    f.write('n_rings,im_neff\n')
    for i in range(len(n_rings_vals)):
        f.write(f'{int(n_rings_vals[i])},{im_neff_rings[i]:.8e}\n')

with open('../data/fig8b_pitch.csv', 'w') as f:
    f.write('# Im(n_eff) vs pitch Lambda [nm]\n')
    f.write('pitch_nm,im_neff\n')
    for i in range(len(pitch_vals)):
        f.write(f'{pitch_vals[i]:.4f},{im_neff_pitch[i]:.8e}\n')

with open('../data/fig8c_surround.csv', 'w') as f:
    f.write('# Im(n_eff) vs core surround thickness t [nm]\n')
    f.write('t_nm,im_neff\n')
    for i in range(len(t_vals)):
        f.write(f'{t_vals[i]:.4f},{im_neff_t[i]:.8e}\n')

with open('../data/fig8d_strut.csv', 'w') as f:
    f.write('# Im(n_eff) vs strut thickness w [nm]\n')
    f.write('w_nm,im_neff\n')
    for i in range(len(w_vals)):
        f.write(f'{w_vals[i]:.4f},{im_neff_w[i]:.8e}\n')

with open('../data/fig8e_radius.csv', 'w') as f:
    f.write('# Im(n_eff) vs hole edge radius r [nm]\n')
    f.write('r_nm,im_neff\n')
    for i in range(len(r_vals)):
        f.write(f'{r_vals[i]:.4f},{im_neff_r[i]:.8e}\n')

# ---- Plot ----
fig, axes = plt.subplots(2, 3, figsize=(16, 10))

# (a) Rings
ax = axes[0, 0]
ax.semilogy(n_rings_vals, im_neff_rings, 'bo-')
ax.set_xlabel('Number of cladding rings')
ax.set_ylabel('Im(n_eff)')
ax.set_title('(a) Cladding rings')
ax.grid(True, alpha=0.3)

# (b) Pitch
ax = axes[0, 1]
ax.semilogy(pitch_vals, im_neff_pitch, 'b-')
ax.set_xlabel(r'Pitch $\Lambda$ [nm]')
ax.set_ylabel('Im(n_eff)')
ax.set_title(r'(b) Pitch $\Lambda$')
ax.grid(True, alpha=0.3)

# (c) Core surround
ax = axes[0, 2]
ax.semilogy(t_vals, im_neff_t, 'b-')
ax.set_xlabel('Core surround thickness t [nm]')
ax.set_ylabel('Im(n_eff)')
ax.set_title('(c) Core surround t')
ax.grid(True, alpha=0.3)

# (d) Strut width
ax = axes[1, 0]
ax.semilogy(w_vals, im_neff_w, 'b-')
ax.set_xlabel('Strut thickness w [nm]')
ax.set_ylabel('Im(n_eff)')
ax.set_title('(d) Strut thickness w')
ax.grid(True, alpha=0.3)

# (e) Hole edge radius
ax = axes[1, 1]
ax.semilogy(r_vals, im_neff_r, 'b-')
ax.set_xlabel('Hole edge radius r [nm]')
ax.set_ylabel('Im(n_eff)')
ax.set_title('(e) Hole edge radius r')
ax.grid(True, alpha=0.3)

# Remove empty subplot
axes[1, 2].axis('off')

plt.suptitle('Im(n_eff) vs geometric parameters (HCPCF)', fontsize=14, y=1.02)
plt.tight_layout()
plt.savefig('../plots/fig8_parameter_scan.png', dpi=150, bbox_inches='tight')
plt.close()

print("Figure 8 complete: plots/fig8_parameter_scan.png, data/fig8[a-e]_*.csv")
