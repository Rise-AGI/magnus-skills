"""
Fig 1: Band structure and asymmetry gap Delta(n).
McCann, Phys. Rev. B 74, 161403(R) (2006)

(a) Band structure near K point for large asymmetry Delta = gamma1
(b) Self-consistent Delta_tilde(n) (solid) vs analytic Eq.(1) (dashed)
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(os.path.realpath(sys.argv[0]))))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from bilayer_graphene import (
    GAMMA1, GAMMA1_EV, VF, C0, e_charge, hbar,
    band_energies, self_consistent_delta, delta_analytic_eq1
)

# ---- Figure 1(a): Band structure ----
print("Computing Fig 1(a): Band structure...")

# Use Delta = gamma1 for illustration (as stated in caption)
delta_large = GAMMA1

# Momentum range: p from -p_max to p_max
p_max = 2.5 * GAMMA1 / VF  # covers well beyond relevant scale
Np = 500
p_arr = np.linspace(0, p_max, Np)

# Compute bands
E = band_energies(p_arr, delta_large) / e_charge  # convert to eV

# Also compute for Delta=0 reference
E0 = band_energies(p_arr, 0.0) / e_charge

# Convert momentum to dimensionless vp/gamma1
vp_norm = VF * p_arr / GAMMA1

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Plot bands with Delta = gamma1
for i in range(4):
    ax1.plot(vp_norm, E[:, i], 'b-', linewidth=1.5)
    ax1.plot(-vp_norm, E[:, i], 'b-', linewidth=1.5)

ax1.axhline(0, color='gray', linewidth=0.5, linestyle='--')
ax1.set_xlabel(r'$vp / \gamma_1$', fontsize=13)
ax1.set_ylabel(r'$\epsilon$ [eV]', fontsize=13)
ax1.set_title(r'(a) Band structure, $\Delta = \gamma_1$', fontsize=13)
ax1.set_xlim(-2, 2)
ax1.set_ylim(-0.8, 0.8)
ax1.grid(True, alpha=0.3)

# ---- Figure 1(b): Delta(n) ----
print("Computing Fig 1(b): Delta(n)...")

# Density range: 0 to 100 x 10^11 cm^-2 = 0 to 10^17 m^-2
n_cm2 = np.linspace(0.5e11, 100e11, 80)  # in cm^-2
n_m2 = n_cm2 * 1e4  # convert to m^-2

# Numerical self-consistent Delta(n) with Delta0 = 0
delta_numerical = np.zeros(len(n_m2))
for i, n in enumerate(n_m2):
    delta_numerical[i] = self_consistent_delta(n, delta0=0.0)
    if (i + 1) % 20 == 0:
        print(f"  {i+1}/{len(n_m2)} densities computed")

# Convert to Delta_tilde = |Delta| * gamma1 / sqrt(gamma1^2 + Delta^2)
delta_tilde_num = np.abs(delta_numerical) * GAMMA1 / np.sqrt(GAMMA1**2 + delta_numerical**2)
delta_tilde_num_meV = delta_tilde_num / e_charge * 1000  # to meV

# Analytic Eq.(1)
delta_analytic = np.array([delta_analytic_eq1(n) for n in n_m2])
delta_tilde_ana = np.abs(delta_analytic) * GAMMA1 / np.sqrt(GAMMA1**2 + delta_analytic**2)
delta_tilde_ana_meV = delta_tilde_ana / e_charge * 1000

n_plot = n_cm2 / 1e11  # in units of 10^11 cm^-2

ax2.plot(n_plot, delta_tilde_num_meV, 'k-', linewidth=2, label='Numerical')
ax2.plot(n_plot, delta_tilde_ana_meV, 'k--', linewidth=2, label='Eq. (1)')
ax2.set_xlabel(r'$n$ [$10^{11}$ cm$^{-2}$]', fontsize=13)
ax2.set_ylabel(r'$\tilde{\Delta}$ [meV]', fontsize=13)
ax2.set_title(r'(b) $\tilde{\Delta}(n)$, $\Delta_0 = 0$', fontsize=13)
ax2.legend(fontsize=11)
ax2.set_xlim(0, 100)
ax2.set_ylim(0, max(delta_tilde_num_meV) * 1.1)
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('../plots/fig1_bands_gap.png', dpi=300, bbox_inches='tight')
print("Saved ../plots/fig1_bands_gap.png")
plt.close()

# ---- Export CSV ----
# Fig 1(a) band structure
with open('../data/fig1a_bandstructure.csv', 'w') as f:
    f.write("# Fig 1(a): Band structure with Delta = gamma1\n")
    f.write("# Columns: vp_over_gamma1, E_band1_eV, E_band2_eV, E_band3_eV, E_band4_eV\n")
    f.write("vp_over_gamma1,E_band1_eV,E_band2_eV,E_band3_eV,E_band4_eV\n")
    for i in range(len(p_arr)):
        f.write(f"{vp_norm[i]:.8f},{E[i,0]:.8f},{E[i,1]:.8f},{E[i,2]:.8f},{E[i,3]:.8f}\n")
print("Exported ../data/fig1a_bandstructure.csv")

# Fig 1(b) gap vs density
with open('../data/fig1b_delta_n.csv', 'w') as f:
    f.write("# Fig 1(b): Self-consistent Delta_tilde(n) with Delta0=0\n")
    f.write("# Columns: n_1e11_cm2, delta_tilde_numerical_meV, delta_tilde_analytic_meV\n")
    f.write("n_1e11_cm2,delta_tilde_numerical_meV,delta_tilde_analytic_meV\n")
    for i in range(len(n_m2)):
        f.write(f"{n_plot[i]:.8f},{delta_tilde_num_meV[i]:.8f},{delta_tilde_ana_meV[i]:.8f}\n")
print("Exported ../data/fig1b_delta_n.csv")
