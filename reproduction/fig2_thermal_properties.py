"""
Figure 2: Thermal properties of FCC-Al.

Reproduces Fig. 2 from Togo & Tanaka, Scripta Materialia 108, 1-5 (2015).
Shows entropy, Cv, Cp, and Helmholtz free energy vs temperature.
"""

import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(os.path.realpath(sys.argv[0]))))
from phonon_model import (AL_A0, AL_MASS, AL_V0, AL_B0, AL_B0P,
                           generate_fcc_shells, fit_al_force_constants,
                           compute_dos, thermal_properties,
                           scale_force_constants, phonon_free_energy_per_atom,
                           vinet_energy)
from scipy.constants import k as kB, hbar, N_A, eV
from scipy.optimize import minimize_scalar

print("Fitting force constants...")
fc_base, _ = fit_al_force_constants()
a0 = AL_A0
mass = AL_MASS
shells = generate_fcc_shells(a0, n_shells=4)

# Compute DOS at equilibrium
print("Computing phonon DOS...")
freq_dos, dos = compute_dos(shells, fc_base, mass, a0, n_grid=35, n_bins=300, sigma_THz=0.12)

# Temperature array
T = np.linspace(0, 800, 801)

# Compute harmonic thermal properties (Cv, S, F)
print("Computing thermal properties...")
props = thermal_properties(freq_dos, dos, T)

# Compute Cp using QHA (simplified: Cp = Cv + T*V*beta^2*B)
# Use Gruneisen relation: Cp - Cv = T * V * beta^2 * B
# where beta = gamma * Cv / (V * B)
# So Cp = Cv * (1 + gamma * beta * T)
# Simplified: Cp = Cv + gamma^2 * Cv^2 * T / (V * B * N_A)

gamma_gr = 2.2  # Gruneisen parameter for Al
V_mol = AL_V0 * N_A  # m^3/mol
B = AL_B0  # Pa

# Cp - Cv = gamma^2 * Cv^2 * T / (V_mol * B)
Cp = np.zeros_like(T)
for i in range(len(T)):
    if T[i] > 0 and props['Cv'][i] > 0:
        Cp[i] = props['Cv'][i] + gamma_gr**2 * props['Cv'][i]**2 * T[i] / (V_mol * B)
    else:
        Cp[i] = props['Cv'][i]

# Experimental Cp data (JANAF tables, Chase 1998)
T_exp = np.array([100, 200, 298.15, 400, 500, 600, 700, 800])
Cp_exp = np.array([12.2, 21.3, 24.2, 25.7, 26.8, 28.2, 30.1, 31.7])

# Save data
outdir = os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(os.path.realpath(sys.argv[0])))), 'data')
os.makedirs(outdir, exist_ok=True)

header = "T_K,Cv_J_per_K_per_mol,S_J_per_K_per_mol,F_kJ_per_mol,Cp_J_per_K_per_mol"
data = np.column_stack([T, props['Cv'], props['S'], props['F'], Cp])
np.savetxt(os.path.join(outdir, 'fig2_thermal.csv'), data,
           delimiter=',', header=header, fmt='%.8f')

# Plot
print("Plotting...")
plotdir = os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(os.path.realpath(sys.argv[0])))),'plots')
os.makedirs(plotdir, exist_ok=True)

fig, ax = plt.subplots(figsize=(8, 6))

ax.plot(T, props['S'], 'b-', linewidth=2, label=r'Entropy (J/K$\cdot$mol)')
ax.plot(T, Cp, 'g-', linewidth=2, label=r'$C_P$ (J/K$\cdot$mol)')
ax.plot(T_exp, Cp_exp, 'k--', linewidth=1.5, label=r'$C_P$ exp. [Chase]')
ax.plot(T, props['Cv'], color='darkgreen', linewidth=2, label=r'$C_V$ (J/K$\cdot$mol)')
ax.plot(T, props['F'], 'r-', linewidth=2, label=r'Free energy (kJ/mol)')

ax.axhline(0, color='gray', linewidth=0.5, linestyle='--')
ax.set_xlabel('Temperature (K)', fontsize=13)
ax.set_ylabel('Thermal properties', fontsize=13)
ax.set_xlim(0, 800)
ax.set_ylim(-30, 60)
ax.legend(fontsize=10, loc='best')
ax.set_title('Thermal properties of Al', fontsize=14)

plt.tight_layout()
plt.savefig(os.path.join(plotdir, 'fig2_thermal_properties.png'), dpi=150, bbox_inches='tight')
print("Done: fig2_thermal_properties.png")
