"""
Figure 3: Quasi-harmonic approximation properties of FCC-Al.

Reproduces Fig. 3 from Togo & Tanaka, Scripta Materialia 108, 1-5 (2015).
Three panels: (a) phonon frequencies vs volume, (b) F(V,T), (c) thermal expansion.
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
                           compute_dos, phonon_free_energy_per_atom,
                           scale_force_constants, vinet_energy,
                           phonon_frequencies_Hz)
from scipy.constants import eV, N_A
from scipy.optimize import minimize_scalar

print("Fitting force constants...")
fc_base, _ = fit_al_force_constants()
a0 = AL_A0
mass = AL_MASS

# Volume range: 10 points around equilibrium
V0 = AL_V0  # m^3/atom
# Al: a=4.05 A -> V=16.61 A^3/atom. Range: ~15-18 A^3 = 15e-30 to 18e-30 m^3
volumes_A3 = np.linspace(59, 79, 10)  # conventional cell volume in A^3
volumes_per_atom = volumes_A3 * 1e-30 / 4  # 4 atoms per conventional cell

# Temperature array for panel (b)
T_panel = np.arange(0, 900, 100)  # 0 to 800 K

# Temperature array for panel (c)
T_fine = np.linspace(1, 800, 200)

print("Computing QHA properties at 10 volumes...")
n_V = len(volumes_per_atom)
n_T_panel = len(T_panel)
n_T_fine = len(T_fine)

freq_X = np.zeros((n_V, 3))
freq_L = np.zeros((n_V, 3))
F_ph_panel = np.zeros((n_V, n_T_panel))
F_ph_fine = np.zeros((n_V, n_T_fine))

FC_POWER = 13.0  # Force constant scaling power (p = 6*gamma, gamma~2.2 for Al)

for iv, V in enumerate(volumes_per_atom):
    a = (4 * V) ** (1.0 / 3.0)
    shells = generate_fcc_shells(a, n_shells=4)
    fc = scale_force_constants(fc_base, a, a0, power=FC_POWER)

    # Frequencies at X and L
    q_X = np.array([1, 0, 0], dtype=float) * 2 * np.pi / a
    q_L = np.array([0.5, 0.5, 0.5]) * 2 * np.pi / a
    freq_X[iv] = phonon_frequencies_Hz(q_X, shells, fc, mass) * 1e-12
    freq_L[iv] = phonon_frequencies_Hz(q_L, shells, fc, mass) * 1e-12

    # DOS for free energy
    freq_grid, dos = compute_dos(shells, fc, mass, a, n_grid=20, n_bins=150, sigma_THz=0.15)

    for it, T in enumerate(T_panel):
        F_ph_panel[iv, it] = phonon_free_energy_per_atom(freq_grid, dos, T)

    for it, T in enumerate(T_fine):
        F_ph_fine[iv, it] = phonon_free_energy_per_atom(freq_grid, dos, T)

    print(f"  Volume {iv+1}/{n_V}: a={a*1e10:.3f} A, "
          f"max_freq_X={np.max(freq_X[iv]):.2f} THz")

# Electronic energy (Vinet EOS)
E_el = np.array([vinet_energy(V, 0, V0, AL_B0, AL_B0P) for V in volumes_per_atom])

# Total free energy for panel (b)
F_total_panel = E_el[:, np.newaxis] + F_ph_panel  # eV/atom

# Total free energy for panel (c)
F_total_fine = E_el[:, np.newaxis] + F_ph_fine

# Compute thermal expansion using Gruneisen relation:
# beta = gamma * Cv / (B * V)
# This is the standard QHA result and avoids numerical differentiation issues
from phonon_model import thermal_properties as tp_func, compute_dos as cd_func
gamma_gr = FC_POWER / 6.0  # Mode Gruneisen parameter from scaling power

# Get Cv at equilibrium
shells_eq = generate_fcc_shells(a0, n_shells=4)
freq_eq, dos_eq = cd_func(shells_eq, fc_base, mass, a0, n_grid=25, n_bins=200, sigma_THz=0.15)
props_fine = tp_func(freq_eq, dos_eq, T_fine)
Cv_fine = props_fine['Cv']  # J/K/mol

# beta = gamma * Cv / (B * V * N_A)
# Cv is in J/K/mol, B in Pa, V in m^3/atom
beta_T = gamma_gr * Cv_fine / (AL_B0 * AL_V0 * N_A)

# Also compute V_eq from thermal expansion
V_eq = np.zeros(n_T_fine)
V_eq[0] = V0
for it in range(1, n_T_fine):
    dT_step = T_fine[it] - T_fine[it - 1]
    V_eq[it] = V_eq[it - 1] * (1 + beta_T[it] * dT_step)

# Experimental thermal expansion data (Wilson 1941, Nix & MacNair 1941)
T_exp = np.array([100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650])
beta_exp = np.array([29, 44, 56, 63, 69, 73, 76, 80, 84, 89, 95, 101]) * 1e-6  # K^-1

# Save data
outdir = os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(os.path.realpath(sys.argv[0])))), 'data')
os.makedirs(outdir, exist_ok=True)

# Panel (a) data
header_a = "volume_A3," + ",".join([f"freq_X_{i+1}_THz" for i in range(3)] +
                                    [f"freq_L_{i+1}_THz" for i in range(3)])
data_a = np.column_stack([volumes_A3, freq_X, freq_L])
np.savetxt(os.path.join(outdir, 'fig3a_freq_vs_volume.csv'), data_a,
           delimiter=',', header=header_a, fmt='%.8f')

# Panel (c) data
header_c = "T_K,V_eq_A3,beta_K_inv"
data_c = np.column_stack([T_fine, V_eq * 4e30, beta_T])
np.savetxt(os.path.join(outdir, 'fig3c_thermal_expansion.csv'), data_c,
           delimiter=',', header=header_c, fmt='%.8f')

# Plot
print("Plotting...")
plotdir = os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(os.path.realpath(sys.argv[0])))), 'plots')
os.makedirs(plotdir, exist_ok=True)

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(14, 5))

# Panel (a): Phonon frequencies vs volume
for i in range(3):
    style_X = 'ro-' if i == 0 else ('rs-' if i == 1 else 'r^-')
    style_L = 'bo-' if i == 0 else ('bs-' if i == 1 else 'b^-')
    label_X = 'X' if i == 0 else None
    label_L = 'L' if i == 0 else None
    ax1.plot(volumes_A3, freq_X[:, i], style_X, markersize=4, linewidth=1,
             label=label_X)
    ax1.plot(volumes_A3, freq_L[:, i], style_L, markersize=4, linewidth=1,
             label=label_L)

ax1.set_xlabel(r'Volume ($\AA^3$)', fontsize=12)
ax1.set_ylabel('Frequency (THz)', fontsize=12)
ax1.legend(fontsize=10)
ax1.set_title('(a)', fontsize=13, loc='left')

# Panel (b): F(V,T)
for it, T in enumerate(T_panel):
    # Offset for clarity (referenced to minimum)
    F_curve = F_total_panel[:, it]
    ax2.plot(volumes_A3, F_curve, 'b.-', markersize=6, linewidth=1)

    # Mark minimum
    idx_min = np.argmin(F_curve)
    if 0 < idx_min < n_V - 1:
        coeffs = np.polyfit(volumes_per_atom, F_curve, 4)
        p = np.poly1d(coeffs)
        res = minimize_scalar(p, bounds=(volumes_per_atom[1], volumes_per_atom[-2]),
                              method='bounded')
        ax2.plot(res.x * 4e30, p(res.x), 'rx', markersize=8, markeredgewidth=2)

ax2.set_xlabel(r'Volume ($\AA^3$)', fontsize=12)
ax2.set_ylabel(r'$U_{el} + F_{ph}$ (eV)', fontsize=12)
ax2.set_title('(b)', fontsize=13, loc='left')

# Add temperature labels
if n_T_panel > 0:
    ax2.text(volumes_A3[-1] + 0.3, F_total_panel[-1, 0], '0K', fontsize=8)
    if n_T_panel > 1:
        ax2.text(volumes_A3[-1] + 0.3, F_total_panel[-1, -1], f'{T_panel[-1]:.0f}K',
                 fontsize=8)

# Panel (c): Thermal expansion coefficient
ax3.plot(T_fine, beta_T * 1e6, 'b-', linewidth=2, label='Calculation')
ax3.plot(T_exp, beta_exp * 1e6, 'ro', markersize=6, label='Experiment')
ax3.set_xlabel('Temperature (K)', fontsize=12)
ax3.set_ylabel(r'Thermal expansion coeff. ($\times 10^{-6}$ K$^{-1}$)', fontsize=12)
ax3.set_title('(c)', fontsize=13, loc='left')
ax3.legend(fontsize=10)
ax3.set_xlim(0, 800)

plt.tight_layout()
plt.savefig(os.path.join(plotdir, 'fig3_qha_properties.png'), dpi=150, bbox_inches='tight')
print("Done: fig3_qha_properties.png")
