"""
Fig 2: Asymmetry effects in bilayer graphene.
McCann, Phys. Rev. B 74, 161403(R) (2006)

(a) Delta(n) for different bare asymmetries Delta0
(b) Layer densities n1, n2 vs total density n
(c) Delta(0) vs bare asymmetry Delta0
(d) Cyclotron mass vs density for different Delta0
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(os.path.realpath(sys.argv[0]))))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from bilayer_graphene import (
    GAMMA1, GAMMA1_EV, VF, C0, e_charge, hbar, m_e,
    self_consistent_delta, compute_layer_densities,
    self_consistent_delta_at_zero, cyclotron_mass_array, pF_from_density
)

# ---- Parameters ----
delta0_values_gamma1 = [0.0, 0.1, 0.2]  # Delta0 in units of gamma1
delta0_J = [x * GAMMA1 for x in delta0_values_gamma1]
labels_a = [r'$\Delta_0 = 0$',
            r'$\Delta_0 = 0.1\gamma_1$',
            r'$\Delta_0 = 0.2\gamma_1$']
styles = ['k-', 'b--', 'r:']

# Density range: -100 to 100 x 10^11 cm^-2
n_cm2 = np.linspace(-100e11, 100e11, 120)
n_m2 = n_cm2 * 1e4
n_plot = n_cm2 / 1e11

# ---- Figure 2(a): Delta(n) for different Delta0 ----
print("Computing Fig 2(a): Delta(n) for different Delta0...")

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

delta_results = {}
for j, (d0, label, style) in enumerate(zip(delta0_J, labels_a, styles)):
    delta_arr = np.zeros(len(n_m2))
    for i, n in enumerate(n_m2):
        delta_arr[i] = self_consistent_delta(n, delta0=d0)
    delta_results[j] = delta_arr
    delta_meV = delta_arr / e_charge * 1000
    axes[0, 0].plot(n_plot, delta_meV, style, linewidth=2, label=label)
    print(f"  Delta0 = {delta0_values_gamma1[j]}*gamma1 done")

axes[0, 0].set_xlabel(r'$n$ [$10^{11}$ cm$^{-2}$]', fontsize=12)
axes[0, 0].set_ylabel(r'$\Delta$ [meV]', fontsize=12)
axes[0, 0].set_title('(a)', fontsize=13)
axes[0, 0].legend(fontsize=10)
axes[0, 0].grid(True, alpha=0.3)
axes[0, 0].set_xlim(-100, 100)

# ---- Figure 2(b): Layer densities for Delta0 = 0.2*gamma1 ----
print("Computing Fig 2(b): Layer densities...")

d0_for_b = 0.2 * GAMMA1
n1_arr = np.zeros(len(n_m2))
n2_arr = np.zeros(len(n_m2))
for i, n in enumerate(n_m2):
    delta_i = self_consistent_delta(n, delta0=d0_for_b)
    n1_arr[i], n2_arr[i] = compute_layer_densities(n, delta_i)

n1_plot = n1_arr / 1e4 / 1e11  # to 10^11 cm^-2
n2_plot = n2_arr / 1e4 / 1e11

axes[0, 1].plot(n_plot, n1_plot, 'b-', linewidth=2, label=r'$n_1$')
axes[0, 1].plot(n_plot, n2_plot, 'r--', linewidth=2, label=r'$n_2$')
axes[0, 1].set_xlabel(r'$n$ [$10^{11}$ cm$^{-2}$]', fontsize=12)
axes[0, 1].set_ylabel(r'$n_{1,2}$ [$10^{11}$ cm$^{-2}$]', fontsize=12)
axes[0, 1].set_title(r'(b) $\Delta_0 = 0.2\gamma_1$', fontsize=13)
axes[0, 1].legend(fontsize=10)
axes[0, 1].grid(True, alpha=0.3)
axes[0, 1].set_xlim(-100, 100)

# ---- Figure 2(c): Delta(0) vs Delta0 ----
print("Computing Fig 2(c): Delta(0) vs Delta0...")

delta0_meV_range = np.linspace(0, 100, 60)
delta0_J_range = delta0_meV_range * 1e-3 * e_charge

# eps_r = 1
delta_0_epsr1 = np.zeros(len(delta0_J_range))
for i, d0 in enumerate(delta0_J_range):
    delta_0_epsr1[i] = self_consistent_delta_at_zero(d0, eps_r=1.0)
delta_0_epsr1_meV = delta_0_epsr1 / e_charge * 1000

# eps_r = 2
delta_0_epsr2 = np.zeros(len(delta0_J_range))
for i, d0 in enumerate(delta0_J_range):
    delta_0_epsr2[i] = self_consistent_delta_at_zero(d0, eps_r=2.0)
delta_0_epsr2_meV = delta_0_epsr2 / e_charge * 1000

axes[1, 0].plot(delta0_meV_range, delta_0_epsr1_meV, 'k-', linewidth=2,
                label=r'$\varepsilon_r = 1$')
axes[1, 0].plot(delta0_meV_range, delta_0_epsr2_meV, 'b--', linewidth=2,
                label=r'$\varepsilon_r = 2$')
axes[1, 0].set_xlabel(r'$\Delta_0$ [meV]', fontsize=12)
axes[1, 0].set_ylabel(r'$\Delta(0)$ [meV]', fontsize=12)
axes[1, 0].set_title('(c)', fontsize=13)
axes[1, 0].legend(fontsize=10)
axes[1, 0].grid(True, alpha=0.3)

# ---- Figure 2(d): Cyclotron mass ----
print("Computing Fig 2(d): Cyclotron mass...")

delta0_for_d = [0.0, 0.1 * GAMMA1, 0.2 * GAMMA1]
labels_d = [r'$\Delta_0 = 0$', r'$\Delta_0 = 0.1\gamma_1$',
            r'$\Delta_0 = 0.2\gamma_1$']
styles_d = ['k-', 'b--', 'r:']

for d0, label, style in zip(delta0_for_d, labels_d, styles_d):
    delta_arr = np.zeros(len(n_m2))
    for i, n in enumerate(n_m2):
        delta_arr[i] = self_consistent_delta(n, delta0=d0)
    mc = cyclotron_mass_array(n_m2, delta_arr)
    axes[1, 1].plot(n_plot, mc, style, linewidth=2, label=label)

axes[1, 1].set_xlabel(r'$n$ [$10^{11}$ cm$^{-2}$]', fontsize=12)
axes[1, 1].set_ylabel(r'$m_c / m_e$', fontsize=12)
axes[1, 1].set_title('(d)', fontsize=13)
axes[1, 1].legend(fontsize=10, loc='upper right')
axes[1, 1].set_ylim(0, 0.15)
axes[1, 1].grid(True, alpha=0.3)
axes[1, 1].set_xlim(-100, 100)

plt.tight_layout()
plt.savefig('../plots/fig2_asymmetry.png', dpi=300, bbox_inches='tight')
print("Saved ../plots/fig2_asymmetry.png")
plt.close()

# ---- Export CSVs ----
# Fig 2(a)
with open('../data/fig2a_delta_n.csv', 'w') as f:
    f.write("# Fig 2(a): Delta(n) for different Delta0\n")
    f.write("# Columns: n_1e11_cm2, delta_D0_0_meV, delta_D0_01g1_meV, delta_D0_02g1_meV\n")
    f.write("n_1e11_cm2,delta_D0_0_meV,delta_D0_01g1_meV,delta_D0_02g1_meV\n")
    for i in range(len(n_m2)):
        vals = [delta_results[j][i] / e_charge * 1000 for j in range(3)]
        f.write(f"{n_plot[i]:.8f},{vals[0]:.8f},{vals[1]:.8f},{vals[2]:.8f}\n")
print("Exported ../data/fig2a_delta_n.csv")

# Fig 2(b)
with open('../data/fig2b_layer_densities.csv', 'w') as f:
    f.write("# Fig 2(b): Layer densities for Delta0 = 0.2*gamma1\n")
    f.write("# Columns: n_1e11_cm2, n1_1e11_cm2, n2_1e11_cm2\n")
    f.write("n_1e11_cm2,n1_1e11_cm2,n2_1e11_cm2\n")
    for i in range(len(n_m2)):
        f.write(f"{n_plot[i]:.8f},{n1_plot[i]:.8f},{n2_plot[i]:.8f}\n")
print("Exported ../data/fig2b_layer_densities.csv")

# Fig 2(c)
with open('../data/fig2c_screening.csv', 'w') as f:
    f.write("# Fig 2(c): Delta(0) vs Delta0\n")
    f.write("# Columns: delta0_meV, delta_0_epsr1_meV, delta_0_epsr2_meV\n")
    f.write("delta0_meV,delta_0_epsr1_meV,delta_0_epsr2_meV\n")
    for i in range(len(delta0_J_range)):
        f.write(f"{delta0_meV_range[i]:.8f},{delta_0_epsr1_meV[i]:.8f},{delta_0_epsr2_meV[i]:.8f}\n")
print("Exported ../data/fig2c_screening.csv")

# Fig 2(d) - compute fresh since we need the data arrays
with open('../data/fig2d_cyclotron_mass.csv', 'w') as f:
    f.write("# Fig 2(d): Cyclotron mass vs density\n")
    f.write("# Columns: n_1e11_cm2, mc_D0_0, mc_D0_01g1, mc_D0_02g1 (in units of m_e)\n")
    f.write("n_1e11_cm2,mc_D0_0,mc_D0_01g1,mc_D0_02g1\n")

    mc_all = {}
    for j, d0 in enumerate(delta0_for_d):
        delta_arr = np.zeros(len(n_m2))
        for i, n in enumerate(n_m2):
            delta_arr[i] = self_consistent_delta(n, delta0=d0)
        mc_all[j] = cyclotron_mass_array(n_m2, delta_arr)

    for i in range(len(n_m2)):
        f.write(f"{n_plot[i]:.8f},{mc_all[0][i]:.8f},{mc_all[1][i]:.8f},{mc_all[2][i]:.8f}\n")
print("Exported ../data/fig2d_cyclotron_mass.csv")
