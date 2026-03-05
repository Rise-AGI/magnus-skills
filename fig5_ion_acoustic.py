"""
Fig 5: Landau damping of the ion-acoustic (sound) mode
Hunana et al. (2019), Figure 5

Plots |Im(zeta)|/Re(zeta) vs tau = T_e/T_p for the ion-acoustic mode.
Compares exact kinetic solution with R_{4,3}, R_{4,2}, R_{5,3} closures.
Electron inertia is neglected in this figure.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, '.')
from plasma_dispersion import (R_exact, R_4_2, R_4_3, R_5_3,
                                solve_ion_acoustic_kinetic, solve_ion_acoustic_fluid)

# tau range (T_e/T_p)
tau_values = np.logspace(-0.5, 2.5, 100)

# Solve exact kinetic
zeta_kinetic = solve_ion_acoustic_kinetic(tau_values)
damping_kinetic = np.abs(np.imag(zeta_kinetic)) / np.abs(np.real(zeta_kinetic))

# Solve fluid models
closures = [
    (R_4_3, r'$R_{4,3}$ (Hammett-Perkins)', 'r--'),
    (R_4_2, r'$R_{4,2}$ (new static)', 'g-.'),
    (R_5_3, r'$R_{5,3}$ (new dynamic)', 'b:'),
]

fluid_results = {}
for R_func, label, style in closures:
    zeta_fluid = solve_ion_acoustic_fluid(tau_values, R_func)
    damping_fluid = np.abs(np.imag(zeta_fluid)) / np.abs(np.real(zeta_fluid))
    fluid_results[label] = (damping_fluid, style)

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

for ax_idx, (ax, xlim, title) in enumerate(zip(axes,
    [(0.3, 100), (0.3, 300)],
    ['Ion-acoustic damping', 'Ion-acoustic damping (extended)'])):

    ax.semilogy(tau_values, damping_kinetic, 'k-', linewidth=2, label='Exact kinetic')
    for label, (damping, style) in fluid_results.items():
        ax.semilogy(tau_values, damping, style, linewidth=1.5, label=label)

    ax.set_xlabel(r'$\tau = T_e/T_p$', fontsize=13)
    ax.set_ylabel(r'$|\mathrm{Im}(\zeta)| / \mathrm{Re}(\zeta)$', fontsize=13)
    ax.set_xlim(xlim)
    ax.set_ylim(1e-3, 2)
    ax.legend(fontsize=10)
    ax.set_title(f'Fig. 5: {title}', fontsize=12)
    ax.grid(True, alpha=0.3, which='both')

plt.tight_layout()
plt.savefig('../plots/fig5_ion_acoustic.png', dpi=200, bbox_inches='tight')
plt.close()
print("Saved ../plots/fig5_ion_acoustic.png")

# Export CSV
with open('../data/fig5_ion_acoustic.csv', 'w') as f:
    f.write("# Figure 5: Ion-acoustic mode Landau damping\n")
    f.write("# Columns: tau, damping_ratio_kinetic, damping_ratio_R4_3, damping_ratio_R4_2, damping_ratio_R5_3\n")
    f.write("tau,damping_kinetic,damping_R4_3,damping_R4_2,damping_R5_3\n")
    labels_ordered = [r'$R_{4,3}$ (Hammett-Perkins)', r'$R_{4,2}$ (new static)', r'$R_{5,3}$ (new dynamic)']
    for i in range(len(tau_values)):
        vals = [tau_values[i], damping_kinetic[i]]
        for label in labels_ordered:
            vals.append(fluid_results[label][0][i])
        f.write(",".join(f"{v:.8f}" for v in vals) + "\n")
print("Exported ../data/fig5_ion_acoustic.csv")
