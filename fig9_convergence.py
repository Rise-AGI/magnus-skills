"""
Fig 9: Convergence of fluid and kinetic descriptions for the ion-acoustic mode
Hunana et al. (2019), Figure 9

Plots |Im(zeta)|/Re(zeta) vs tau for increasingly precise closures:
R_{4,2}, R_{5,3}, R_{6,4}, R_{7,5}
Demonstrates convergence to exact kinetic result.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, '.')
from plasma_dispersion import (R_exact, R_4_2, R_5_3, R_6_4, R_7_5,
                                solve_ion_acoustic_kinetic, solve_ion_acoustic_fluid)

tau_values = np.logspace(-0.5, 2.5, 100)

# Exact kinetic
zeta_kinetic = solve_ion_acoustic_kinetic(tau_values)
damping_kinetic = np.abs(np.imag(zeta_kinetic)) / np.abs(np.real(zeta_kinetic))

closures = [
    (R_4_2, r'$R_{4,2}$', 'g-'),
    (R_5_3, r'$R_{5,3}$', 'b--'),
    (R_6_4, r'$R_{6,4}$', 'orange'),
    (R_7_5, r'$R_{7,5}$', 'r-.'),
]

fluid_results = {}
for R_func, label, style in closures:
    zeta_fluid = solve_ion_acoustic_fluid(tau_values, R_func)
    damping_fluid = np.abs(np.imag(zeta_fluid)) / np.abs(np.real(zeta_fluid))
    fluid_results[label] = (damping_fluid, style)

fig, ax = plt.subplots(figsize=(8, 6))
ax.semilogy(tau_values, damping_kinetic, 'k-', linewidth=2, label='Exact kinetic')
for label, (damping, style) in fluid_results.items():
    if isinstance(style, str) and not style.startswith('#'):
        ax.semilogy(tau_values, damping, style, linewidth=1.5, label=label)
    else:
        ax.semilogy(tau_values, damping, color=style, linewidth=1.5, label=label)

ax.set_xlabel(r'$\tau = T_e/T_p$', fontsize=13)
ax.set_ylabel(r'$|\mathrm{Im}(\zeta)| / \mathrm{Re}(\zeta)$', fontsize=13)
ax.legend(fontsize=11)
ax.set_title('Fig. 9: Convergence of fluid and kinetic descriptions', fontsize=13)
ax.grid(True, alpha=0.3, which='both')
ax.set_xlim(0.3, 300)
ax.set_ylim(1e-3, 2)

plt.tight_layout()
plt.savefig('../plots/fig9_convergence.png', dpi=200, bbox_inches='tight')
plt.close()
print("Saved ../plots/fig9_convergence.png")

# Export CSV
with open('../data/fig9_convergence.csv', 'w') as f:
    f.write("# Figure 9: Convergence of fluid and kinetic descriptions\n")
    f.write("tau,damping_kinetic,damping_R4_2,damping_R5_3,damping_R6_4,damping_R7_5\n")
    labels_ordered = [r'$R_{4,2}$', r'$R_{5,3}$', r'$R_{6,4}$', r'$R_{7,5}$']
    for i in range(len(tau_values)):
        vals = [tau_values[i], damping_kinetic[i]]
        for label in labels_ordered:
            vals.append(fluid_results[label][0][i])
        f.write(",".join(f"{v:.8f}" for v in vals) + "\n")
print("Exported ../data/fig9_convergence.csv")
