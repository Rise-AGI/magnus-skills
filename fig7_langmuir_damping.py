"""
Fig 7: Landau damping of the Langmuir mode
Hunana et al. (2019), Figure 7

Plots |Im(zeta)|/Re(zeta) vs k*lambda_D for the Langmuir mode.
Log-log scale. Compares exact kinetic with several Pade closures.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, '.')
from plasma_dispersion import (R_exact, R_4_2, R_4_3, R_5_3,
                                solve_langmuir_kinetic, solve_langmuir_fluid)

# k*lambda_D range
klD = np.logspace(-1, 0.1, 80)

# Exact kinetic
zeta_kin = solve_langmuir_kinetic(klD)
damping_kin = np.abs(np.imag(zeta_kin)) / np.abs(np.real(zeta_kin))

# Asymptotic kinetic formula (Eq. 459)
# gamma/omega_r ~ sqrt(pi/8) / (klD)^3 * exp(-(1+3*klD^2)/(2*klD^2))
zeta_r = np.sqrt((1 + 3 * klD**2) / (2 * klD**2))
damping_asym = np.sqrt(np.pi / 8) / klD**3 * np.exp(-(1 + 3 * klD**2) / (2 * klD**2))

# Fluid closures
closures = [
    (R_4_3, r'$R_{4,3}$', 'r--'),
    (R_4_2, r'$R_{4,2}$', 'g-.'),
    (R_5_3, r'$R_{5,3}$', 'b:'),
]

fluid_results = {}
for R_func, label, style in closures:
    zeta_fl = solve_langmuir_fluid(klD, R_func)
    damping_fl = np.abs(np.imag(zeta_fl)) / np.abs(np.real(zeta_fl))
    fluid_results[label] = (damping_fl, style)

fig, ax = plt.subplots(figsize=(8, 6))
ax.loglog(klD, damping_kin, 'k-', linewidth=2, label='Exact kinetic')
ax.loglog(klD, damping_asym, 'k:', linewidth=1.5, label='Asymptotic (Eq. 459)')
for label, (damping, style) in fluid_results.items():
    ax.loglog(klD, damping, style, linewidth=1.5, label=label)

ax.set_xlabel(r'$k_\parallel \lambda_D$', fontsize=13)
ax.set_ylabel(r'$|\gamma| / \omega_r$', fontsize=13)
ax.legend(fontsize=11)
ax.set_title('Fig. 7: Landau damping of the Langmuir mode', fontsize=13)
ax.grid(True, alpha=0.3, which='both')
ax.set_xlim(0.1, 1.2)
ax.set_ylim(1e-8, 2)

plt.tight_layout()
plt.savefig('../plots/fig7_langmuir_damping.png', dpi=200, bbox_inches='tight')
plt.close()
print("Saved ../plots/fig7_langmuir_damping.png")

# Export CSV
with open('../data/fig7_langmuir_damping.csv', 'w') as f:
    f.write("# Figure 7: Langmuir mode Landau damping\n")
    f.write("k_lambda_D,damping_kinetic,damping_asymptotic,damping_R4_3,damping_R4_2,damping_R5_3\n")
    labels_ordered = [r'$R_{4,3}$', r'$R_{4,2}$', r'$R_{5,3}$']
    for i in range(len(klD)):
        vals = [klD[i], damping_kin[i], damping_asym[i]]
        for label in labels_ordered:
            vals.append(fluid_results[label][0][i])
        f.write(",".join(f"{v:.8f}" for v in vals) + "\n")
print("Exported ../data/fig7_langmuir_damping.csv")
