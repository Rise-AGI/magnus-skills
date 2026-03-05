"""
Fig 8: Real frequency of the Langmuir mode
Hunana et al. (2019), Figure 8

Plots omega_r/omega_pe vs k*lambda_D. Linear scale.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, '.')
from plasma_dispersion import (R_exact, R_4_2, R_4_3, R_5_3,
                                solve_langmuir_kinetic, solve_langmuir_fluid)

klD = np.linspace(0.1, 1.0, 80)

# Exact kinetic
zeta_kin = solve_langmuir_kinetic(klD)
# omega_r/omega_pe = Re(zeta)*k*v_th/omega_pe = Re(zeta)*k*lambda_D*sqrt(2)
# since lambda_D = v_th/(sqrt(2)*omega_pe)
omega_r_kin = np.real(zeta_kin) * klD * np.sqrt(2)

# Asymptotic: omega^2/omega_pe^2 = 1 + 3*klD^2
omega_r_asym = np.sqrt(1 + 3 * klD**2)

# Fluid closures
closures = [
    (R_4_3, r'$R_{4,3}$', 'r--'),
    (R_4_2, r'$R_{4,2}$', 'g-.'),
    (R_5_3, r'$R_{5,3}$', 'b:'),
]

fluid_results = {}
for R_func, label, style in closures:
    zeta_fl = solve_langmuir_fluid(klD, R_func)
    omega_r_fl = np.real(zeta_fl) * klD * np.sqrt(2)
    fluid_results[label] = (omega_r_fl, style)

fig, ax = plt.subplots(figsize=(8, 6))
ax.plot(klD, omega_r_kin, 'k-', linewidth=2, label='Exact kinetic')
ax.plot(klD, omega_r_asym, 'k:', linewidth=1.5, label='Asymptotic')
for label, (omega_r, style) in fluid_results.items():
    ax.plot(klD, omega_r, style, linewidth=1.5, label=label)

ax.set_xlabel(r'$k_\parallel \lambda_D$', fontsize=13)
ax.set_ylabel(r'$\omega_r / \omega_{pe}$', fontsize=13)
ax.legend(fontsize=11)
ax.set_title('Fig. 8: Real frequency of the Langmuir mode', fontsize=13)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('../plots/fig8_langmuir_freq.png', dpi=200, bbox_inches='tight')
plt.close()
print("Saved ../plots/fig8_langmuir_freq.png")

# Export CSV
with open('../data/fig8_langmuir_freq.csv', 'w') as f:
    f.write("# Figure 8: Langmuir mode real frequency\n")
    f.write("k_lambda_D,omega_r_kinetic,omega_r_asymptotic,omega_r_R4_3,omega_r_R4_2,omega_r_R5_3\n")
    labels_ordered = [r'$R_{4,3}$', r'$R_{4,2}$', r'$R_{5,3}$']
    for i in range(len(klD)):
        vals = [klD[i], omega_r_kin[i], omega_r_asym[i]]
        for label in labels_ordered:
            vals.append(fluid_results[label][0][i])
        f.write(",".join(f"{v:.8f}" for v in vals) + "\n")
print("Exported ../data/fig8_langmuir_freq.csv")
