"""
Fig 3: Growth rate vs wave number at theta_k = 85 deg
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from plasma_params import (PlasmaParams, growth_rate_eq22, growth_rate_eq20,
                           growth_rate_numerical)

params = PlasmaParams()
theta = 85.0 * np.pi / 180.0

k_range = np.linspace(20, 900, 200)

gamma_eq22 = np.array([growth_rate_eq22(k, theta, params) for k in k_range])
gamma_eq20 = np.array([growth_rate_eq20(k, theta, params) for k in k_range])

k_num = np.linspace(50, 900, 40)
print("Computing numerical growth rates...")
gamma_num = np.array([growth_rate_numerical(k, theta, params) for k in k_num])
print("Done.")

fig, ax = plt.subplots(figsize=(7, 5))

ax.plot(k_range, gamma_eq22 / params.omega_ci, 'k-', linewidth=2,
        label='Eq. (22)')
ax.plot(k_range, gamma_eq20 / params.omega_ci, 'r--', linewidth=2,
        label='Eq. (20)')
ax.plot(k_num, gamma_num / params.omega_ci, 'bD', markersize=5,
        markerfacecolor='none', markeredgewidth=1.2,
        label='Numerical')

ax.set_xlabel(r'$k$ [m$^{-1}$]', fontsize=13)
ax.set_ylabel(r'$\gamma_i / \omega_{ci}$', fontsize=13)
ax.set_xlim(0, 900)
ax.set_ylim(0, 1.8)
ax.legend(fontsize=11)
ax.set_title(r'Fig. 3: Growth rate ($\theta_k = 85°$)', fontsize=13)
ax.grid(True, alpha=0.3, linestyle='--')

plt.tight_layout()
plt.savefig('../plots/fig3_growth_rate.png', dpi=300, bbox_inches='tight')
print("Saved ../plots/fig3_growth_rate.png")
plt.close()

with open('../data/fig3_growth_rate.csv', 'w') as f:
    f.write("# Figure 3: Growth rate vs wave number at theta_k=85deg\n")
    f.write("# Analytical columns: k [m^-1], gamma_eq22/omega_ci, gamma_eq20/omega_ci\n")
    f.write("k,gamma_eq22,gamma_eq20\n")
    for i in range(len(k_range)):
        f.write(f"{k_range[i]:.8f},{gamma_eq22[i]/params.omega_ci:.8f},{gamma_eq20[i]/params.omega_ci:.8f}\n")
    f.write("# Numerical data points: k [m^-1], gamma_numerical/omega_ci\n")
    for i in range(len(k_num)):
        f.write(f"# num,{k_num[i]:.8f},{gamma_num[i]/params.omega_ci:.8f}\n")
print("Exported ../data/fig3_growth_rate.csv")
