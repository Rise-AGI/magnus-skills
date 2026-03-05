"""
Fig 2: Dispersion relation comparison at theta_k = 85 deg
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from plasma_params import (PlasmaParams, omega0_dispersion,
                           omega_magnetosonic, omega_whistler,
                           solve_cold_dispersion_sweep)

params = PlasmaParams()
theta = 85.0 * np.pi / 180.0

k_range = np.linspace(5, 900, 300)
k_par = k_range * np.cos(theta)

w_analytic = omega0_dispersion(k_range, k_par, params)
w_ms = omega_magnetosonic(k_range, params)
w_wh = omega_whistler(k_range, k_par, params)

k_num = np.linspace(30, 900, 60)
w_num = solve_cold_dispersion_sweep(k_num, theta, params)

fig, ax = plt.subplots(figsize=(7, 5))

ax.plot(k_range, w_ms / params.omega_ci, 'r-.', linewidth=1.5,
        label='Magnetosonic: $\\omega = k v_A$')
ax.plot(k_range, w_wh / params.omega_ci, 'g--', linewidth=1.5,
        label=r'Whistler: $\omega = k k_\parallel v_A^2/\omega_{ci}$')
ax.plot(k_range, w_analytic / params.omega_ci, 'b-', linewidth=2,
        label=r'Analytical: $\omega_0 = k v_A\sqrt{1+k_\parallel^2 c^2/\omega_{pi}^2}$')
ax.plot(k_num, w_num / params.omega_ci, 'kD', markersize=5,
        markerfacecolor='none', markeredgewidth=1.2,
        label='Numerical (Stix)')

ax.set_xlabel(r'$k$ [m$^{-1}$]', fontsize=13)
ax.set_ylabel(r'$\omega_0 / \omega_{ci}$', fontsize=13)
ax.set_xlim(0, 900)
ax.set_ylim(0, 150)
ax.legend(fontsize=10, loc='upper left')
ax.set_title(r'Fig. 2: Dispersion relations ($\theta_k = 85°$)', fontsize=13)
ax.grid(True, alpha=0.3, linestyle='--')

plt.tight_layout()
plt.savefig('../plots/fig2_dispersion.png', dpi=300, bbox_inches='tight')
print("Saved ../plots/fig2_dispersion.png")
plt.close()

with open('../data/fig2_dispersion.csv', 'w') as f:
    f.write("# Figure 2: Dispersion relations at theta_k=85deg\n")
    f.write("# Columns: k [m^-1], omega_magnetosonic/omega_ci, omega_whistler/omega_ci, omega_analytic/omega_ci\n")
    f.write("k,omega_magnetosonic,omega_whistler,omega_analytic\n")
    for i in range(len(k_range)):
        f.write(f"{k_range[i]:.8f},{w_ms[i]/params.omega_ci:.8f},{w_wh[i]/params.omega_ci:.8f},{w_analytic[i]/params.omega_ci:.8f}\n")
    f.write("# Numerical (Stix) data points\n")
    f.write("# k_numerical [m^-1], omega_numerical/omega_ci\n")
    for i in range(len(k_num)):
        f.write(f"# num,{k_num[i]:.8f},{w_num[i]/params.omega_ci:.8f}\n")
print("Exported ../data/fig2_dispersion.csv")
