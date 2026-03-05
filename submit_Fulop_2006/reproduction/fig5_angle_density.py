"""
Fig 5: Growth rate vs propagation angle for different densities
n_e = 2e19, 5e19, 1e20 m^-3; B = 2 T, T = 10 eV, n_r/n_e = 5e-3
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from plasma_params import (PlasmaParams, find_max_growth_rate,
                           growth_rate_numerical)

n_values = [2e19, 5e19, 1e20]
colors_an = ['r', 'g', 'b']
markers = ['^', 's', 'D']
labels_n = [r'$n_e = 2\times10^{19}$ m$^{-3}$',
            r'$n_e = 5\times10^{19}$ m$^{-3}$',
            r'$n_e = 10^{20}$ m$^{-3}$']

theta_range = np.linspace(1, 89.5, 80) * np.pi / 180.0

fig, ax = plt.subplots(figsize=(8, 5.5))

all_data = {}

for ne, col, mk, lab in zip(n_values, colors_an, markers, labels_n):
    params = PlasmaParams(n_e=ne)
    gamma_d_abs = abs(params.gamma_d)

    ratio_analytic = np.zeros(len(theta_range))
    for i, th in enumerate(theta_range):
        _, gamma_max = find_max_growth_rate(th, params, k_min=10, k_max=5000,
                                            n_scan=400, use_eq22=True)
        ratio_analytic[i] = gamma_max / gamma_d_abs - 1

    theta_num = np.linspace(75, 89, 12) * np.pi / 180.0
    ratio_num = np.zeros(len(theta_num))
    print(f"Computing numerical growth rates for n_e = {ne:.0e} ...")
    for i, th in enumerate(theta_num):
        k_scan = np.linspace(50, 5000, 100)
        gamma_vals = np.array([growth_rate_numerical(k, th, params) for k in k_scan])
        gamma_max = np.max(gamma_vals)
        ratio_num[i] = gamma_max / gamma_d_abs - 1
    print(f"  Done n_e = {ne:.0e}")

    all_data[ne] = {'ratio_analytic': ratio_analytic,
                    'theta_num': theta_num, 'ratio_num': ratio_num}

    ax.plot(np.degrees(theta_range), ratio_analytic, col + '--', linewidth=1.5,
            label=lab + ' (analytic)')
    ax.plot(np.degrees(theta_num), ratio_num, col + mk, markersize=6,
            markerfacecolor='none', markeredgewidth=1.2,
            label=lab + ' (numerical)')

ax.axhline(0, color='k', linewidth=0.8, linestyle='-')
ax.set_xlabel(r'$\theta_k$ [deg]', fontsize=13)
ax.set_ylabel(r'$\gamma_i / |\gamma_d| - 1$', fontsize=13)
ax.set_xlim(0, 90)
ax.set_ylim(-1.1, 1.0)
ax.legend(fontsize=9, ncol=2, loc='upper left')
ax.set_title(r'Fig. 5: Growth rate vs angle (different $n_e$)', fontsize=13)
ax.grid(True, alpha=0.3, linestyle='--')

plt.tight_layout()
plt.savefig('../plots/fig5_angle_density.png', dpi=300, bbox_inches='tight')
print("Saved ../plots/fig5_angle_density.png")
plt.close()

with open('../data/fig5_angle_density.csv', 'w') as f:
    f.write("# Figure 5: Normalized growth rate vs angle for different n_e\n")
    f.write("# Analytical: theta_deg, ratio_ne_2e19, ratio_ne_5e19, ratio_ne_1e20\n")
    f.write("theta_deg,ratio_analytic_ne_2e19,ratio_analytic_ne_5e19,ratio_analytic_ne_1e20\n")
    for i in range(len(theta_range)):
        f.write(f"{np.degrees(theta_range[i]):.8f}")
        for ne in n_values:
            f.write(f",{all_data[ne]['ratio_analytic'][i]:.8f}")
        f.write("\n")
    f.write("# Numerical: theta_deg, ratio_ne_2e19, ratio_ne_5e19, ratio_ne_1e20\n")
    for i in range(len(all_data[2e19]['theta_num'])):
        f.write(f"# num,{np.degrees(all_data[2e19]['theta_num'][i]):.8f}")
        for ne in n_values:
            f.write(f",{all_data[ne]['ratio_num'][i]:.8f}")
        f.write("\n")
print("Exported ../data/fig5_angle_density.csv")
