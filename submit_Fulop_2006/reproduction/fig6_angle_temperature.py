"""
Fig 6: Growth rate vs propagation angle for different temperatures
T = 5, 10, 15 eV; B = 2 T, n_e = 5e19, n_r/n_e = 5e-3
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from plasma_params import (PlasmaParams, find_max_growth_rate,
                           growth_rate_numerical)

T_values = [5, 10, 15]
colors_an = ['r', 'g', 'b']
markers = ['^', 's', 'D']
labels_T = [r'$T = 5$ eV', r'$T = 10$ eV', r'$T = 15$ eV']

theta_range = np.linspace(1, 89.5, 80) * np.pi / 180.0

fig, ax = plt.subplots(figsize=(8, 5.5))

all_data = {}

for T, col, mk, lab in zip(T_values, colors_an, markers, labels_T):
    params = PlasmaParams(T_eV=T)
    gamma_d_abs = abs(params.gamma_d)

    ratio_analytic = np.zeros(len(theta_range))
    for i, th in enumerate(theta_range):
        _, gamma_max = find_max_growth_rate(th, params, k_min=10, k_max=5000,
                                            n_scan=400, use_eq22=True)
        ratio_analytic[i] = gamma_max / gamma_d_abs - 1

    theta_num = np.linspace(75, 89, 12) * np.pi / 180.0
    ratio_num = np.zeros(len(theta_num))
    print(f"Computing numerical growth rates for T = {T} eV ...")
    for i, th in enumerate(theta_num):
        k_scan = np.linspace(50, 5000, 100)
        gamma_vals = np.array([growth_rate_numerical(k, th, params) for k in k_scan])
        gamma_max = np.max(gamma_vals)
        ratio_num[i] = gamma_max / gamma_d_abs - 1
    print(f"  Done T = {T} eV")

    all_data[T] = {'ratio_analytic': ratio_analytic,
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
ax.legend(fontsize=9, ncol=2, loc='upper left')
ax.set_title('Fig. 6: Growth rate vs angle (different T)', fontsize=13)
ax.grid(True, alpha=0.3, linestyle='--')

plt.tight_layout()
plt.savefig('../plots/fig6_angle_temperature.png', dpi=300, bbox_inches='tight')
print("Saved ../plots/fig6_angle_temperature.png")
plt.close()

with open('../data/fig6_angle_temperature.csv', 'w') as f:
    f.write("# Figure 6: Normalized growth rate vs angle for different T\n")
    f.write("# Analytical: theta_deg, ratio_T5, ratio_T10, ratio_T15\n")
    f.write("theta_deg,ratio_analytic_T5,ratio_analytic_T10,ratio_analytic_T15\n")
    for i in range(len(theta_range)):
        f.write(f"{np.degrees(theta_range[i]):.8f}")
        for T in T_values:
            f.write(f",{all_data[T]['ratio_analytic'][i]:.8f}")
        f.write("\n")
    f.write("# Numerical: theta_deg, ratio_T5, ratio_T10, ratio_T15\n")
    for i in range(len(all_data[5]['theta_num'])):
        f.write(f"# num,{np.degrees(all_data[5]['theta_num'][i]):.8f}")
        for T in T_values:
            f.write(f",{all_data[T]['ratio_num'][i]:.8f}")
        f.write("\n")
print("Exported ../data/fig6_angle_temperature.csv")
