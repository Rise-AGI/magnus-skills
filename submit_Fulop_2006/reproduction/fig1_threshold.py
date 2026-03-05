"""
Fig 1: Stability threshold from Eq.(25)
T_crit = (Z^2 B_T / (20 n_r/n_e))^(2/3)
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

Z = 1
B_T = np.linspace(1, 4, 200)

ratios = [5e-4, 1e-3, 5e-3]
labels = [r'$n_r/n_e = 5\times10^{-4}$',
          r'$n_r/n_e = 10^{-3}$',
          r'$n_r/n_e = 5\times10^{-3}$']
styles = ['r-.', 'b--', 'k-']

fig, ax = plt.subplots(figsize=(7, 5))

all_T_crit = {}
for ratio, label, style in zip(ratios, labels, styles):
    T_crit = (Z**2 * B_T / (20 * ratio))**(2.0 / 3.0)
    all_T_crit[ratio] = T_crit
    ax.plot(B_T, T_crit, style, linewidth=2, label=label)

ax.text(2.5, 40, 'Unstable', fontsize=14, ha='center')
ax.text(2.5, 8, 'Stable', fontsize=14, ha='center')

ax.set_xlabel(r'$B$ [T]', fontsize=13)
ax.set_ylabel(r'$T_{\mathrm{eV}}$', fontsize=13)
ax.set_xlim(1, 4)
ax.set_ylim(0, 50)
ax.legend(fontsize=11, loc='upper left')
ax.set_title('Fig. 1: Stability threshold (Eq. 25)', fontsize=13)
ax.grid(True, alpha=0.3, linestyle='--')

plt.tight_layout()
plt.savefig('../plots/fig1_threshold.png', dpi=300, bbox_inches='tight')
print("Saved ../plots/fig1_threshold.png")
plt.close()

with open('../data/fig1_threshold.csv', 'w') as f:
    f.write("# Figure 1: Stability threshold from Eq.(25), Z=1\n")
    f.write("# Columns: B_T [T], T_crit_nr_5e-4 [eV], T_crit_nr_1e-3 [eV], T_crit_nr_5e-3 [eV]\n")
    f.write("B_T,T_crit_nr_5e-4,T_crit_nr_1e-3,T_crit_nr_5e-3\n")
    for i in range(len(B_T)):
        f.write(f"{B_T[i]:.8f},{all_T_crit[5e-4][i]:.8f},{all_T_crit[1e-3][i]:.8f},{all_T_crit[5e-3][i]:.8f}\n")
print("Exported ../data/fig1_threshold.csv")
