"""
Figure 3: Dissipation scale r_diss/l as a function of gamma_max for different zeta.

Eq. 44: r_diss = 12 * (gamma_max / zeta)^2 * l

Also marks the astrophysical regimes (pulsars, AGN, GRBs).
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(os.path.realpath(sys.argv[0]))))
from ks_instability import dissipation_scale

data_dir = os.path.join(os.path.dirname(os.path.abspath(os.path.realpath(sys.argv[0]))), '..', 'data')
os.makedirs(data_dir, exist_ok=True)
plots_dir = os.path.join(os.path.dirname(os.path.abspath(os.path.realpath(sys.argv[0]))), '..', 'plots')
os.makedirs(plots_dir, exist_ok=True)

gamma_max = np.logspace(0.5, 6.5, 300)
zeta_values = [0.05, 0.1, 0.3, 0.5]

fig, ax = plt.subplots(1, 1, figsize=(8, 6))
colors = ['blue', 'red', 'green', 'purple']

all_r_diss = []
for zeta, color in zip(zeta_values, colors):
    l = 1.0  # normalized
    r_diss = dissipation_scale(gamma_max, zeta, l)
    r_diss_over_l = r_diss / l
    all_r_diss.append(r_diss_over_l)
    ax.loglog(gamma_max, r_diss_over_l, color=color, linewidth=2,
              label=rf'$\zeta = {zeta}$')

# Mark astrophysical regimes
ax.axvspan(5, 30, alpha=0.1, color='blue', label='AGN')
ax.axvspan(100, 3000, alpha=0.1, color='red', label='GRB')
ax.axvspan(1e4, 1e6, alpha=0.1, color='green', label='Pulsar')

# Save data
header = "gamma_max," + ",".join([f"r_diss_over_l_zeta_{z}" for z in zeta_values])
out = np.column_stack([gamma_max] + all_r_diss)
np.savetxt(os.path.join(data_dir, 'fig3_dissipation_scale.csv'), out,
           delimiter=',', header=header, comments='', fmt='%.8e')

ax.set_xlabel(r'$\gamma_{\mathrm{max}}$', fontsize=14)
ax.set_ylabel(r'$r_{\mathrm{diss}} / l$', fontsize=14)
ax.set_title('Dissipation Scale (Eq. 44)', fontsize=14)
ax.legend(fontsize=10, loc='upper left')
ax.grid(True, alpha=0.3, which='both')
ax.set_xlim(gamma_max[0], gamma_max[-1])
fig.tight_layout()
fig.savefig(os.path.join(plots_dir, 'fig3_dissipation_scale.png'), dpi=150)
plt.close(fig)
print("Figure 3 saved.")
