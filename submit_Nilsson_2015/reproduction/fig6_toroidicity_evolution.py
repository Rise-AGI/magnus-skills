"""
Figure 6: Evolution of runaway population for different inverse aspect ratios.
E/Ec = 40, Te = 0.5 keV
"""
import sys
sys.path.insert(0, '.')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from physics import solve_runaway_evolution

n_e = 1e19
T_eV = 500.0
ln_Lambda = 15.0
E_over_Ec = 40.0
t_max = 300.0

epsilons = [0.0, 0.05, 0.1, 0.2, 0.3, 0.5]
colors = ['black', 'blue', 'green', 'orange', 'red', 'purple']

fig, ax = plt.subplots(figsize=(8, 6))

all_t = []
all_nr = []

for eps, col in zip(epsilons, colors):
    t, nr_frac = solve_runaway_evolution(E_over_Ec, T_eV, n_e, ln_Lambda,
                                          t_max, n_points=3000,
                                          include_avalanche=True, epsilon=eps)
    ax.semilogy(t, nr_frac, color=col, linewidth=2,
                label=rf'$\epsilon = {eps}$')
    all_t.append(t)
    all_nr.append(nr_frac)

# Save data
t_common = np.linspace(0, t_max, 1000)
data_cols = [t_common]
header_parts = ["t_over_tau_th"]
for i, eps in enumerate(epsilons):
    nr_interp = np.interp(t_common, all_t[i], all_nr[i])
    data_cols.append(nr_interp)
    header_parts.append(f"nr_frac_eps_{eps}")

np.savetxt('../data/fig6_toroidicity_evolution.csv',
           np.column_stack(data_cols),
           delimiter=',', header=','.join(header_parts), comments='', fmt='%.8e')

ax.set_xlabel(r'$t / \tau_{th}$', fontsize=14)
ax.set_ylabel(r'$n_r / n_{tot}$', fontsize=14)
ax.set_title(r'Figure 6: Runaway evolution vs inverse aspect ratio ($E/E_c=40$, $T_e=0.5$ keV)',
             fontsize=12)
ax.legend(fontsize=11)
ax.grid(True, alpha=0.3)
ax.set_ylim([1e-18, 1])
ax.set_xlim([0, t_max])

plt.tight_layout()
plt.savefig('../plots/fig6_toroidicity_evolution.png', dpi=150)
print("Figure 6 saved.")
