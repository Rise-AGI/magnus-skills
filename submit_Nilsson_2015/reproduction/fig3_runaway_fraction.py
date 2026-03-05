"""
Figure 3: Fraction of runaway electrons vs time
With and without avalanche effect.
E/Ec = 30, Te = 0.5 keV
"""
import sys
sys.path.insert(0, '.')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from physics import solve_runaway_evolution

n_e = 1e19
T_eV = 500.0  # 0.5 keV
ln_Lambda = 15.0
E_over_Ec = 30.0
t_max = 400.0  # in tau_th units

# With avalanche
t1, nr1 = solve_runaway_evolution(E_over_Ec, T_eV, n_e, ln_Lambda,
                                   t_max, n_points=5000, include_avalanche=True)

# Without avalanche
t2, nr2 = solve_runaway_evolution(E_over_Ec, T_eV, n_e, ln_Lambda,
                                   t_max, n_points=5000, include_avalanche=False)

# Save data
header = "t_over_tau_th, nr_over_ntot_with_avalanche, nr_over_ntot_without_avalanche"
# Interpolate to common time grid
t_common = np.linspace(0, t_max, 2000)
nr1_interp = np.interp(t_common, t1, nr1)
nr2_interp = np.interp(t_common, t2, nr2)

np.savetxt('../data/fig3_runaway_fraction.csv',
           np.column_stack([t_common, nr1_interp, nr2_interp]),
           delimiter=',', header=header, comments='', fmt='%.8e')

# Plot
fig, ax = plt.subplots(figsize=(8, 6))
ax.semilogy(t1, nr1, 'b-', linewidth=2, label='With knock-on (avalanche)')
ax.semilogy(t2, nr2, 'r--', linewidth=2, label='Without knock-on (Dreicer only)')

ax.set_xlabel(r'$t / \tau_{th}$', fontsize=14)
ax.set_ylabel(r'$n_r / n_{tot}$', fontsize=14)
ax.set_title(r'Figure 3: Runaway fraction ($E/E_c=30$, $T_e=0.5$ keV)', fontsize=14)
ax.legend(fontsize=12)
ax.grid(True, alpha=0.3)
ax.set_ylim([1e-18, 1])
ax.set_xlim([0, t_max])

plt.tight_layout()
plt.savefig('../plots/fig3_runaway_fraction.png', dpi=150)
print("Figure 3 saved.")
