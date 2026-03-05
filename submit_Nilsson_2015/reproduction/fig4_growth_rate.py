"""
Figure 4: Growth rate vs runaway electron density
Showing affine function: (1/(ntot-nr)) * dnr/dt = gamma_D + nr * gamma_A_bar
For E/Ec = 40 and 60, Te = 0.5 keV
"""
import sys
sys.path.insert(0, '.')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from physics import (dreicer_rate, avalanche_growth_rate,
                     thermal_collision_freq, solve_runaway_evolution)

n_e = 1e19
T_eV = 500.0
ln_Lambda = 15.0
nu_th = thermal_collision_freq(n_e, T_eV, ln_Lambda)

fig, ax = plt.subplots(figsize=(8, 6))

all_data = []
for E_over_Ec, color, ls in [(40, 'blue', '-'), (60, 'red', '-')]:
    gamma_D_total = dreicer_rate(E_over_Ec, n_e, T_eV, ln_Lambda)
    gamma_D_per_ne = gamma_D_total / n_e
    gamma_A_bar = avalanche_growth_rate(E_over_Ec, n_e, ln_Lambda)
    gamma_A_normalized = gamma_A_bar / n_e  # per ne per nr

    # Solve evolution to get trajectory
    t_max_tau = 300.0
    t, nr_frac = solve_runaway_evolution(E_over_Ec, T_eV, n_e, ln_Lambda,
                                          t_max_tau, n_points=5000,
                                          include_avalanche=True)
    nr = nr_frac * n_e
    ne_local = n_e - nr

    # Compute growth rate: (1/(ntot-nr)) * dnr/dt
    # = gamma_D_per_ne + nr * gamma_A_normalized
    growth = gamma_D_per_ne + nr * gamma_A_normalized

    # Normalize to nu_th
    growth_norm = growth / nu_th
    nr_norm = nr / n_e

    # Filter reasonable range
    mask = (nr_norm > 1e-15) & (nr_norm < 0.1)
    ax.plot(nr_norm[mask], growth_norm[mask], color=color, linewidth=2,
            label=f'With avalanche, $E/E_c={E_over_Ec}$')

    # Without avalanche (constant Dreicer)
    growth_no_av = np.full_like(nr_norm, gamma_D_per_ne / nu_th)
    ax.plot(nr_norm[mask], growth_no_av[mask], color=color, linewidth=1.5,
            linestyle='--', label=f'Dreicer only, $E/E_c={E_over_Ec}$')

# Save data
header = "nr_over_ne, growth_rate_Ec40_with_av, growth_rate_Ec40_dreicer, growth_rate_Ec60_with_av, growth_rate_Ec60_dreicer"

ax.set_xlabel(r'$n_r / n_e$', fontsize=14)
ax.set_ylabel(r'Growth rate $/ \nu_{th}$', fontsize=14)
ax.set_title(r'Figure 4: Growth rate vs $n_r$ ($T_e=0.5$ keV)', fontsize=14)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)
ax.set_xscale('log')

plt.tight_layout()
plt.savefig('../plots/fig4_growth_rate.png', dpi=150)
print("Figure 4 saved.")
