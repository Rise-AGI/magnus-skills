"""
Figure 4: Momentum evolution structure in (p_parallel, p_perpendicular) space.

Reproduces the relaxation-equation version of the momentum trajectory.
The full-orbit version (with neoclassical pitch-angle scattering) requires
~10^12 time steps and is not feasible here.

Default parameters: R0=1.7m, a=0.4m, q=2, B0=2T, E_l=0.2V/m
Initial: p_par0=5 m0c, p_perp0=1 m0c, R=1.8m
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) if os.path.dirname(os.path.abspath(__file__)) else '.')

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from runaway_physics import solve_momentum_evolution, M0C2_MEV

# Parameters (from paper)
R0 = 1.7       # major radius [m]
q = 2.0        # safety factor
B0 = 2.0       # toroidal field [T]
E_l = 0.2      # loop electric field [V/m]
p_par0 = 5.0   # initial parallel momentum [m0c]
p_perp0 = 1.0  # initial perpendicular momentum [m0c]
r_orbit = 0.1  # orbit minor radius [m]
t_max = 3.5    # simulation time [s]

# Solve
result = solve_momentum_evolution(p_par0, p_perp0, E_l, B0, R0, q,
                                  r_orbit, t_max, n_points=10000)

# Save data
data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data')
os.makedirs(data_dir, exist_ok=True)

header = ("# Figure 4: Momentum evolution (relaxation equations)\n"
          "# Columns: time[s], p_parallel[m0c], p_perpendicular[m0c], "
          "gamma, energy_MeV\n"
          "time,p_parallel,p_perpendicular,gamma,energy_MeV")

np.savetxt(os.path.join(data_dir, 'fig4_momentum.csv'),
           np.column_stack([result['t'], result['p_par'], result['p_perp'],
                            result['gamma'], result['energy_MeV']]),
           delimiter=',', header=header, comments='', fmt='%.8e')

# Plot
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# (a) p_par vs p_perp trajectory
ax = axes[0]
ax.plot(result['p_par'], result['p_perp'], 'b-', linewidth=0.5)
ax.set_xlabel(r'$p_\parallel / m_0 c$', fontsize=12)
ax.set_ylabel(r'$p_\perp / m_0 c$', fontsize=12)
ax.set_title('(a) Momentum trajectory (relaxation)', fontsize=12)
ax.set_xlim(0, None)
ax.set_ylim(0, None)

# (b) Energy vs time
ax = axes[1]
ax.plot(result['t'], result['energy_MeV'], 'r-', linewidth=1)
ax.set_xlabel('Time [s]', fontsize=12)
ax.set_ylabel('Kinetic Energy [MeV]', fontsize=12)
ax.set_title('(b) Energy evolution', fontsize=12)
ax.axhline(y=result['energy_MeV'][-1], color='k', linestyle='--',
           alpha=0.5, label=f"$E_{{max}} = {result['energy_MeV'][-1]:.1f}$ MeV")
ax.legend(fontsize=10)

plt.tight_layout()
plot_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'plots')
os.makedirs(plot_dir, exist_ok=True)
plt.savefig(os.path.join(plot_dir, 'fig4_momentum.png'), dpi=150, bbox_inches='tight')
plt.close()

print(f"Energy limit: {result['energy_MeV'][-1]:.2f} MeV")
print(f"Max p_parallel: {np.max(result['p_par']):.2f} m0c")
print(f"Max p_perp: {np.max(result['p_perp']):.4f} m0c")
print("Saved: data/fig4_momentum.csv, plots/fig4_momentum.png")
