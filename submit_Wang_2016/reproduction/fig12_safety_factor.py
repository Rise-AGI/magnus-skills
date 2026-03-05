"""
Figure 12: Energy limit and balance time vs safety factor q.

Parameters: E_l=0.2V/m, B0=2T, R0=1.7m
Initial: p_par=5m0c, p_perp=1m0c
Varied: q from 0.5 to 6.0
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) if os.path.dirname(os.path.abspath(__file__)) else '.')

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from runaway_physics import find_energy_limit

# Parameters
E_l = 0.2
B0 = 2.0
R0 = 1.7
r_orbit = 0.1
q_values = np.linspace(0.5, 6.0, 40)
t_max = 8.0

# Compute
E_max_list = []
t_blc_list = []
for q in q_values:
    E_max, t_blc = find_energy_limit(E_l, B0, R0, q, r_orbit, t_max=t_max)
    E_max_list.append(E_max)
    t_blc_list.append(t_blc)

E_max_arr = np.array(E_max_list)
t_blc_arr = np.array(t_blc_list)
print(f"q range: [{q_values[0]:.1f}, {q_values[-1]:.1f}]")
print(f"E_max range: [{E_max_arr.min():.1f}, {E_max_arr.max():.1f}] MeV")

# Save data
data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data')
os.makedirs(data_dir, exist_ok=True)

header = ("# Figure 12: Energy limit and balance time vs safety factor\n"
          "# Columns: q, E_max[MeV], t_blc[s]\n"
          "q,E_max,t_blc")
np.savetxt(os.path.join(data_dir, 'fig12_safety_factor.csv'),
           np.column_stack([q_values, E_max_arr, t_blc_arr]),
           delimiter=',', header=header, comments='', fmt='%.8e')

# Plot
fig, ax1 = plt.subplots(figsize=(8, 5))

color1 = 'blue'
ax1.plot(q_values, E_max_arr, color=color1, linewidth=2, label='$E_{max}$')
ax1.set_xlabel('Safety Factor $q$', fontsize=12)
ax1.set_ylabel('Maximum Energy $E_{max}$ [MeV]', fontsize=12, color=color1)
ax1.tick_params(axis='y', labelcolor=color1)

ax2 = ax1.twinx()
color2 = 'red'
ax2.plot(q_values, t_blc_arr, color=color2, linewidth=2, linestyle='--',
         label='$t_{blc}$')
ax2.set_ylabel('Balance Time $t_{blc}$ [s]', fontsize=12, color=color2)
ax2.tick_params(axis='y', labelcolor=color2)

fig.legend(loc='upper left', bbox_to_anchor=(0.15, 0.95), fontsize=11)
plt.title('Energy Limit and Balance Time vs Safety Factor', fontsize=12)
plt.tight_layout()

plot_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'plots')
os.makedirs(plot_dir, exist_ok=True)
plt.savefig(os.path.join(plot_dir, 'fig12_safety_factor.png'), dpi=150, bbox_inches='tight')
plt.close()

print("Saved: data/fig12_safety_factor.csv, plots/fig12_safety_factor.png")
