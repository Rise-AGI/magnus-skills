"""
Figure 9: Energy limit and balance time vs loop electric field
for different magnetic field strengths.

Parameters: R0=1.7m, q=2, initial p_par=5m0c, p_perp=1m0c
Varied: E_l from 0.05 to 5 V/m, B0 = 1.5, 2.0, 3.0, 5.0 T
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
R0 = 1.7
q = 2.0
r_orbit = 0.1
B_values = [1.5, 2.0, 3.0, 5.0]
E_l_values = np.linspace(0.05, 5.0, 40)
colors = ['blue', 'red', 'green', 'purple']
t_max = 8.0

# Compute
results = {}
for B0 in B_values:
    E_max_list = []
    t_blc_list = []
    for E_l in E_l_values:
        E_max, t_blc = find_energy_limit(E_l, B0, R0, q, r_orbit, t_max=t_max)
        E_max_list.append(E_max)
        t_blc_list.append(t_blc)
    results[B0] = {'E_max': np.array(E_max_list), 't_blc': np.array(t_blc_list)}
    print(f"B0={B0}T: E_max range [{min(E_max_list):.1f}, {max(E_max_list):.1f}] MeV")

# Save data
data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data')
os.makedirs(data_dir, exist_ok=True)

cols = [E_l_values]
col_names = ['E_l']
for B0 in B_values:
    cols.append(results[B0]['E_max'])
    cols.append(results[B0]['t_blc'])
    col_names.extend([f'E_max_B{B0}', f't_blc_B{B0}'])

header = ("# Figure 9: Energy limit and balance time vs E_l for different B0\n"
          f"# Columns: {','.join(col_names)}")
np.savetxt(os.path.join(data_dir, 'fig9_energy_field.csv'),
           np.column_stack(cols), delimiter=',', header=header,
           comments='', fmt='%.8e')

# Plot
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# (a) Energy limit vs E_l
ax = axes[0]
for i, B0 in enumerate(B_values):
    ax.plot(E_l_values, results[B0]['E_max'], color=colors[i],
            linewidth=1.5, label=f'$B_0 = {B0}$ T')
ax.set_xlabel('Loop Electric Field $E_l$ [V/m]', fontsize=12)
ax.set_ylabel('Maximum Energy $E_{max}$ [MeV]', fontsize=12)
ax.set_title('(a) Energy Limit', fontsize=12)
ax.legend(fontsize=10)
ax.set_xlim(0, 5)
ax.set_ylim(0, None)

# (b) Balance time vs E_l
ax = axes[1]
for i, B0 in enumerate(B_values):
    ax.plot(E_l_values, results[B0]['t_blc'], color=colors[i],
            linewidth=1.5, label=f'$B_0 = {B0}$ T')
ax.set_xlabel('Loop Electric Field $E_l$ [V/m]', fontsize=12)
ax.set_ylabel('Energy Balance Time $t_{blc}$ [s]', fontsize=12)
ax.set_title('(b) Balance Time', fontsize=12)
ax.legend(fontsize=10)
ax.set_xlim(0, 5)
ax.set_ylim(0, None)

plt.tight_layout()
plot_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'plots')
os.makedirs(plot_dir, exist_ok=True)
plt.savefig(os.path.join(plot_dir, 'fig9_energy_field.png'), dpi=150, bbox_inches='tight')
plt.close()

print("Saved: data/fig9_energy_field.csv, plots/fig9_energy_field.png")
