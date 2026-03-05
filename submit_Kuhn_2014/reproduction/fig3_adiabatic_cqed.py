"""
Figure 3: Adiabatic preparation of the ground state - truncated cQED model.

Overlap with exact ground state after adiabatic evolution as a function of
total evolution time T.
"""
import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(os.path.realpath(sys.argv[0]))))
from schwinger_ed import adiabatic_evolution

# Parameters - use small system sizes for ED
N_values = [4, 6]
d_values = [3, 5]
xF = 10.0  # final x value (paper uses ~100 but we use smaller for ED)
T_values = np.linspace(10, 100, 10)
n_steps = 200  # time steps

results = {}

for N in N_values:
    for d in d_values:
        key = (N, d)
        overlaps = []
        errors = []
        for T in T_values:
            try:
                overlap, energy_error = adiabatic_evolution(N, d, xF, T, n_steps, 'cqed')
                overlaps.append(overlap)
                errors.append(energy_error)
                print(f"  cQED N={N} d={d} T={T:.1f}: overlap={overlap:.6f}, eps={energy_error:.2e}")
            except Exception as e:
                overlaps.append(np.nan)
                errors.append(np.nan)
                print(f"  cQED N={N} d={d} T={T:.1f}: FAILED ({e})")
        results[key] = {'overlaps': overlaps, 'errors': errors}

# Save data
os.makedirs('../data', exist_ok=True)

with open('../data/fig3_adiabatic_cqed.csv', 'w') as f:
    header = 'T'
    for N in N_values:
        for d in d_values:
            header += f',overlap_N{N}_d{d},error_N{N}_d{d}'
    f.write(header + '\n')

    for i, T in enumerate(T_values):
        line = f'{T:.8f}'
        for N in N_values:
            for d in d_values:
                key = (N, d)
                ov = results[key]['overlaps'][i]
                er = results[key]['errors'][i]
                line += f',{ov:.8f},{er:.8e}'
        f.write(line + '\n')

print("Data saved to ../data/fig3_adiabatic_cqed.csv")

# Plot
fig, ax = plt.subplots(1, 1, figsize=(8, 6))

colors = {3: 'blue', 5: 'red'}
markers = {4: 'x', 6: '^'}

for N in N_values:
    for d in d_values:
        key = (N, d)
        ax.plot(T_values, results[key]['overlaps'],
                marker=markers[N], color=colors[d],
                label=f'N={N}, d={d}', markersize=6, linewidth=1)

ax.set_xlabel('T (total evolution time)')
ax.set_ylabel('Overlap')
ax.set_title('Truncated cQED: Adiabatic preparation overlap')
ax.legend()
ax.set_ylim(0.9, 1.01)
ax.axhline(y=0.99, color='gray', linestyle='--', alpha=0.5)

# Inset for energy error
ax_inset = ax.inset_axes([0.55, 0.1, 0.4, 0.35])
for N in N_values:
    for d in d_values:
        key = (N, d)
        ax_inset.plot(T_values, results[key]['errors'],
                      marker=markers[N], color=colors[d],
                      markersize=4, linewidth=1)
ax_inset.set_xlabel('T', fontsize=8)
ax_inset.set_ylabel(r'$\epsilon$', fontsize=8)
ax_inset.tick_params(labelsize=7)

plt.tight_layout()
plt.savefig('../plots/fig3_adiabatic_cqed.png', dpi=150, bbox_inches='tight')
print("Plot saved to ../plots/fig3_adiabatic_cqed.png")
