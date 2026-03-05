"""
Figure 2: Energy gap between ground and first excited state vs x.

Computes the gap in the Gauss law sector for both models at different N and d.
"""
import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(os.path.realpath(sys.argv[0]))))
from schwinger_ed import compute_ground_state

# Parameters
system_sizes = [4, 6]
d_values = [3, 5]
x_values = np.array([0.5, 1, 2, 4, 8, 12, 16, 20, 25])

results = {}

for model in ['cqed', 'zd']:
    results[model] = {}
    for d in d_values:
        results[model][d] = {}
        for N in system_sizes:
            gaps = []
            for x in x_values:
                try:
                    evals, _ = compute_ground_state(N, d, x, model, n_states=2)
                    if len(evals) >= 2:
                        gap = evals[1] - evals[0]
                    else:
                        gap = np.nan
                    gaps.append(gap)
                    print(f"  {model} d={d} N={N} x={x:.1f}: gap={gap:.8f}")
                except Exception as e:
                    gaps.append(np.nan)
                    print(f"  {model} d={d} N={N} x={x:.1f}: FAILED ({e})")
            results[model][d][N] = gaps

# Save data
os.makedirs('../data', exist_ok=True)

with open('../data/fig2_energy_gap.csv', 'w') as f:
    header = 'x'
    for model in ['cqed', 'zd']:
        for d in d_values:
            for N in system_sizes:
                header += f',gap_{model}_d{d}_N{N}'
    f.write(header + '\n')

    for i, x in enumerate(x_values):
        line = f'{x:.8f}'
        for model in ['cqed', 'zd']:
            for d in d_values:
                for N in system_sizes:
                    val = results[model][d][N][i]
                    line += f',{val:.8f}'
        f.write(line + '\n')

print("Data saved to ../data/fig2_energy_gap.csv")

# Plot
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

colors = {'cqed': {3: 'blue', 5: 'red'}, 'zd': {3: 'blue', 5: 'red'}}
markers_N = {4: 'v', 6: 'o'}

for model, ax, title in [('cqed', ax1, 'Truncated cQED'), ('zd', ax2, r'$Z_d$ model')]:
    for d in d_values:
        for N in system_sizes:
            gaps = results[model][d][N]
            ax.plot(x_values, gaps, marker=markers_N[N],
                    color=colors[model][d],
                    label=f'd={d}, N={N}',
                    markersize=4, linewidth=1)
    ax.set_xlabel('x')
    ax.set_ylabel(r'$|E_0 - E_1|$')
    ax.set_title(title)
    ax.legend(fontsize=8)

plt.tight_layout()
plt.savefig('../plots/fig2_energy_gap.png', dpi=150, bbox_inches='tight')
print("Plot saved to ../plots/fig2_energy_gap.png")
