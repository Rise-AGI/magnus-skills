"""
Figure 1: Energy density in the thermodynamic limit for various x values.

Compares truncated cQED and Zd models at different d values.
Uses exact diagonalization for small systems (N=4,6,8) and
extrapolates to the thermodynamic limit.
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
system_sizes = [4, 6, 8]
d_values = [3, 5, 7]
x_values = [1, 4, 9, 16, 25]  # Smaller x values since ED is limited

# Exact Schwinger model value in massless continuum: -1/pi
exact_continuum = -1.0 / np.pi

results = {}

for model in ['cqed', 'zd']:
    results[model] = {}
    for d in d_values:
        results[model][d] = {}
        for x in x_values:
            energies = []
            sizes_used = []
            for N in system_sizes:
                try:
                    evals, _ = compute_ground_state(N, d, x, model, n_states=1)
                    omega = evals[0] / N
                    energies.append(omega)
                    sizes_used.append(N)
                    print(f"  {model} d={d} x={x} N={N}: omega={omega:.8f}")
                except Exception as e:
                    print(f"  {model} d={d} x={x} N={N}: FAILED ({e})")
            if len(energies) >= 2:
                # Extrapolate to thermodynamic limit
                inv_N = 1.0 / np.array(sizes_used)
                coeffs = np.polyfit(inv_N, energies, 1)
                omega_inf = coeffs[1]
            elif len(energies) == 1:
                omega_inf = energies[0]
            else:
                omega_inf = np.nan
            results[model][d][x] = omega_inf

# Save data
os.makedirs('../data', exist_ok=True)

ga_values = [1.0 / np.sqrt(x) for x in x_values]

with open('../data/fig1_energy_density.csv', 'w') as f:
    header = 'ga'
    for model in ['cqed', 'zd']:
        for d in d_values:
            header += f',omega_{model}_d{d}'
    f.write(header + '\n')

    for i, x in enumerate(x_values):
        ga = ga_values[i]
        line = f'{ga:.8f}'
        for model in ['cqed', 'zd']:
            for d in d_values:
                val = results[model][d].get(x, np.nan)
                line += f',{val:.8f}'
        f.write(line + '\n')

print("Data saved to ../data/fig1_energy_density.csv")

# Plot
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

colors_d = {3: 'blue', 5: 'green', 7: 'red'}
markers_model = {'cqed': 'x', 'zd': 'o'}

for model, ax, title in [('cqed', ax1, 'Truncated cQED'), ('zd', ax2, r'$Z_d$ model')]:
    for d in d_values:
        omegas = [results[model][d].get(x, np.nan) for x in x_values]
        ax.plot(ga_values, omegas, marker=markers_model[model],
                color=colors_d[d], label=f'd={d}', linewidth=1)
    ax.axhline(y=exact_continuum, color='gray', linestyle='--',
               label=r'Schwinger ($-1/\pi$)', alpha=0.7)
    ax.set_xlabel(r'$ga = 1/\sqrt{x}$')
    ax.set_ylabel(r'$\omega = E_0/N$')
    ax.set_title(title)
    ax.legend()

plt.tight_layout()
plt.savefig('../plots/fig1_energy_density.png', dpi=150, bbox_inches='tight')
print("Plot saved to ../plots/fig1_energy_density.png")
