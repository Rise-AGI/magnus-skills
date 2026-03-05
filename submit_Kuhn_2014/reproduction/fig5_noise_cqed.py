"""
Figure 5: Effect of noise on adiabatic preparation - truncated cQED model.

Penalty energy per site and overlap vs noise strength lambda.
"""
import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(os.path.realpath(sys.argv[0]))))
from schwinger_ed import noisy_adiabatic_evolution

# Parameters - small systems only (full Hilbert space needed for noise)
N_values = [4]
d_values = [3, 5]
xF = 10.0
T = 100.0
n_steps = 300
lambda_values = [0, 1e-4, 5e-4, 1e-3, 3e-3]

results = {}

for N in N_values:
    for d in d_values:
        key = (N, d)
        overlaps = []
        errors = []
        penalties = []
        for lam in lambda_values:
            try:
                overlap, energy_error, penalty = noisy_adiabatic_evolution(
                    N, d, xF, T, n_steps, lam, 'cqed')
                if overlap is None:
                    raise ValueError("System too large")
                overlaps.append(overlap)
                errors.append(energy_error)
                penalties.append(penalty)
                print(f"  cQED N={N} d={d} lam={lam:.1e}: overlap={overlap:.6f}, P/N={penalty:.2e}")
            except Exception as e:
                overlaps.append(np.nan)
                errors.append(np.nan)
                penalties.append(np.nan)
                print(f"  cQED N={N} d={d} lam={lam:.1e}: FAILED ({e})")
        results[key] = {'overlaps': overlaps, 'errors': errors, 'penalties': penalties}

# Save data
os.makedirs('../data', exist_ok=True)

with open('../data/fig5_noise_cqed.csv', 'w') as f:
    header = 'lambda'
    for N in N_values:
        for d in d_values:
            header += f',overlap_N{N}_d{d},error_N{N}_d{d},penalty_N{N}_d{d}'
    f.write(header + '\n')

    for i, lam in enumerate(lambda_values):
        line = f'{lam:.8e}'
        for N in N_values:
            for d in d_values:
                key = (N, d)
                ov = results[key]['overlaps'][i]
                er = results[key]['errors'][i]
                pen = results[key]['penalties'][i]
                line += f',{ov:.8f},{er:.8e},{pen:.8e}'
        f.write(line + '\n')

print("Data saved to ../data/fig5_noise_cqed.csv")

# Plot
fig, ax = plt.subplots(1, 1, figsize=(8, 6))

colors = {3: 'blue', 5: 'red'}
for N in N_values:
    for d in d_values:
        key = (N, d)
        # Filter out lambda=0
        lam_nonzero = [l for l in lambda_values if l > 0]
        pen_nonzero = [results[key]['penalties'][i] for i, l in enumerate(lambda_values) if l > 0]
        if any(not np.isnan(p) for p in pen_nonzero):
            ax.plot(lam_nonzero, pen_nonzero,
                    marker='o', color=colors[d],
                    label=f'N={N}, d={d}', linewidth=1)

ax.set_xlabel(r'$\lambda$ (noise strength)')
ax.set_ylabel('P/N (penalty per site)')
ax.set_title('Truncated cQED: Effect of noise')
ax.legend()
ax.set_xscale('log')
ax.set_yscale('log')

# Inset for overlap
ax_inset = ax.inset_axes([0.55, 0.55, 0.4, 0.35])
for N in N_values:
    for d in d_values:
        key = (N, d)
        lam_nonzero = [l for l in lambda_values if l > 0]
        ov_nonzero = [results[key]['overlaps'][i] for i, l in enumerate(lambda_values) if l > 0]
        if any(not np.isnan(o) for o in ov_nonzero):
            ax_inset.plot(lam_nonzero, ov_nonzero,
                          marker='o', color=colors[d], markersize=4, linewidth=1)
ax_inset.set_xlabel(r'$\lambda$', fontsize=8)
ax_inset.set_ylabel('Overlap', fontsize=8)
ax_inset.set_xscale('log')
ax_inset.tick_params(labelsize=7)

plt.tight_layout()
plt.savefig('../plots/fig5_noise_cqed.png', dpi=150, bbox_inches='tight')
print("Plot saved to ../plots/fig5_noise_cqed.png")
