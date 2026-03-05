"""
Figure 2: Wilson loop area law demonstration.

Shows Wilson loops W(S x S) vs loop area S^2 on a semi-log plot,
demonstrating the area-law behavior W ~ exp(-sigma * S^2) that
signals quark confinement.

Related to Greensite (2003) Section 3.1 (Wilson loops and confinement).
"""
import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) if os.path.dirname(os.path.abspath(__file__)) else '.')
from su2_lattice import SU2Lattice

# Parameters
L = 8
n_therm = 20
n_meas = 5
seed = 42
beta_values = [2.0, 2.3, 2.6]
max_loop = min(L // 2, 4)

print("Figure 2: Wilson loop area law")
print(f"Lattice: {L}^4, thermalization: {n_therm}, measurements: {n_meas}")

all_data = {}

for beta in beta_values:
    print(f"  beta = {beta:.1f} ... ", flush=True)
    lattice = SU2Lattice(L, beta, seed=seed)
    lattice.initialize_cold()

    for _ in range(n_therm):
        lattice.sweep()

    wilson_measurements = {S: [] for S in range(1, max_loop + 1)}

    for m in range(n_meas):
        lattice.sweep()
        for S in range(1, max_loop + 1):
            w = lattice.square_wilson_loop(S)
            wilson_measurements[S].append(w)
            print(f"    meas {m+1}, W({S}x{S}) = {w:.8f}")

    all_data[beta] = {}
    for S in range(1, max_loop + 1):
        vals = np.array(wilson_measurements[S])
        all_data[beta][S] = {
            'mean': np.mean(vals),
            'std': np.std(vals) / np.sqrt(len(vals)) if len(vals) > 1 else 0.0
        }

# Save data
data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data')
os.makedirs(data_dir, exist_ok=True)

with open(os.path.join(data_dir, 'fig2_wilson_area_law.csv'), 'w') as f:
    f.write("# Figure 2: Wilson loops W(SxS) vs loop size\n")
    f.write("# L=8, n_therm=20, n_meas=5\n")
    f.write("beta,loop_size,area,wilson_mean,wilson_std\n")
    for beta in beta_values:
        for S in range(1, max_loop + 1):
            d = all_data[beta][S]
            f.write(f"{beta:.1f},{S},{S*S},{d['mean']:.8f},{d['std']:.8f}\n")

# Plot
plot_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'plots')
os.makedirs(plot_dir, exist_ok=True)

fig, ax = plt.subplots(figsize=(8, 6))
markers = ['o', 's', '^']
colors = ['blue', 'red', 'green']

for i, beta in enumerate(beta_values):
    areas = []
    means = []
    errs = []
    for S in range(1, max_loop + 1):
        d = all_data[beta][S]
        areas.append(S * S)
        means.append(d['mean'])
        errs.append(d['std'])
    ax.errorbar(areas, means, yerr=errs, fmt=markers[i]+'-', color=colors[i],
                label=fr'$\beta = {beta}$', markersize=6, capsize=3)

ax.set_yscale('log')
ax.set_xlabel(r'Loop area $S^2$', fontsize=14)
ax.set_ylabel(r'$W(S \times S)$', fontsize=14)
ax.set_title('Wilson Loop Area Law', fontsize=14)
ax.legend(fontsize=12)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(plot_dir, 'fig2_wilson_area_law.png'), dpi=150)
print("Saved fig2_wilson_area_law.png")
