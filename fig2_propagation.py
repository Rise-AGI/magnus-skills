"""
Figure 2: Propagation intensity maps of HGB on curved surface.

Reproduces the flattened 2D intensity maps for orders n=0,1,4,10 on CGCS,
showing two full periods, for two different divergence coefficient values.
Also generates flat-space comparison.
"""
import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) if os.path.dirname(os.path.abspath(__file__)) else '.')
from hgb_core import abcd_matrix, compute_intensity_map, compute_intensity_map_flat, divergence_coefficient

# Parameters
r0 = 1.0
r_over_r0 = 2.0
r = r_over_r0 * r0

# Wave number and beam parameters
k = 100.0  # wave number

# Two different gamma_d values: one > 1 (diverge first), one < 1 (converge first)
# gamma_d = 2r / (k * sigma^2)
# For gamma_d = 0.5: sigma^2 = 2r/(k*0.5) = 4r/k
# For gamma_d = 2.0: sigma^2 = 2r/(k*2.0) = r/k
gamma_d_values = [0.5, 2.0]
sigma_values = [np.sqrt(2*r / (k * gd)) for gd in gamma_d_values]

orders = [0, 1, 4, 10]
period = np.pi * r / r0
num_theta = 300
num_h = 200

# Two full periods
theta_array = np.linspace(0.01, 2 * period - 0.01, num_theta)
h_max = 0.3  # transverse range
h_array = np.linspace(-h_max, h_max, num_h)

# Save CSV data
data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data')
os.makedirs(data_dir, exist_ok=True)
plot_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'plots')
os.makedirs(plot_dir, exist_ok=True)

# Figure 2a: Curved surface propagation
fig, axes = plt.subplots(len(orders), len(gamma_d_values), figsize=(12, 14))
fig.suptitle('Fig. 2(a): HGB Propagation on Curved Surface (Two Periods)', fontsize=14)

all_data = {}

for j, gamma_d in enumerate(gamma_d_values):
    sigma = sigma_values[j]
    for i, n in enumerate(orders):
        print(f"Computing curved surface: gamma_d={gamma_d}, n={n}")
        I_map = compute_intensity_map(n, sigma, k, r0, r, theta_array, h_array)

        # Normalize each map to its max for visualization
        I_max = I_map.max()
        if I_max > 0:
            I_plot = I_map / I_max
        else:
            I_plot = I_map

        ax = axes[i, j]
        extent = [theta_array[0] / period, theta_array[-1] / period, h_array[0], h_array[-1]]
        im = ax.imshow(I_plot.T, extent=extent, aspect='auto', origin='lower',
                       cmap='hot', vmin=0, vmax=1)
        ax.set_ylabel('h')
        if i == 0:
            ax.set_title(f'$\\gamma_d = {gamma_d}$')
        if i == len(orders) - 1:
            ax.set_xlabel('$\\theta$ / period')
        ax.text(0.02, 0.95, f'n={n}', transform=ax.transAxes,
                color='white', fontsize=10, va='top')

        key = f"curved_gd{gamma_d}_n{n}"
        all_data[key] = I_map

plt.tight_layout()
plt.savefig(os.path.join(plot_dir, 'fig2a_curved_propagation.png'), dpi=150, bbox_inches='tight')
plt.close()

# Figure 2b: Flat space comparison
fig2, axes2 = plt.subplots(len(orders), 1, figsize=(8, 14))
fig2.suptitle('Fig. 2(b): HGB Propagation in Flat Space', fontsize=14)

sigma_flat = sigma_values[1]  # Use gamma_d=2 sigma for comparison
z_max = 2 * period * r0  # equivalent propagation distance
z_array = np.linspace(0.01, z_max, num_theta)

for i, n in enumerate(orders):
    print(f"Computing flat space: n={n}")
    I_map = compute_intensity_map_flat(n, sigma_flat, k, z_array, h_array)

    I_max = I_map.max()
    if I_max > 0:
        I_plot = I_map / I_max
    else:
        I_plot = I_map

    ax = axes2[i]
    extent = [z_array[0], z_array[-1], h_array[0], h_array[-1]]
    im = ax.imshow(I_plot.T, extent=extent, aspect='auto', origin='lower',
                   cmap='hot', vmin=0, vmax=1)
    ax.set_ylabel('h')
    if i == len(orders) - 1:
        ax.set_xlabel('z')
    ax.text(0.02, 0.95, f'n={n}', transform=ax.transAxes,
            color='white', fontsize=10, va='top')

    key = f"flat_n{n}"
    all_data[key] = I_map

plt.tight_layout()
plt.savefig(os.path.join(plot_dir, 'fig2b_flat_propagation.png'), dpi=150, bbox_inches='tight')
plt.close()

# Save data as CSV (save a representative slice at several theta positions)
theta_indices = np.linspace(0, num_theta - 1, 20, dtype=int)
header_parts = ['h']
rows = []
for idx in theta_indices:
    for gd_idx, gamma_d in enumerate(gamma_d_values):
        for n in orders:
            header_parts.append(f'I_gd{gamma_d}_n{n}_theta_idx{idx}')

# Save transverse profiles at selected theta values
with open(os.path.join(data_dir, 'fig2_propagation.csv'), 'w') as f:
    f.write('# Fig 2: HGB intensity profiles on curved surface at selected propagation angles\n')
    f.write(f'# r0={r0}, r={r}, k={k}\n')
    f.write(f'# gamma_d values: {gamma_d_values}\n')
    f.write(f'# Orders n: {orders}\n')
    f.write(f'# theta indices (out of {num_theta}): {list(theta_indices)}\n')
    f.write(f'# Period = {period:.8f}\n')

    # Header
    cols = ['h']
    for idx in theta_indices:
        theta_val = theta_array[idx]
        for gd in gamma_d_values:
            for n in orders:
                cols.append(f'I_gd{gd}_n{n}_theta{theta_val:.4f}')
    f.write(','.join(cols) + '\n')

    for hi in range(num_h):
        row = [f'{h_array[hi]:.8f}']
        for idx in theta_indices:
            for gd_idx, gd in enumerate(gamma_d_values):
                for n in orders:
                    key = f"curved_gd{gd}_n{n}"
                    row.append(f'{all_data[key][idx, hi]:.8e}')
        f.write(','.join(row) + '\n')

print("Fig 2 completed.")
