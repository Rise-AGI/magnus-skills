"""
Figure 4: Fractional order HGB propagation on curved surface.

Reproduces transmission images for fractional and half-integer order beams,
showing that half-integer orders remain hollow on axis.
"""
import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) if os.path.dirname(os.path.abspath(__file__)) else '.')
from hgb_core import compute_intensity_map

# Parameters
r0 = 1.0
r_over_r0 = 2.0
r = r_over_r0 * r0
k = 100.0
gamma_d = 2.0
sigma = np.sqrt(2 * r / (k * gamma_d))

period = np.pi * r / r0
num_theta = 300
num_h = 200

# Two full periods
theta_array = np.linspace(0.01, 2 * period - 0.01, num_theta)
h_max = 0.3
h_array = np.linspace(-h_max, h_max, num_h)

data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data')
os.makedirs(data_dir, exist_ok=True)
plot_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'plots')
os.makedirs(plot_dir, exist_ok=True)

# Fractional orders: half-integers and other fractions
frac_orders = [0.5, 1.5, 2.5, 0.3, 0.7, 1.3]

fig, axes = plt.subplots(2, 3, figsize=(15, 8))
fig.suptitle('Fig. 4: Fractional Order HGB on Curved Surface ($\\gamma_d=2.0$)', fontsize=14)

all_data = {}

for idx, n in enumerate(frac_orders):
    print(f"Computing fractional order n={n}")
    I_map = compute_intensity_map(n, sigma, k, r0, r, theta_array, h_array)

    I_max = I_map.max()
    if I_max > 0:
        I_plot = I_map / I_max
    else:
        I_plot = I_map

    ax = axes[idx // 3, idx % 3]
    extent = [theta_array[0] / period, theta_array[-1] / period, h_array[0], h_array[-1]]
    im = ax.imshow(I_plot.T, extent=extent, aspect='auto', origin='lower',
                   cmap='hot', vmin=0, vmax=1)
    ax.set_ylabel('h')
    ax.set_xlabel('$\\theta$ / period')
    is_half_int = abs(n - round(n * 2) / 2) < 0.01 and abs(n - round(n)) > 0.01
    suffix = ' (half-int)' if is_half_int else ''
    ax.set_title(f'n={n}{suffix}')

    all_data[f"n{n}"] = I_map

plt.tight_layout()
plt.savefig(os.path.join(plot_dir, 'fig4_fractional.png'), dpi=150, bbox_inches='tight')
plt.close()

# Save CSV: on-axis intensity for fractional orders
with open(os.path.join(data_dir, 'fig4_fractional.csv'), 'w') as f:
    f.write('# Fig 4: On-axis intensity for fractional order HGB\n')
    f.write(f'# r0={r0}, r={r}, k={k}, gamma_d={gamma_d}\n')

    cols = ['theta']
    for n in frac_orders:
        cols.append(f'I_axis_n{n}')
    f.write(','.join(cols) + '\n')

    # h=0 corresponds to the middle index
    h_center_idx = num_h // 2

    for i in range(num_theta):
        row = [f'{theta_array[i]:.8f}']
        for n in frac_orders:
            row.append(f'{all_data[f"n{n}"][i, h_center_idx]:.8e}')
        f.write(','.join(row) + '\n')

print("Fig 4 completed.")
