"""
Figure 3: A-type and B-type axial intensity distributions.

Reproduces the on-axis intensity along propagation for different orders n,
showing the two types of intensity profiles distinguished by gamma_d vs sqrt(2n+1).
"""
import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) if os.path.dirname(os.path.abspath(__file__)) else '.')
from hgb_core import axial_intensity_integer_n, divergence_coefficient

# Parameters
r0 = 1.0
r_over_r0 = 2.0
r = r_over_r0 * r0
k = 100.0

period = np.pi * r / r0
num_theta = 500
theta_array = np.linspace(0.02, period - 0.02, num_theta)

data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data')
os.makedirs(data_dir, exist_ok=True)
plot_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'plots')
os.makedirs(plot_dir, exist_ok=True)

# For n=1: boundary is gamma_d = sqrt(3) ~ 1.73
# A-type: gamma_d < sqrt(2n+1), B-type: gamma_d > sqrt(2n+1)

# Demonstrate with n=1: A-type (gamma_d=1.0) and B-type (gamma_d=3.0)
# Also show n=4: A-type (gamma_d=1.0) and B-type (gamma_d=5.0)
cases = [
    {'n': 1, 'gamma_d': 1.0, 'type': 'A', 'label': 'n=1, $\\gamma_d$=1.0 (A-type)'},
    {'n': 1, 'gamma_d': 3.0, 'type': 'B', 'label': 'n=1, $\\gamma_d$=3.0 (B-type)'},
    {'n': 4, 'gamma_d': 1.0, 'type': 'A', 'label': 'n=4, $\\gamma_d$=1.0 (A-type)'},
    {'n': 4, 'gamma_d': 5.0, 'type': 'B', 'label': 'n=4, $\\gamma_d$=5.0 (B-type)'},
]

fig, axes = plt.subplots(2, 2, figsize=(12, 8))
fig.suptitle('Fig. 3: A-type and B-type Axial Intensity Distributions', fontsize=14)

csv_data = {'theta': theta_array}

for idx, case in enumerate(cases):
    n = case['n']
    gamma_d = case['gamma_d']
    label = case['label']

    I_axis = axial_intensity_integer_n(n, gamma_d, r0, r, theta_array)

    # Normalize
    I_max = I_axis.max()
    if I_max > 0:
        I_norm = I_axis / I_max
    else:
        I_norm = I_axis

    ax = axes[idx // 2, idx % 2]
    ax.plot(theta_array / period, I_norm, 'b-', linewidth=1.5)
    ax.set_xlabel('$\\theta$ / period')
    ax.set_ylabel('Normalized Intensity')
    ax.set_title(label)
    ax.set_xlim(0, 1)
    ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.3)

    # Mark half-period
    ax.axvline(x=0.5, color='r', linestyle='--', alpha=0.5, label='$\\theta = \\pi r / 2r_0$')
    ax.legend(fontsize=8)

    csv_key = f"I_n{n}_gd{gamma_d}"
    csv_data[csv_key] = I_norm

plt.tight_layout()
plt.savefig(os.path.join(plot_dir, 'fig3_axial_intensity.png'), dpi=150, bbox_inches='tight')
plt.close()

# Save CSV
with open(os.path.join(data_dir, 'fig3_axial_intensity.csv'), 'w') as f:
    f.write('# Fig 3: On-axis intensity for A-type and B-type distributions\n')
    f.write(f'# r0={r0}, r={r}, k={k}\n')
    f.write(f'# Period = {period:.8f}\n')
    f.write('# Boundary: gamma_d = sqrt(2n+1) separates A-type from B-type\n')

    cols = ['theta']
    for case in cases:
        cols.append(f"I_n{case['n']}_gd{case['gamma_d']}")
    f.write(','.join(cols) + '\n')

    for i in range(num_theta):
        row = [f'{theta_array[i]:.8f}']
        for case in cases:
            key = f"I_n{case['n']}_gd{case['gamma_d']}"
            row.append(f'{csv_data[key][i]:.8e}')
        f.write(','.join(row) + '\n')

print("Fig 3 completed.")
