"""
Figure 5: Skewness of HGB intensity distribution vs order n.

Shows that skewness oscillates with n, reaching zero at integer and half-integer values.
The skewness in two adjacent periods is complementary.
"""
import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) if os.path.dirname(os.path.abspath(__file__)) else '.')
from hgb_core import compute_skewness

# Parameters
r0 = 1.0
r_over_r0 = 2.0
r = r_over_r0 * r0
k = 100.0
gamma_d = 2.0
sigma = np.sqrt(2 * r / (k * gamma_d))

period = np.pi * r / r0
# Half-period position within two adjacent periods
theta_half_1 = period / 2  # half-period in first period
theta_half_2 = period + period / 2  # half-period in second period

h_max = 0.4
num_h = 300
h_array = np.linspace(-h_max, h_max, num_h)

data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data')
os.makedirs(data_dir, exist_ok=True)
plot_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'plots')
os.makedirs(plot_dir, exist_ok=True)

# Scan n from 0.1 to 5.0
n_values = np.linspace(0.1, 5.0, 100)
skewness_period1 = np.zeros_like(n_values)
skewness_period2 = np.zeros_like(n_values)

for i, n in enumerate(n_values):
    if i % 10 == 0:
        print(f"Computing skewness: n={n:.2f} ({i+1}/{len(n_values)})")
    skewness_period1[i] = compute_skewness(n, sigma, k, r0, r, theta_half_1, h_array)
    skewness_period2[i] = compute_skewness(n, sigma, k, r0, r, theta_half_2, h_array)

# Plot
fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(n_values, skewness_period1, 'r-', linewidth=1.5, label='Period 1 (solid)')
ax.plot(n_values, skewness_period2, 'k--', linewidth=1.5, label='Period 2 (dashed)')
ax.axhline(y=0, color='gray', linestyle='-', alpha=0.3)

# Mark integer and half-integer positions
for ni in range(1, 6):
    ax.axvline(x=ni, color='blue', linestyle=':', alpha=0.3)
for ni in [0.5, 1.5, 2.5, 3.5, 4.5]:
    ax.axvline(x=ni, color='green', linestyle=':', alpha=0.3)

ax.set_xlabel('n', fontsize=12)
ax.set_ylabel('Skewness', fontsize=12)
ax.set_title('Fig. 5: Skewness vs Order n ($\\gamma_d=2.0$)', fontsize=14)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(plot_dir, 'fig5_skewness.png'), dpi=150, bbox_inches='tight')
plt.close()

# Save CSV
with open(os.path.join(data_dir, 'fig5_skewness.csv'), 'w') as f:
    f.write('# Fig 5: Skewness of HGB intensity distribution vs order n\n')
    f.write(f'# r0={r0}, r={r}, k={k}, gamma_d={gamma_d}\n')
    f.write(f'# theta_half_1={theta_half_1:.8f}, theta_half_2={theta_half_2:.8f}\n')
    f.write('n,skewness_period1,skewness_period2\n')

    for i in range(len(n_values)):
        f.write(f'{n_values[i]:.8f},{skewness_period1[i]:.8e},{skewness_period2[i]:.8e}\n')

print("Fig 5 completed.")
