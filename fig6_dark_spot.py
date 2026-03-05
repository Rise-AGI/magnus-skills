"""
Figure 6: Dark spot size analysis.

Computes the longitudinal dark spot size for type A (Eq. 19) and type B (Eq. 20)
as a function of r/r0 and n. Shows how curved surface modifies the dark region
compared to flat space.
"""
import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import brentq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) if os.path.dirname(os.path.abspath(__file__)) else '.')
from hgb_core import dark_spot_size_type_B

# Parameters
r0 = 1.0
k = 100.0

data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data')
os.makedirs(data_dir, exist_ok=True)
plot_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'plots')
os.makedirs(plot_dir, exist_ok=True)


def dark_spot_type_A_eq19(n, gamma_d, r0, r):
    """
    Solve Eq. 19 numerically for w_theta (type A dark spot size).
    (gamma_d^2 + cot^2(r0*w/r))^(-n-0.5) / sin(r0*w/r) = gamma_d^(-2n-1) / 2
    """
    target = gamma_d**(-2*n - 1) / 2.0

    def equation(w):
        sin_val = np.sin(r0 * w / r)
        if abs(sin_val) < 1e-15:
            return 1e10
        cot_val = np.cos(r0 * w / r) / sin_val
        lhs = (gamma_d**2 + cot_val**2)**(-n - 0.5) / abs(sin_val)
        return lhs - target

    # Search in (0, pi*r/(2*r0))
    period_half = np.pi * r / (2 * r0)
    try:
        w_sol = brentq(equation, 0.01, period_half - 0.01, maxiter=200)
        return w_sol
    except (ValueError, RuntimeError):
        return np.nan


# Compute dark spot sizes for different n and gamma_d
orders = [1, 2, 3, 4, 5, 10]
gamma_d_range = np.linspace(0.5, 8.0, 100)

fig, axes = plt.subplots(1, 2, figsize=(14, 6))
fig.suptitle('Fig. 6: Dark Spot Size Analysis', fontsize=14)

# Panel (a): w_theta vs gamma_d for different n
ax = axes[0]
csv_data_a = {'gamma_d': gamma_d_range}

for n in orders:
    boundary = np.sqrt(2 * n + 1)
    w_theta = np.zeros_like(gamma_d_range)

    for i, gd in enumerate(gamma_d_range):
        if gd <= boundary:
            # Type A
            r_val = 2.0 * r0  # r/r0 = 2
            w_theta[i] = dark_spot_type_A_eq19(n, gd, r0, r_val)
        else:
            # Type B (Eq. 20)
            r_val = 2.0 * r0
            w_theta[i] = dark_spot_size_type_B(n, gd, r0, r_val)

    ax.plot(gamma_d_range, w_theta, linewidth=1.5, label=f'n={n}')
    csv_data_a[f'w_theta_n{n}'] = w_theta

ax.set_xlabel('$\\gamma_d$', fontsize=12)
ax.set_ylabel('$w_\\theta$ (dark spot size)', fontsize=12)
ax.set_title('(a) Dark spot size vs $\\gamma_d$ ($r/r_0=2$)')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

# Panel (b): w_theta vs n for different r/r0 values (type B, gamma_d=5)
ax2 = axes[1]
gamma_d_fixed = 5.0
r_over_r0_values = [1.5, 2.0, 3.0, 5.0, 10.0]
n_range = np.arange(1, 11)

csv_data_b = {'n': n_range}

for r_ratio in r_over_r0_values:
    r_val = r_ratio * r0
    w_theta_n = np.zeros_like(n_range, dtype=float)

    for i, n in enumerate(n_range):
        boundary = np.sqrt(2 * n + 1)
        if gamma_d_fixed > boundary:
            w_theta_n[i] = dark_spot_size_type_B(n, gamma_d_fixed, r0, r_val)
        else:
            w_theta_n[i] = dark_spot_type_A_eq19(n, gamma_d_fixed, r0, r_val)

    ax2.plot(n_range, w_theta_n, 'o-', linewidth=1.5, markersize=4,
             label=f'$r/r_0={r_ratio}$')
    csv_data_b[f'w_theta_r{r_ratio}'] = w_theta_n

# Flat space limit: r0 * w_theta -> sqrt(2n) * k * sigma^2 / 2
sigma_flat = np.sqrt(2 * r0 * 2.0 / (k * gamma_d_fixed))  # r/r0=2
w_flat = np.sqrt(2.0 * n_range) * k * sigma_flat**2 / (2 * r0)
ax2.plot(n_range, w_flat, 'k--', linewidth=1.5, label='Flat space limit')
csv_data_b['w_theta_flat'] = w_flat

ax2.set_xlabel('n (order)', fontsize=12)
ax2.set_ylabel('$w_\\theta$ (dark spot size)', fontsize=12)
ax2.set_title(f'(b) Dark spot size vs n ($\\gamma_d={gamma_d_fixed}$)')
ax2.legend(fontsize=9)
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(plot_dir, 'fig6_dark_spot.png'), dpi=150, bbox_inches='tight')
plt.close()

# Save CSV
with open(os.path.join(data_dir, 'fig6_dark_spot_vs_gamma.csv'), 'w') as f:
    f.write('# Fig 6a: Dark spot size vs gamma_d for different n\n')
    f.write(f'# r0={r0}, r/r0=2.0, k={k}\n')
    cols = ['gamma_d'] + [f'w_theta_n{n}' for n in orders]
    f.write(','.join(cols) + '\n')
    for i in range(len(gamma_d_range)):
        row = [f'{gamma_d_range[i]:.8f}']
        for n in orders:
            row.append(f'{csv_data_a[f"w_theta_n{n}"][i]:.8e}')
        f.write(','.join(row) + '\n')

with open(os.path.join(data_dir, 'fig6_dark_spot_vs_n.csv'), 'w') as f:
    f.write('# Fig 6b: Dark spot size vs n for different r/r0\n')
    f.write(f'# r0={r0}, gamma_d={gamma_d_fixed}, k={k}\n')
    cols = ['n'] + [f'w_theta_r{rv}' for rv in r_over_r0_values] + ['w_theta_flat']
    f.write(','.join(cols) + '\n')
    for i in range(len(n_range)):
        row = [f'{n_range[i]}']
        for rv in r_over_r0_values:
            row.append(f'{csv_data_b[f"w_theta_r{rv}"][i]:.8e}')
        row.append(f'{csv_data_b["w_theta_flat"][i]:.8e}')
        f.write(','.join(row) + '\n')

print("Fig 6 completed.")
