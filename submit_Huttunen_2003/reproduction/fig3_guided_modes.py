"""
Figure 3: Slab waveguide guided mode analysis.

Computes guided mode effective indices for the photonic crystal slab
waveguide as a function of boundary material dielectric constant.
This demonstrates the core mechanism: the guided mode properties
change with boundary material, enabling light control.

Parameters: eps_rod=13, eps_air=1, l/P=0.2, h/P=0.5
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import brentq
import os
import sys

script_dir = os.path.dirname(os.path.abspath(sys.argv[0])) if sys.argv[0] else '.'
data_dir = os.path.join(script_dir, '..', 'data')
plots_dir = os.path.join(script_dir, '..', 'plots')
os.makedirs(data_dir, exist_ok=True)
os.makedirs(plots_dir, exist_ok=True)

# Parameters
eps_rod = 13.0
eps_air = 1.0
l_frac = 0.2
h_frac = 0.5

# Average slab dielectric constant
eps_avg = l_frac * eps_rod + (1.0 - l_frac) * eps_air  # = 3.4
n_avg = np.sqrt(eps_avg)

print(f"Average slab dielectric: eps_avg = {eps_avg:.2f}, n_avg = {n_avg:.4f}")


def find_guided_modes(omega_norm, n_slab, n_boundary, h):
    """
    Find guided mode effective indices for a symmetric slab waveguide.

    Solves the TE eigenvalue equation:
    tan(kz * h/2) = gamma / kz  (even modes)
    -cot(kz * h/2) = gamma / kz  (odd modes)

    where kz = k0 * sqrt(n_slab^2 - n_eff^2)
          gamma = k0 * sqrt(n_eff^2 - n_boundary^2)
    """
    k0 = omega_norm * 2.0 * np.pi
    modes = []

    if n_slab <= n_boundary:
        return modes

    # Search for even modes
    n_search = np.linspace(n_boundary + 1e-6, n_slab - 1e-6, 10000)
    for mode_type, sign in [('even', 1), ('odd', -1)]:
        prev_val = None
        for n_eff in n_search:
            kz = k0 * np.sqrt(n_slab ** 2 - n_eff ** 2)
            gamma = k0 * np.sqrt(n_eff ** 2 - n_boundary ** 2)
            if mode_type == 'even':
                val = np.tan(kz * h / 2.0) - gamma / kz
            else:
                val = -1.0 / np.tan(kz * h / 2.0 + 1e-15) - gamma / kz

            if prev_val is not None and prev_val * val < 0 and abs(prev_val) < 100 and abs(val) < 100:
                try:
                    def eq(ne):
                        kz_ = k0 * np.sqrt(n_slab ** 2 - ne ** 2)
                        gamma_ = k0 * np.sqrt(ne ** 2 - n_boundary ** 2)
                        if mode_type == 'even':
                            return np.tan(kz_ * h / 2.0) - gamma_ / kz_
                        else:
                            return -1.0 / np.tan(kz_ * h / 2.0 + 1e-15) - gamma_ / kz_

                    n_low = n_search[np.where(n_search == n_eff)[0][0] - 1]
                    n_eff_sol = brentq(eq, n_low, n_eff)
                    modes.append((n_eff_sol, mode_type))
                except (ValueError, RuntimeError):
                    pass
            prev_val = val

    return modes


# Scan boundary material
eps_b_values = np.linspace(1.0, 12.0, 50)
omega_values = np.linspace(0.1, 0.5, 5)

# For each omega, find guided modes as function of eps_b
print("\nGuided mode analysis:")
results = {}

for omega_norm in omega_values:
    n_eff_vs_epsb = []
    for eps_b in eps_b_values:
        n_b = np.sqrt(eps_b)
        modes = find_guided_modes(omega_norm, n_avg, n_b, h_frac)
        if modes:
            n_eff_vs_epsb.append((eps_b, modes[0][0], modes[0][1]))
        else:
            n_eff_vs_epsb.append((eps_b, np.nan, 'none'))
    results[omega_norm] = n_eff_vs_epsb

# Also compute cutoff boundary eps for each frequency
cutoff_eps = []
for omega_norm in np.linspace(0.1, 0.5, 100):
    # Guided mode cutoff when n_boundary >= n_slab
    # Actually cutoff when V = k0*h*sqrt(n_slab^2 - n_b^2)/2 < pi/2 (first mode)
    k0 = omega_norm * 2.0 * np.pi
    # V = pi/2 -> k0*h*sqrt(eps_avg - eps_b)/2 = pi/2
    # eps_b_cutoff = eps_avg - (pi/(k0*h))^2
    eps_b_cutoff = eps_avg - (np.pi / (k0 * h_frac)) ** 2
    if eps_b_cutoff > 0:
        cutoff_eps.append((omega_norm, eps_b_cutoff))

cutoff_eps = np.array(cutoff_eps) if cutoff_eps else np.zeros((0, 2))

# Save data
for omega_norm, data in results.items():
    arr = np.array([(d[0], d[1] if not np.isnan(d[1]) else -1) for d in data])
    np.savetxt(os.path.join(data_dir, f'fig3_neff_omega{omega_norm:.2f}.csv'),
               arr, delimiter=',',
               header=f'eps_b,n_eff (omega_norm={omega_norm:.2f})',
               comments='# ')

if len(cutoff_eps) > 0:
    np.savetxt(os.path.join(data_dir, 'fig3_cutoff.csv'),
               cutoff_eps, delimiter=',',
               header='omega_normalized,eps_b_cutoff',
               comments='# Guided mode cutoff boundary eps\n# ')

# Plot 1: n_eff vs eps_b for different frequencies
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(omega_values)))
for i, omega_norm in enumerate(omega_values):
    data = results[omega_norm]
    eps_b_arr = [d[0] for d in data]
    n_eff_arr = [d[1] for d in data]
    valid = [not np.isnan(n) for n in n_eff_arr]
    eps_valid = [e for e, v in zip(eps_b_arr, valid) if v]
    n_valid = [n for n, v in zip(n_eff_arr, valid) if v]
    if eps_valid:
        ax1.plot(eps_valid, n_valid, '-o', color=colors[i], markersize=2,
                 label=rf'$\omega P/(2\pi c) = {omega_norm:.2f}$')

ax1.plot(eps_b_values, np.sqrt(eps_b_values), 'k--', alpha=0.5, label=r'$n_b = \sqrt{\varepsilon_b}$ (cutoff)')
ax1.set_xlabel(r'Boundary $\varepsilon_b$')
ax1.set_ylabel(r'Effective mode index $n_{eff}$')
ax1.set_title('Guided mode index vs boundary material')
ax1.legend(fontsize=8)
ax1.grid(True, alpha=0.3)
ax1.set_xlim(1, 12)

# Plot 2: Cutoff condition
if len(cutoff_eps) > 0:
    ax2.fill_between(cutoff_eps[:, 0], 0, cutoff_eps[:, 1], alpha=0.3, color='green', label='Guided mode exists')
    ax2.fill_between(cutoff_eps[:, 0], cutoff_eps[:, 1], eps_avg, alpha=0.3, color='red', label='No guided mode (cutoff)')
    ax2.plot(cutoff_eps[:, 0], cutoff_eps[:, 1], 'k-', linewidth=2, label='Cutoff boundary')
    ax2.axhline(y=eps_avg, color='gray', linestyle=':', label=rf'$\varepsilon_{{avg}} = {eps_avg:.1f}$')
    ax2.axvline(x=0.3, color='green', linestyle='--', alpha=0.5, label=r'$\omega = 0.3$')

ax2.set_xlabel(r'Normalized frequency $\omega P/(2\pi c)$')
ax2.set_ylabel(r'Boundary $\varepsilon_b$')
ax2.set_title('Guided mode existence map')
ax2.legend(fontsize=8)
ax2.grid(True, alpha=0.3)
ax2.set_xlim(0.1, 0.5)
ax2.set_ylim(0, 5)

fig.suptitle('Slab waveguide guided mode analysis\n'
             rf'$\varepsilon_{{avg}} = {eps_avg:.1f}$, $h/P = {h_frac}$',
             fontsize=12)
fig.tight_layout()
fig.savefig(os.path.join(plots_dir, 'fig3_guided_mode_analysis.png'), dpi=150, bbox_inches='tight')
plt.close(fig)

print("Done: fig3_guided_mode_analysis.png")
