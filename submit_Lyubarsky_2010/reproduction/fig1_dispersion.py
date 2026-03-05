"""
Figure 1: Dispersion relation and growth rate of the K-S instability.

Plots the exact growth rate (Eq. 26) vs. wavenumber k*Delta,
compared with the short wavelength (Eq. 28) and long wavelength (Eq. 29) limits.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(os.path.realpath(sys.argv[0]))))
from ks_instability import (
    dispersion_relation_omega4,
    growth_rate_exact,
    growth_rate_short_wavelength,
    growth_rate_long_wavelength,
)

# Normalized units: set g = 1, Delta = 1
g = 1.0
Delta = 1.0

# Wavenumber range
kDelta = np.linspace(0.01, 10.0, 500)
k = kDelta / Delta

# Compute growth rates
eta_exact = growth_rate_exact(k, g, Delta)
eta_short = growth_rate_short_wavelength(k, g)
eta_long = growth_rate_long_wavelength(k, g, Delta)

# Normalize: eta / sqrt(g/Delta) to make dimensionless
norm = np.sqrt(g / Delta)
eta_exact_norm = eta_exact / norm
eta_short_norm = eta_short / norm
eta_long_norm = eta_long / norm

# Save data
data_dir = os.path.join(os.path.dirname(os.path.abspath(os.path.realpath(sys.argv[0]))), '..', 'data')
os.makedirs(data_dir, exist_ok=True)

header = "kDelta,eta_exact,eta_short_wavelength,eta_long_wavelength"
data = np.column_stack([kDelta, eta_exact_norm, eta_short_norm, eta_long_norm])
np.savetxt(os.path.join(data_dir, 'fig1_dispersion.csv'), data,
           delimiter=',', header=header, comments='', fmt='%.8e')

# Plot
fig, ax = plt.subplots(1, 1, figsize=(8, 6))
ax.plot(kDelta, eta_exact_norm, 'k-', linewidth=2, label='Exact (Eq. 26)')
ax.plot(kDelta, eta_short_norm, 'r--', linewidth=1.5,
        label=r'Short $\lambda$: $\eta = \sqrt{kg/3}$ (Eq. 28)')
ax.plot(kDelta, eta_long_norm, 'b-.', linewidth=1.5,
        label=r'Long $\lambda$: $\eta = (g/2)^{1/2}\Delta^{1/4}k^{3/4}$ (Eq. 29)')

ax.set_xlabel(r'$k\Delta$', fontsize=14)
ax.set_ylabel(r'$\eta / \sqrt{g/\Delta}$', fontsize=14)
ax.set_title('Kruskal-Schwarzschild Instability Growth Rate', fontsize=14)
ax.legend(fontsize=11)
ax.set_xlim(0, 10)
ax.set_ylim(0, None)
ax.grid(True, alpha=0.3)

plots_dir = os.path.join(os.path.dirname(os.path.abspath(os.path.realpath(sys.argv[0]))), '..', 'plots')
os.makedirs(plots_dir, exist_ok=True)
fig.tight_layout()
fig.savefig(os.path.join(plots_dir, 'fig1_dispersion.png'), dpi=150)
plt.close(fig)
print("Figure 1 saved.")
