"""
Figure 2: Lorentz factor evolution gamma(r) for different astrophysical scenarios.

Computes Eq. 43: gamma = (9 * zeta^2 * gamma_max * r / (4l))^(1/3)
for pulsars, AGN jets, and GRB outflows.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(os.path.realpath(sys.argv[0]))))
from ks_instability import lorentz_factor, dissipation_scale

# Astrophysical parameters from Section 4
# zeta is a free parameter; use zeta = 0.1 as a representative value
scenarios = {
    'AGN': {
        'gamma_max': 10.0,
        'l_label': 'r_g',
        'l_cm': 1.5e14,       # r_g for 10^9 M_sun BH
        'zeta': 0.1,
        'description': 'AGN jet',
    },
    'GRB': {
        'gamma_max': 1000.0,
        'l_label': 'r_g',
        'l_cm': 4.4e5,         # r_g for ~3 M_sun
        'zeta': 0.1,
        'description': 'GRB outflow',
    },
    'Pulsar': {
        'gamma_max': 1e4,
        'l_label': 'r_L',
        'l_cm': 1.6e6 * np.pi,  # pi * r_L for Crab-like (P~33ms)
        'zeta': 0.1,
        'description': 'Pulsar wind',
    },
}

data_dir = os.path.join(os.path.dirname(os.path.abspath(os.path.realpath(sys.argv[0]))), '..', 'data')
os.makedirs(data_dir, exist_ok=True)
plots_dir = os.path.join(os.path.dirname(os.path.abspath(os.path.realpath(sys.argv[0]))), '..', 'plots')
os.makedirs(plots_dir, exist_ok=True)

fig, ax = plt.subplots(1, 1, figsize=(8, 6))
colors = {'AGN': 'blue', 'GRB': 'red', 'Pulsar': 'green'}
styles = {'AGN': '-', 'GRB': '--', 'Pulsar': '-.'}

all_data = {}
for name, params in scenarios.items():
    gmax = params['gamma_max']
    l = params['l_cm']
    zeta = params['zeta']
    r_diss = dissipation_scale(gmax, zeta, l)

    # Plot from some initial radius to dissipation scale
    r = np.logspace(np.log10(l), np.log10(2 * r_diss), 500)
    gamma = lorentz_factor(r, zeta, gmax, l)
    # Cap at gamma_max
    gamma = np.minimum(gamma, gmax)

    r_over_l = r / l
    all_data[name] = (r_over_l, gamma)

    ax.loglog(r_over_l, gamma, color=colors[name], linestyle=styles[name],
              linewidth=2, label=f'{name} ($\\gamma_{{max}}={gmax:.0f}$)')

    # Mark dissipation scale
    ax.axvline(r_diss / l, color=colors[name], alpha=0.3, linestyle=':')

# Save data
max_len = max(len(v[0]) for v in all_data.values())
header = "r_over_l_AGN,gamma_AGN,r_over_l_GRB,gamma_GRB,r_over_l_Pulsar,gamma_Pulsar"
out = np.full((max_len, 6), np.nan)
for i, name in enumerate(['AGN', 'GRB', 'Pulsar']):
    r_l, gam = all_data[name]
    out[:len(r_l), 2*i] = r_l
    out[:len(r_l), 2*i+1] = gam
np.savetxt(os.path.join(data_dir, 'fig2_lorentz_evolution.csv'), out,
           delimiter=',', header=header, comments='', fmt='%.8e')

ax.set_xlabel(r'$r / l$', fontsize=14)
ax.set_ylabel(r'$\gamma$', fontsize=14)
ax.set_title(r'Lorentz Factor Evolution (Eq. 43, $\zeta=0.1$)', fontsize=14)
ax.legend(fontsize=11)
ax.grid(True, alpha=0.3, which='both')
fig.tight_layout()
fig.savefig(os.path.join(plots_dir, 'fig2_lorentz_evolution.png'), dpi=150)
plt.close(fig)
print("Figure 2 saved.")
