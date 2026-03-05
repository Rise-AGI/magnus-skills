"""
Topological susceptibility and delta from index distribution (Table 1).
Also computes eta' mass and verifies consistency of delta.

Reproduces the analysis in Section 2 of Chiu & Hsieh (2003).
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(".")))

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from qxpt_analysis import (
    INDEX_DISTRIBUTION, N_CONFIGS, N_SITES,
    compute_topological_susceptibility, delta_from_chi_t, eta_prime_mass,
    FPI_F0, A_INV_GEV, DELTA, N_F
)

# Compute topological susceptibility from Table 1
a4_chi_t, a4_chi_t_err, mean_Q2 = compute_topological_susceptibility()

# Compute delta from chi_t
delta_topo = delta_from_chi_t(a4_chi_t, FPI_F0)

# Compute eta' mass
meta_GeV = eta_prime_mass(a4_chi_t, FPI_F0, A_INV_GEV)

print("=" * 60)
print("Topological Susceptibility Analysis")
print("=" * 60)
print(f"<(n+ - n-)^2> = {mean_Q2:.2f}")
print(f"a^4 * chi_t = {a4_chi_t:.2e} +/- {a4_chi_t_err:.2e}")
print(f"  (published: {6.03e-5:.2e} +/- {0.75e-5:.2e})")
print(f"\nchi_t = ({(a4_chi_t * A_INV_GEV**4 * 1e12)**(0.25):.0f} MeV)^4")
print(f"  (published: (175 +/- 6 MeV)^4)")
print(f"\ndelta from chi_t = {delta_topo:.3f}")
print(f"  (published: 0.16 +/- 0.02)")
print(f"delta from pion mass fit = {DELTA}")
print(f"  (published: 0.164 +/- 0.013)")
print(f"\nm_eta' = {meta_GeV*1000:.0f} MeV")
print(f"  (published: 813 +/- 51 MeV)")

# Plot: Histogram of topological charges
fig, ax = plt.subplots(figsize=(8, 6))
charges = sorted(INDEX_DISTRIBUTION.keys())
counts = [INDEX_DISTRIBUTION[q] for q in charges]
ax.bar(charges, counts, color='steelblue', edgecolor='navy', alpha=0.8)
ax.set_xlabel('Topological Charge $Q = n_+ - n_-$', fontsize=14)
ax.set_ylabel('Number of Configurations', fontsize=14)
ax.set_title('Distribution of Topological Charges\n'
             f'($\\beta=6.0$, $16^3 \\times 32$, 100 configs)', fontsize=14)
ax.set_xticks(charges)
ax.grid(True, alpha=0.3, axis='y')

# Add annotation
ax.text(0.98, 0.95,
        f'$\\langle Q^2 \\rangle$ = {mean_Q2:.2f}\n'
        f'$\\chi_t$ = (175 MeV)$^4$\n'
        f'$\\delta$ = {delta_topo:.3f}',
        transform=ax.transAxes, fontsize=11,
        verticalalignment='top', horizontalalignment='right',
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

plt.tight_layout()
plt.savefig('../plots/fig_topology.png', dpi=150)
print("\nSaved plots/fig_topology.png")

# Save data
data = np.column_stack([charges, counts])
np.savetxt('../data/topology_distribution.csv', data, delimiter=',',
           header='topological_charge,n_configurations', fmt='%d', comments='')
print("Saved data/topology_distribution.csv")
