"""
Fig 33: Singlet and triplet gaps vs system size for the Heisenberg chain.

Reproduces Figure 33 from Sandvik (2010).
Shows how the singlet and triplet gaps scale with 1/N.
Both gaps close as ~1/N (with logarithmic corrections) for the gapless Heisenberg chain.
"""
import sys
import os
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) if os.path.dirname(os.path.abspath(__file__)) else '.')
from heisenberg_ed import compute_gaps_sz_sectors

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    chain_lengths = [8, 10, 12, 14, 16, 18, 20, 24]
    J = 1.0

    gs_energies = []
    singlet_gaps = []
    triplet_gaps = []

    for N in chain_lengths:
        print(f"Computing N={N}...")
        gs_e, sg, tg = compute_gaps_sz_sectors(N, J=J)
        gs_energies.append(gs_e)
        singlet_gaps.append(sg)
        triplet_gaps.append(tg)
        print(f"  E0/N = {gs_e/N:.8f}, Delta_s = {sg:.8f}, Delta_t = {tg:.8f}")

    # Save data
    os.makedirs('../data', exist_ok=True)
    header = 'N,E0_per_site,singlet_gap,triplet_gap'
    data_array = np.column_stack([
        chain_lengths,
        np.array(gs_energies) / np.array(chain_lengths),
        singlet_gaps,
        triplet_gaps
    ])
    np.savetxt('../data/fig33_gaps.csv', data_array, delimiter=',',
               header=header, fmt='%.8f', comments='')

    # Plot
    os.makedirs('../plots', exist_ok=True)
    inv_N = 1.0 / np.array(chain_lengths)

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(inv_N, singlet_gaps, 'o-', label='Singlet gap $\\Delta_s$', markersize=6)
    ax.plot(inv_N, triplet_gaps, 's-', label='Triplet gap $\\Delta_t$', markersize=6)
    ax.set_xlabel('$1/N$', fontsize=14)
    ax.set_ylabel('Gap / $J$', fontsize=14)
    ax.set_title('Singlet and Triplet Gaps of the Heisenberg Chain', fontsize=14)
    ax.legend(fontsize=12)
    ax.set_xlim(0, 0.15)
    ax.set_ylim(bottom=0)
    plt.tight_layout()
    plt.savefig('../plots/fig33_gaps.png', dpi=150)
    print("Saved fig33_gaps.png")


if __name__ == '__main__':
    main()
