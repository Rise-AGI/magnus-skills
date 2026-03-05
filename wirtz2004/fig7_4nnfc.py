"""
Figure 7 (a-d): 4NNFC phonon dispersion of graphene.
Reproduces the 4 panels of Fig. 7 showing different 4NNFC parameterizations
compared to the GGA ab-initio calculation.

Reference: Wirtz & Rubio, cond-mat/0404637, Table 3 and Fig. 7.
"""

import sys
import os
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(".")))
sys.path.insert(0, os.path.dirname(__file__) if '__file__' in dir() else '.')
sys.path.insert(0, 'reproduction')

from graphene_phonons import (
    compute_dispersion, A_GGA,
    params_4nnfc_jishi1993, params_4nnfc_gruneis2002,
    params_4nnfc_gga_diag, params_4nnfc_gga_offdiag,
)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    a = A_GGA
    n_points = 200

    parameterizations = [
        ('a', 'Jishi et al. (1993)', params_4nnfc_jishi1993(), 51.5),
        ('b', 'Gruneis et al. (2002)', params_4nnfc_gruneis2002(), 69.5),
        ('c', '4NNFC fit to GGA (this work)', params_4nnfc_gga_diag(), 15.4),
        ('d', '4NNFC + off-diag fit to GGA', params_4nnfc_gga_offdiag(), 13.5),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.flatten()

    for idx, (label, title, params, sigma) in enumerate(parameterizations):
        ax = axes[idx]
        q_dist, freqs, ticks, tick_labels = compute_dispersion(a, params, '4nnfc', n_points)

        # Save data to CSV
        csv_path = os.path.join('data', f'fig7{label}_4nnfc.csv')
        header = 'q_dist,ZA,TA,LA,ZO,TO,LO'
        data = np.column_stack([q_dist, freqs])
        np.savetxt(csv_path, data, delimiter=',', header=header, comments='',
                   fmt='%.6f')

        # Plot all 6 branches
        branch_labels = ['ZA', 'TA', 'LA', 'ZO', 'TO', 'LO']
        for i in range(6):
            ax.plot(q_dist, freqs[:, i], 'k-', linewidth=0.8)

        for t in ticks:
            ax.axvline(t, color='gray', linewidth=0.5, linestyle='--')

        ax.set_xticks(ticks)
        ax.set_xticklabels(tick_labels)
        ax.set_ylim(0, 1700)
        ax.set_ylabel(r'$\omega$ (cm$^{-1}$)')
        ax.set_title(f'{label}) {title} ($\\sigma$ = {sigma} cm$^{{-1}}$)')

    fig.tight_layout()
    fig.savefig(os.path.join('plots', 'fig7_4nnfc.png'), dpi=150)
    print('Saved plots/fig7_4nnfc.png')

    # Print high-symmetry point frequencies for verification
    for label, title, params, sigma in parameterizations:
        q_dist, freqs, ticks, _ = compute_dispersion(a, params, '4nnfc', n_points)
        print(f'\n{label}) {title}:')
        print(f'  Gamma: {np.round(freqs[0], 1)}')
        # Find M and K points
        m_idx = np.argmin(np.abs(q_dist - ticks[1]))
        k_idx = np.argmin(np.abs(q_dist - ticks[2]))
        print(f'  M: {np.round(freqs[m_idx], 1)}')
        print(f'  K: {np.round(freqs[k_idx], 1)}')


main()
