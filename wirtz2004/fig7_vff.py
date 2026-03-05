"""
Figure 7 (e-f): VFF phonon dispersion of graphene.
Reproduces VFF model fits from Fig. 7 panels e and f.

Reference: Wirtz & Rubio, cond-mat/0404637, Table 3 and Fig. 7.
"""

import sys
import os
import numpy as np

sys.path.insert(0, os.path.dirname(__file__) if '__file__' in dir() else '.')
sys.path.insert(0, 'reproduction')

from graphene_phonons import (
    compute_dispersion, A_GGA,
    params_vff_aizawa1990, params_vff_siebentritt1997, params_vff_gga,
)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    a = A_GGA
    n_points = 200

    parameterizations = [
        ('e', 'VFF: Aizawa (1990, solid) & Siebentritt (1997, dashed)',
         [('Aizawa (1990)', params_vff_aizawa1990(), '-', 47.3),
          ('Siebentritt (1997)', params_vff_siebentritt1997(), '--', 55.0)]),
        ('f', 'VFF fit to GGA (this work)',
         [('VFF fit to GGA', params_vff_gga(), '-', 33.6)]),
    ]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    for idx, (label, title, param_list) in enumerate(parameterizations):
        ax = axes[idx]

        for name, params, ls, sigma in param_list:
            q_dist, freqs, ticks, tick_labels = compute_dispersion(
                a, params, 'vff', n_points)

            # Save data
            safe_name = name.replace(' ', '_').replace('(', '').replace(')', '')
            csv_path = os.path.join('data', f'fig7{label}_vff_{safe_name}.csv')
            header = 'q_dist,ZA,TA,LA,ZO,TO,LO'
            data = np.column_stack([q_dist, freqs])
            np.savetxt(csv_path, data, delimiter=',', header=header,
                       comments='', fmt='%.6f')

            for i in range(6):
                ax.plot(q_dist, freqs[:, i], f'k{ls}', linewidth=0.8)

        for t in ticks:
            ax.axvline(t, color='gray', linewidth=0.5, linestyle='--')

        ax.set_xticks(ticks)
        ax.set_xticklabels(tick_labels)
        ax.set_ylim(0, 1700)
        ax.set_ylabel(r'$\omega$ (cm$^{-1}$)')
        ax.set_title(f'{label}) {title}')

    fig.tight_layout()
    fig.savefig(os.path.join('plots', 'fig7_vff.png'), dpi=150)
    print('Saved plots/fig7_vff.png')

    # Verification
    for label, title, param_list in parameterizations:
        for name, params, _, sigma in param_list:
            q_dist, freqs, ticks, _ = compute_dispersion(a, params, 'vff', n_points)
            print(f'\n{name} (sigma={sigma}):')
            print(f'  Gamma: {np.round(freqs[0], 1)}')
            m_idx = np.argmin(np.abs(q_dist - ticks[1]))
            k_idx = np.argmin(np.abs(q_dist - ticks[2]))
            print(f'  M: {np.round(freqs[m_idx], 1)}')
            print(f'  K: {np.round(freqs[k_idx], 1)}')


main()
