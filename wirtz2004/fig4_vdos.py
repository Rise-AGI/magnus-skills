"""
Figure 4: Vibrational density of states (vDOS) of graphene.
Computes vDOS from the 4NNFC model and the VFF model of Aizawa et al.

Reference: Wirtz & Rubio, cond-mat/0404637, Fig. 4.
"""

import sys
import os
import numpy as np

sys.path.insert(0, os.path.dirname(__file__) if '__file__' in dir() else '.')
sys.path.insert(0, 'reproduction')

from graphene_phonons import (
    compute_vdos, A_GGA,
    params_4nnfc_gga_diag, params_4nnfc_jishi1993,
    params_vff_aizawa1990,
)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    a = A_GGA
    n_grid = 80
    n_bins = 600
    freq_max = 1700

    models = [
        ('4NNFC fit to GGA', params_4nnfc_gga_diag(), '4nnfc', 'k-'),
        ('4NNFC Jishi (1993)', params_4nnfc_jishi1993(), '4nnfc', 'b--'),
        ('VFF Aizawa (1990)', params_vff_aizawa1990(), 'vff', 'r-.'),
    ]

    fig, ax = plt.subplots(figsize=(8, 5))

    all_data = []
    headers = ['freq_cm-1']

    for name, params, model, ls in models:
        print(f'Computing vDOS for {name}...')
        freq_bins, dos = compute_vdos(a, params, model, n_grid, n_bins, freq_max)
        ax.plot(freq_bins, dos, ls, linewidth=1.0, label=name)
        all_data.append(dos)
        headers.append(name.replace(' ', '_').replace('(', '').replace(')', ''))

    ax.set_xlabel(r'Phonon energy (cm$^{-1}$)')
    ax.set_ylabel('DOS (arb. units)')
    ax.set_xlim(0, freq_max)
    ax.legend()
    ax.set_title('Vibrational density of states of graphene')

    fig.tight_layout()
    fig.savefig(os.path.join('plots', 'fig4_vdos.png'), dpi=150)
    print('Saved plots/fig4_vdos.png')

    # Save data
    data = np.column_stack([freq_bins] + all_data)
    header_str = ','.join(headers)
    np.savetxt(os.path.join('data', 'fig4_vdos.csv'), data,
               delimiter=',', header=header_str, comments='', fmt='%.8f')
    print('Saved data/fig4_vdos.csv')


main()
