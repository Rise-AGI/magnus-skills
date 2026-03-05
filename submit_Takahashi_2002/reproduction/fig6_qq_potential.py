"""
Figure 6: Q-Qbar potential at beta=5.7, 5.8, 6.0.

Plots the quark-antiquark potential from the Y-ansatz fit parameters,
showing the Coulomb + linear + constant form.
"""
import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from three_quark_potential import FIT_PARAMS, vqq_potential, LATTICE_PARAMS


def main():
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=False)
    csv_lines = ['# Q-Qbar potential fit at various beta\n']
    csv_lines.append('beta,r_lattice,V_QQ\n')

    for ax_idx, beta in enumerate([5.7, 5.8, 6.0]):
        params = FIT_PARAMS[beta]["QQ"]
        A, sigma, C = params["A"], params["sigma"], params["C"]
        a_fm = LATTICE_PARAMS[beta]["a_fm"]

        r = np.linspace(0.5, 8, 200)
        v = vqq_potential(r, A, sigma, C)

        axes[ax_idx].plot(r, v, 'b-', linewidth=2, label=r'$V_{Q\bar{Q}}$ fit')
        axes[ax_idx].set_xlabel(r'$r$ (lattice units)', fontsize=12)
        axes[ax_idx].set_ylabel(r'$V_{Q\bar{Q}}$ (lattice units)', fontsize=12)
        axes[ax_idx].set_title(r'$\beta={:.1f}$, $a={:.2f}$ fm'.format(beta, a_fm), fontsize=12)
        axes[ax_idx].grid(True, alpha=0.3)
        axes[ax_idx].legend(fontsize=10)

        for ri in r:
            csv_lines.append('{:.1f},{:.8f},{:.8f}\n'.format(beta, ri, vqq_potential(ri, A, sigma, C)))

    os.makedirs('../data', exist_ok=True)
    with open('../data/fig6_qq_potential.csv', 'w') as f:
        f.writelines(csv_lines)

    os.makedirs('../plots', exist_ok=True)
    plt.tight_layout()
    plt.savefig('../plots/fig6_qq_potential.png', dpi=150)
    plt.close()
    print('Fig 6 saved.')


if __name__ == '__main__':
    main()
