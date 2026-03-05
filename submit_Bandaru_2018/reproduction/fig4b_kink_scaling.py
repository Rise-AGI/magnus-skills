"""
Fig 4b: Linear growth rate of the resistive internal kink mode
as a function of normalized resistivity for different RE current fractions.

Reproduces the kink mode scaling from Section 3.3.
Shows how RE current fraction affects the S^(-1/3) scaling recovery.
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(os.path.join(os.getcwd(), __file__))) if '__file__' not in dir() else os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from re_fluid import kink_growth_rate_scaling

def main():
    # Range of inverse Lundquist numbers
    S_inv = np.logspace(-7, -2, 200)

    # Three RE current fractions
    fractions = [0.0, 0.5, 1.0]
    labels = [r'$I_r/I_p = 0$', r'$I_r/I_p = 0.5$', r'$I_r/I_p = 1.0$']
    colors = ['blue', 'red', 'green']
    styles = ['-', '--', ':']

    data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data')
    os.makedirs(data_dir, exist_ok=True)
    plots_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'plots')
    os.makedirs(plots_dir, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8, 6))

    all_data = [S_inv]
    col_names = ['S_inv']

    for f_RE, label, color, style in zip(fractions, labels, colors, styles):
        gamma = kink_growth_rate_scaling(S_inv, f_RE=f_RE)
        ax.loglog(S_inv, gamma, color=color, linestyle=style,
                 linewidth=2, label=label)
        all_data.append(gamma)
        col_names.append(f'gamma_fRE_{f_RE}')

    # Reference scaling lines
    ax.loglog(S_inv, 0.5 * S_inv**(1./3.), 'k--', alpha=0.3, linewidth=1,
             label=r'$S^{-1/3}$ scaling')

    ax.set_xlabel(r'$S^{-1}$ (normalized resistivity)', fontsize=12)
    ax.set_ylabel(r'$\gamma \tau_A$', fontsize=12)
    ax.set_title('Fig 4b: Kink mode growth rate vs resistivity', fontsize=13)
    ax.legend(fontsize=11)
    ax.set_xlim(1e-7, 1e-2)
    ax.grid(True, alpha=0.3, which='both')
    fig.tight_layout()
    fig.savefig(os.path.join(plots_dir, 'fig4b_kink_scaling.png'), dpi=150)
    plt.close(fig)

    header = ','.join(col_names)
    np.savetxt(os.path.join(data_dir, 'fig4b_kink_scaling.csv'),
               np.column_stack(all_data), delimiter=',',
               header=header, comments='', fmt='%.8e')
    print("Saved fig4b_kink_scaling.csv and fig4b_kink_scaling.png")


if __name__ == '__main__':
    main()
