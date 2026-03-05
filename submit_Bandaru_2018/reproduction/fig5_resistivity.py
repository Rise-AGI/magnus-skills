"""
Fig 5: Normalized resistivity profile as a function of normalized poloidal flux.

Reproduces the resistivity profile used for the ITER VDE simulation (Section 4).
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(os.path.join(os.getcwd(), __file__))) if '__file__' not in dir() else os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from re_fluid import resistivity_profile_vde

def main():
    psi_N = np.linspace(0, 2.0, 500)
    eta_ratio = resistivity_profile_vde(psi_N)

    data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data')
    os.makedirs(data_dir, exist_ok=True)
    plots_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'plots')
    os.makedirs(plots_dir, exist_ok=True)

    header = "psi_N,eta_over_eta_axis"
    np.savetxt(os.path.join(data_dir, 'fig5_resistivity.csv'),
               np.column_stack([psi_N, eta_ratio]),
               delimiter=',', header=header, comments='', fmt='%.8e')

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.semilogy(psi_N, eta_ratio, 'b-', linewidth=2)
    ax.axvline(x=1.0, color='gray', linestyle='--', alpha=0.5, label='LCFS')
    ax.set_xlabel(r'$\psi_N$', fontsize=12)
    ax.set_ylabel(r'$\eta / \eta_{\mathrm{axis}}$', fontsize=12)
    ax.set_title('Fig 5: Normalized resistivity profile', fontsize=13)
    ax.set_xlim(0, 2.0)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)

    # Annotate regions
    ax.annotate('Core', xy=(0.4, 1.2), fontsize=11, ha='center')
    ax.annotate('Halo', xy=(1.25, 4), fontsize=11, ha='center')
    ax.annotate('Vacuum', xy=(1.75, 20), fontsize=11, ha='center')

    fig.tight_layout()
    fig.savefig(os.path.join(plots_dir, 'fig5_resistivity.png'), dpi=150)
    plt.close(fig)
    print("Saved fig5_resistivity.csv and fig5_resistivity.png")


if __name__ == '__main__':
    main()
