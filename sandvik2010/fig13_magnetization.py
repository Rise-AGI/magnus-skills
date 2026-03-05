"""
Fig 13: Squared magnetization vs temperature for 2D Ising model.

Reproduces Figure 13(a) from Sandvik (2010).
Shows <m^2> as a function of temperature for several lattice sizes.

The exact critical temperature is T_c/J = 2/ln(1+sqrt(2)) ~ 2.26919.
"""
import sys
import os
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) if os.path.dirname(os.path.abspath(__file__)) else '.')
from ising_mc import run_ising_simulation

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    T_c = 2.0 / np.log(1.0 + np.sqrt(2.0))
    lattice_sizes = [8, 16, 32, 64]
    T_values = np.concatenate([
        np.linspace(1.5, 2.0, 6),
        np.linspace(2.05, 2.5, 10),
        np.linspace(2.6, 3.5, 6),
    ])
    n_therm = 300
    n_meas = 3000

    # Output data
    all_data = {}
    for L in lattice_sizes:
        print(f"Running L={L}...")
        m2_vals = []
        for T in T_values:
            result = run_ising_simulation(L, T, n_therm=n_therm, n_meas=n_meas,
                                          seed=42 + L * 1000 + int(T * 100))
            m2_vals.append(result['m2'])
            print(f"  T={T:.3f}, <m^2>={result['m2']:.6f}")
        all_data[L] = np.array(m2_vals)

    # Save data
    os.makedirs('../data', exist_ok=True)
    header = 'T,' + ','.join([f'm2_L{L}' for L in lattice_sizes])
    data_array = np.column_stack([T_values] + [all_data[L] for L in lattice_sizes])
    np.savetxt('../data/fig13_magnetization.csv', data_array, delimiter=',',
               header=header, fmt='%.8f', comments='')

    # Plot
    os.makedirs('../plots', exist_ok=True)
    fig, ax = plt.subplots(figsize=(8, 6))
    markers = ['o', 's', '^', 'D']
    for i, L in enumerate(lattice_sizes):
        ax.plot(T_values, all_data[L], markers[i] + '-', label=f'L={L}',
                markersize=4, linewidth=1)
    ax.axvline(T_c, color='gray', linestyle='--', alpha=0.7, label=f'$T_c$={T_c:.4f}')
    ax.set_xlabel('$T/J$', fontsize=14)
    ax.set_ylabel('$\\langle m^2 \\rangle$', fontsize=14)
    ax.set_title('Squared Magnetization vs Temperature (2D Ising)', fontsize=14)
    ax.legend(fontsize=12)
    ax.set_xlim(1.5, 3.5)
    plt.tight_layout()
    plt.savefig('../plots/fig13_magnetization.png', dpi=150)
    print("Saved fig13_magnetization.png")


if __name__ == '__main__':
    main()
