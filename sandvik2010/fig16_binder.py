"""
Fig 16: Binder cumulant vs temperature for 2D Ising model.

Reproduces Figure 16 from Sandvik (2010).
U_2 = 1 - <m^4> / (3 <m^2>^2)

The Binder cumulant is dimensionless and curves for different L cross at T_c.
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
        np.linspace(1.8, 2.15, 8),
        np.linspace(2.16, 2.38, 12),
        np.linspace(2.4, 2.8, 6),
    ])
    n_therm = 300
    n_meas = 3000

    all_data = {}
    for L in lattice_sizes:
        print(f"Running L={L}...")
        binder_vals = []
        for T in T_values:
            result = run_ising_simulation(L, T, n_therm=n_therm, n_meas=n_meas,
                                          seed=42 + L * 1000 + int(T * 100))
            binder_vals.append(result['binder'])
            print(f"  T={T:.3f}, U2={result['binder']:.6f}")
        all_data[L] = np.array(binder_vals)

    # Save data
    os.makedirs('../data', exist_ok=True)
    header = 'T,' + ','.join([f'binder_L{L}' for L in lattice_sizes])
    data_array = np.column_stack([T_values] + [all_data[L] for L in lattice_sizes])
    np.savetxt('../data/fig16_binder.csv', data_array, delimiter=',',
               header=header, fmt='%.8f', comments='')

    # Plot
    os.makedirs('../plots', exist_ok=True)
    fig, ax = plt.subplots(figsize=(8, 6))
    markers = ['o', 's', '^', 'D']
    for i, L in enumerate(lattice_sizes):
        ax.plot(T_values, all_data[L], markers[i] + '-', label=f'L={L}',
                markersize=4, linewidth=1)
    ax.axvline(T_c, color='gray', linestyle='--', alpha=0.7, label=f'$T_c$')
    ax.axhline(0, color='lightgray', linestyle='-', alpha=0.5)
    ax.set_xlabel('$T/J$', fontsize=14)
    ax.set_ylabel('$U_2$', fontsize=14)
    ax.set_title('Binder Cumulant vs Temperature (2D Ising)', fontsize=14)
    ax.legend(fontsize=12)
    ax.set_xlim(1.8, 2.8)
    plt.tight_layout()
    plt.savefig('../plots/fig16_binder.png', dpi=150)
    print("Saved fig16_binder.png")


if __name__ == '__main__':
    main()
