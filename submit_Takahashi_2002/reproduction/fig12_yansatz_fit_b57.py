"""
Figure 12: Comparison of lattice QCD data V_3Q^latt and Y-ansatz fit at beta=5.7.

Plots V_3Q as a function of i for each fixed (j,k), showing both discrete
lattice data points and continuous fit curves.
"""
import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from three_quark_potential import (DATA_BETA_57, FIT_PARAMS, v3q_y_ansatz,
                                   compute_lmin, quark_distances_type1)


def main():
    data = DATA_BETA_57
    params = FIT_PARAMS[5.7]["3QY"]
    A, sigma, C = params["A"], params["sigma"], params["C"]

    # Group by (j, k)
    groups = {}
    for entry in data:
        cfg = entry[0]
        i, j, k = cfg
        key = (j, k)
        if key not in groups:
            groups[key] = []
        groups[key].append((i, entry[1], entry[2]))

    fig, ax = plt.subplots(1, 1, figsize=(10, 7))
    colors = plt.cm.tab10(np.linspace(0, 1, len(groups)))

    csv_lines = ['# V_3Q lattice data and Y-ansatz fit at beta=5.7\n']
    csv_lines.append('i,j,k,V_3Q_latt,error,V_3Q_fit\n')

    for idx, ((j, k), points) in enumerate(sorted(groups.items())):
        i_vals = [p[0] for p in points]
        v_vals = [p[1] for p in points]
        e_vals = [p[2] for p in points]

        # Plot data points
        ax.errorbar(i_vals, v_vals, yerr=e_vals, fmt='o', color=colors[idx],
                     markersize=6, capsize=3, label=r'$(j,k)=({},{})$'.format(j, k))

        # Continuous fit curve
        i_cont = np.linspace(0, max(max(i_vals), 3) + 0.5, 50)
        v_fit = []
        for ic in i_cont:
            v_fit.append(v3q_y_ansatz((ic, j, k), A, sigma, C))
        ax.plot(i_cont, v_fit, '-', color=colors[idx], linewidth=1.0, alpha=0.7)

        for p in points:
            vf = v3q_y_ansatz((p[0], j, k), A, sigma, C)
            csv_lines.append('{},{},{},{:.8f},{:.8f},{:.8f}\n'.format(
                p[0], j, k, p[1], p[2], vf))

    # Save data
    os.makedirs('../data', exist_ok=True)
    with open('../data/fig12_yansatz_fit_b57.csv', 'w') as f:
        f.writelines(csv_lines)

    ax.set_xlabel(r'$i$ (lattice units)', fontsize=14)
    ax.set_ylabel(r'$V_{3Q}$ (lattice units)', fontsize=14)
    ax.set_title(r'Fig. 12: $V_{3Q}$ data vs Y-ansatz fit ($\beta=5.7$)', fontsize=14)
    ax.legend(fontsize=9, ncol=2, loc='upper left')
    ax.grid(True, alpha=0.3)

    os.makedirs('../plots', exist_ok=True)
    plt.tight_layout()
    plt.savefig('../plots/fig12_yansatz_fit_b57.png', dpi=150)
    plt.close()
    print('Fig 12 saved.')


if __name__ == '__main__':
    main()
