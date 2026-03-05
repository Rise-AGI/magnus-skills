"""
Figure 8: V_3Q vs L_min at beta=5.8.

Plots the three-quark potential as a function of the minimal flux-tube length
at beta=5.8, showing linearity that supports the Y-ansatz.
"""
import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from three_quark_potential import (DATA_BETA_58, compute_lmin, FIT_PARAMS)


def main():
    data = DATA_BETA_58
    params = FIT_PARAMS[5.8]["3QY"]

    lmins = []
    v_latt = []
    errors = []
    for entry in data:
        cfg, v, err = entry[0], entry[1], entry[2]
        lmins.append(compute_lmin(*cfg))
        v_latt.append(v)
        errors.append(err)

    lmins = np.array(lmins)
    v_latt = np.array(v_latt)
    errors = np.array(errors)

    idx = np.argsort(lmins)
    lmins = lmins[idx]
    v_latt = v_latt[idx]
    errors = errors[idx]

    lmin_range = np.linspace(0, max(lmins) * 1.1, 100)

    os.makedirs('../data', exist_ok=True)
    with open('../data/fig8_v3q_lmin_b58.csv', 'w') as f:
        f.write('# V_3Q vs L_min at beta=5.8\n')
        f.write('L_min,V_3Q_latt,error\n')
        for i in range(len(lmins)):
            f.write('{:.8f},{:.8f},{:.8f}\n'.format(lmins[i], v_latt[i], errors[i]))

    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    ax.errorbar(lmins, v_latt, yerr=errors, fmt='bo', markersize=4, capsize=2,
                label=r'$V_{3Q}^{latt}$ ($\beta=5.8$)', alpha=0.7)

    sigma = params["sigma"]
    C = params["C"]
    ax.plot(lmin_range, sigma * lmin_range + C, 'r--', linewidth=1.5,
            label=r'$\sigma_{3Q} L_{min} + C_{3Q}$')

    ax.set_xlabel(r'$L_{min}$ (lattice units)', fontsize=14)
    ax.set_ylabel(r'$V_{3Q}$ (lattice units)', fontsize=14)
    ax.set_title(r'Fig. 8: $V_{3Q}$ vs $L_{min}$ at $\beta=5.8$', fontsize=14)
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)

    os.makedirs('../plots', exist_ok=True)
    plt.tight_layout()
    plt.savefig('../plots/fig8_v3q_lmin_b58.png', dpi=150)
    plt.close()
    print('Fig 8 saved.')


if __name__ == '__main__':
    main()
