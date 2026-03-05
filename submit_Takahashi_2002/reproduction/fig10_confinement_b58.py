"""
Figure 10: Semi-quantitative test on confinement in 3Q potential at beta=5.8.

Plots V_3Q - V_3Q^Coul vs L_min, demonstrating linearity after subtracting
the Coulomb contribution (using A_QQ/2 from the Q-Qbar potential).
"""
import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from three_quark_potential import (DATA_BETA_58, compute_lmin, FIT_PARAMS,
                                   quark_distances_type1)


def main():
    data = DATA_BETA_58
    params_qq = FIT_PARAMS[5.8]["QQ"]
    A_qq = params_qq["A"]
    A_3Q_coul = A_qq / 2.0  # OGE relation

    lmins = []
    v_sub = []
    errors = []

    for entry in data:
        cfg = entry[0]
        v_latt = entry[1]
        err = entry[2]
        i, j, k = cfg

        # Compute Coulomb subtraction
        d12, d13, d23 = quark_distances_type1(i, j, k)
        coul = 0
        for d in [d12, d13, d23]:
            if d > 1e-10:
                coul += 1.0 / d
        v_coul = -A_3Q_coul * coul

        lmin = compute_lmin(i, j, k)
        lmins.append(lmin)
        v_sub.append(v_latt - v_coul)
        errors.append(err)

    lmins = np.array(lmins)
    v_sub = np.array(v_sub)
    errors = np.array(errors)

    # Sort
    idx = np.argsort(lmins)
    lmins = lmins[idx]
    v_sub = v_sub[idx]
    errors = errors[idx]

    # Linear fit
    from numpy.polynomial import polynomial as P
    coeffs = np.polyfit(lmins, v_sub, 1)
    lmin_range = np.linspace(0, max(lmins) * 1.1, 100)

    # Save data
    os.makedirs('../data', exist_ok=True)
    with open('../data/fig10_confinement_b58.csv', 'w') as f:
        f.write('# V_3Q - V_3Q^Coul vs L_min at beta=5.8\n')
        f.write('# A_QQ/2 = {:.6f} used for Coulomb subtraction\n'.format(A_3Q_coul))
        f.write('L_min,V_3Q_minus_Coul,error\n')
        for i in range(len(lmins)):
            f.write('{:.8f},{:.8f},{:.8f}\n'.format(lmins[i], v_sub[i], errors[i]))

    # Plot
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    ax.errorbar(lmins, v_sub, yerr=errors, fmt='bo', markersize=4, capsize=2,
                label=r'$V_{3Q} - V_{3Q}^{Coul}$', alpha=0.7)
    ax.plot(lmin_range, np.polyval(coeffs, lmin_range), 'r-', linewidth=1.5,
            label=r'Linear fit: $\sigma={:.4f}$'.format(coeffs[0]))

    ax.set_xlabel(r'$L_{min}$ (lattice units)', fontsize=14)
    ax.set_ylabel(r'$V_{3Q} - V_{3Q}^{Coul}$ (lattice units)', fontsize=14)
    ax.set_title(r'Fig. 10: Confinement Test at $\beta=5.8$', fontsize=14)
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)

    os.makedirs('../plots', exist_ok=True)
    plt.tight_layout()
    plt.savefig('../plots/fig10_confinement_b58.png', dpi=150)
    plt.close()
    print('Fig 10 saved.')


if __name__ == '__main__':
    main()
