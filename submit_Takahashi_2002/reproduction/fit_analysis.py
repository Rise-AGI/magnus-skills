"""
Fit Analysis: Comparison of Y-ansatz and Delta-ansatz for the 3Q potential.

Performs independent fits at beta=5.7 and 5.8, computing chi^2/NDF and
comparing with the paper's reported values.
"""
import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from three_quark_potential import (
    DATA_BETA_57, DATA_BETA_58, FIT_PARAMS, DELTA_FIT_PARAMS,
    fit_y_ansatz, fit_delta_ansatz, chi2_per_dof,
    compute_lmin, quark_distances_type1, v3q_y_ansatz, v3q_delta_ansatz
)


def main():
    csv_lines = ['# Fit analysis: Y-ansatz vs Delta-ansatz\n']
    csv_lines.append('beta,ansatz,sigma,A,C,chi2_NDF\n')

    for beta, data in [(5.7, DATA_BETA_57), (5.8, DATA_BETA_58)]:
        print(f'\n=== beta = {beta} ===')

        # Y-ansatz fit
        (A_y, sigma_y, C_y), _ = fit_y_ansatz(data)
        v_latt = np.array([d[1] for d in data])
        errors = np.array([d[2] for d in data])
        v_fit_y = np.array([v3q_y_ansatz(d[0], A_y, sigma_y, C_y) for d in data])
        chi2_y = chi2_per_dof(v_latt, v_fit_y, errors)

        print(f'Y-ansatz: sigma={sigma_y:.4f}, A={A_y:.4f}, C={C_y:.4f}, chi2/NDF={chi2_y:.2f}')
        print(f'  Paper:  sigma={FIT_PARAMS[beta]["3QY"]["sigma"]:.4f}, '
              f'A={FIT_PARAMS[beta]["3QY"]["A"]:.4f}, C={FIT_PARAMS[beta]["3QY"]["C"]:.4f}')

        csv_lines.append(f'{beta},Y-ansatz,{sigma_y:.6f},{A_y:.6f},{C_y:.6f},{chi2_y:.4f}\n')

        # Delta-ansatz fit
        (A_d, sigma_d, C_d), _ = fit_delta_ansatz(data)
        v_fit_d = np.array([v3q_delta_ansatz(d[0], A_d, sigma_d, C_d) for d in data])
        chi2_d = chi2_per_dof(v_latt, v_fit_d, errors)

        print(f'Delta:    sigma={sigma_d:.4f}, A={A_d:.4f}, C={C_d:.4f}, chi2/NDF={chi2_d:.2f}')
        csv_lines.append(f'{beta},Delta-ansatz,{sigma_d:.6f},{A_d:.6f},{C_d:.6f},{chi2_d:.4f}\n')

        # Key physics: string tension universality
        sigma_qq = FIT_PARAMS[beta]["QQ"]["sigma"]
        A_qq = FIT_PARAMS[beta]["QQ"]["A"]
        print(f'\nString tension: sigma_3Q/sigma_QQ = {sigma_y/sigma_qq:.3f} (expect ~1.0)')
        print(f'Coulomb coeff:  A_3Q/A_QQ = {A_y/A_qq:.3f} (expect ~0.5)')
        print(f'Delta sigma_D/sigma_QQ = {sigma_d/sigma_qq:.3f} (expect ~0.53)')

    # Save data
    os.makedirs('../data', exist_ok=True)
    with open('../data/fit_analysis_comparison.csv', 'w') as f:
        f.writelines(csv_lines)

    # Residual plot: Y-ansatz vs Delta-ansatz at beta=5.7
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    for ax_idx, (beta, data, label) in enumerate([(5.7, DATA_BETA_57, r'$\beta=5.7$'),
                                                    (5.8, DATA_BETA_58, r'$\beta=5.8$')]):
        (A_y, sigma_y, C_y), _ = fit_y_ansatz(data)
        (A_d, sigma_d, C_d), _ = fit_delta_ansatz(data)

        v_latt = np.array([d[1] for d in data])
        lmins = np.array([compute_lmin(*d[0]) for d in data])

        v_fit_y = np.array([v3q_y_ansatz(d[0], A_y, sigma_y, C_y) for d in data])
        v_fit_d = np.array([v3q_delta_ansatz(d[0], A_d, sigma_d, C_d) for d in data])

        res_y = v_latt - v_fit_y
        res_d = v_latt - v_fit_d

        idx = np.argsort(lmins)
        axes[ax_idx].plot(lmins[idx], res_y[idx], 'ro-', markersize=5,
                          label='Y-ansatz residuals')
        axes[ax_idx].plot(lmins[idx], res_d[idx], 'b^-', markersize=5,
                          label=r'$\Delta$-ansatz residuals')
        axes[ax_idx].axhline(0, color='k', linestyle='--', linewidth=0.5)
        axes[ax_idx].set_xlabel(r'$L_{min}$ (lattice units)', fontsize=12)
        axes[ax_idx].set_ylabel(r'$V_{3Q}^{latt} - V_{3Q}^{fit}$', fontsize=12)
        axes[ax_idx].set_title(label, fontsize=14)
        axes[ax_idx].legend(fontsize=10)
        axes[ax_idx].grid(True, alpha=0.3)

    os.makedirs('../plots', exist_ok=True)
    plt.tight_layout()
    plt.savefig('../plots/fit_analysis_residuals.png', dpi=150)
    plt.close()
    print('\nFit analysis and residual plots saved.')


if __name__ == '__main__':
    main()
