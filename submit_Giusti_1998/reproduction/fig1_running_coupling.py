"""
Figure 1: Running coupling alpha_s(mu) at NNLO in quenched QCD.

Computes the running coupling constant alpha_s^MSbar(mu) at NNLO
for quenched QCD (N_f=0, N_c=3) using Eq. B.7 of the paper.
Also shows the sensitivity to Lambda_QCD uncertainty.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, '.')
from qcd_perturbative import (
    alpha_s_value, beta_coefficients, gamma_m_coefficients,
    LAMBDA_QCD_MSBAR, LAMBDA_QCD_MSBAR_ERR, N_C, N_F
)

plt.rcParams.update({
    'font.size': 12,
    'axes.labelsize': 14,
    'axes.titlesize': 14,
    'legend.fontsize': 10,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'text.usetex': False,
    'mathtext.fontset': 'cm',
})


def main():
    print("=" * 60)
    print("Figure 1: Running coupling alpha_s(mu) at NNLO")
    print("=" * 60)

    # Print beta-function coefficients
    b0, b1, b2 = beta_coefficients()
    print(f"\nBeta-function coefficients (N_c={N_C}, N_f={N_F}):")
    print(f"  beta_0 = {b0:.4f}")
    print(f"  beta_1 = {b1:.4f}")
    print(f"  beta_2 = {b2:.4f}")

    g0, g1, g2 = gamma_m_coefficients()
    print(f"\nAnomalous dimension coefficients:")
    print(f"  gamma_m^(0) = {g0:.4f}")
    print(f"  gamma_m^(1) = {g1:.4f}")
    print(f"  gamma_m^(2) = {g2:.4f}")

    # Compute alpha_s over range of scales
    mu_values = np.linspace(0.8, 10.0, 500)

    # Central value
    alpha_central = np.array([alpha_s_value(mu) for mu in mu_values])

    # Error band from Lambda_QCD uncertainty
    lam_hi = LAMBDA_QCD_MSBAR + LAMBDA_QCD_MSBAR_ERR
    lam_lo = LAMBDA_QCD_MSBAR - LAMBDA_QCD_MSBAR_ERR
    alpha_hi = np.array([alpha_s_value(mu, lam_hi) for mu in mu_values])
    alpha_lo = np.array([alpha_s_value(mu, lam_lo) for mu in mu_values])

    # Print values at key scales
    print(f"\nLambda_QCD^MSbar = {LAMBDA_QCD_MSBAR} +/- {LAMBDA_QCD_MSBAR_ERR} GeV")
    for mu_ref in [1.0, 2.0, 3.0, 4.0]:
        a_s = alpha_s_value(mu_ref)
        a_s_hi = alpha_s_value(mu_ref, lam_hi)
        a_s_lo = alpha_s_value(mu_ref, lam_lo)
        print(f"  alpha_s({mu_ref:.0f} GeV) = {a_s:.4f} [{a_s_lo:.4f}, {a_s_hi:.4f}]")

    # Plot
    fig, ax = plt.subplots(figsize=(8, 6))

    ax.fill_between(mu_values, alpha_lo, alpha_hi, alpha=0.3, color='blue',
                    label=r'$\Lambda_{QCD}$ uncertainty')
    ax.plot(mu_values, alpha_central, 'b-', linewidth=2,
            label=r'$\alpha_s^{\overline{MS}}(\mu)$ NNLO, $N_f=0$')

    # Mark lattice scales a^{-1} from Table 1
    lattice_scales = {
        r'$\beta=6.0$ (C)': 2.12,
        r'$\beta=6.2$ (C)': 2.7,
        r'$\beta=6.4$ (C)': 4.0,
    }
    for label, mu_lat in lattice_scales.items():
        a_s_lat = alpha_s_value(mu_lat)
        ax.plot(mu_lat, a_s_lat, 'ro', markersize=8)
        ax.annotate(label, (mu_lat, a_s_lat), textcoords="offset points",
                   xytext=(10, 5), fontsize=9)

    # Reference line at mu = 2 GeV
    ax.axvline(x=2.0, color='gray', linestyle='--', alpha=0.5)
    ax.text(2.05, 0.05, r'$\mu = 2$ GeV', fontsize=10, color='gray')

    ax.set_xlabel(r'$\mu$ [GeV]')
    ax.set_ylabel(r'$\alpha_s^{\overline{MS}}(\mu)$')
    ax.set_title(r'Running coupling at NNLO (quenched, $N_f = 0$)')
    ax.set_xlim(0.8, 10.0)
    ax.set_ylim(0.0, 0.8)
    ax.legend(loc='upper right')
    ax.grid(True, linestyle='--', alpha=0.3)

    plt.tight_layout()
    plt.savefig('../plots/fig1_running_coupling.png', bbox_inches='tight')
    print("\nFigure saved to plots/fig1_running_coupling.png")
    plt.close()

    # Export data
    with open('../data/fig1_running_coupling.csv', 'w') as f:
        f.write("# Figure 1: Running coupling alpha_s^MSbar(mu) at NNLO\n")
        f.write("# Quenched QCD: N_c=3, N_f=0, Lambda_QCD=0.251 GeV\n")
        f.write("# Columns: mu_GeV, alpha_s_central, alpha_s_low, alpha_s_high\n")
        f.write("mu_GeV,alpha_s_central,alpha_s_low,alpha_s_high\n")
        for i in range(len(mu_values)):
            f.write(f"{mu_values[i]:.8f},{alpha_central[i]:.8f},"
                    f"{alpha_lo[i]:.8f},{alpha_hi[i]:.8f}\n")
    print("Data exported to data/fig1_running_coupling.csv")


main()
