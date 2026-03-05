"""
Figure 4: RG running of the chiral condensate and scale dependence.

Computes how the chiral condensate value changes with the
renormalization scale mu, comparing NNLO RG running to LO and NLO.
Also shows the perturbative vs non-perturbative renormalization comparison.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, '.')
from qcd_perturbative import (
    alpha_s_value, alpha_s_nnlo, beta_coefficients, gamma_m_coefficients,
    evolution_coefficient, run_condensate,
    LAMBDA_QCD_MSBAR, N_C, N_F
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


def evolution_coefficient_lo(mu):
    """LO evolution coefficient (leading power only)."""
    beta_0, _, _ = beta_coefficients()
    gamma_0, _, _ = gamma_m_coefficients()
    gamma_s0_bar = -gamma_0 / (2.0 * beta_0)
    a_s = alpha_s_value(mu)
    return a_s**gamma_s0_bar


def evolution_coefficient_nlo(mu):
    """NLO evolution coefficient."""
    beta_0, beta_1, _ = beta_coefficients()
    gamma_0, gamma_1, _ = gamma_m_coefficients()
    beta_1_bar = beta_1 / beta_0
    gamma_s0_bar = -gamma_0 / (2.0 * beta_0)
    gamma_s1_bar = -gamma_1 / (2.0 * beta_0)

    a_s = alpha_s_value(mu)
    a_s_4pi = a_s / (4.0 * np.pi)

    c_s = a_s**gamma_s0_bar * (1.0 + a_s_4pi * (gamma_s1_bar - beta_1_bar * gamma_s0_bar))
    return c_s


def main():
    print("=" * 60)
    print("Figure 4: Scale dependence of the chiral condensate")
    print("=" * 60)

    # Best estimate at 2 GeV from paper
    condensate_2gev = 0.0147  # GeV^3

    mu_values = np.linspace(1.0, 10.0, 400)

    # NNLO running
    cs_2gev = evolution_coefficient(2.0)
    condensate_nnlo = np.array([
        condensate_2gev * evolution_coefficient(mu) / cs_2gev
        for mu in mu_values
    ])

    # NLO running
    cs_2gev_nlo = evolution_coefficient_nlo(2.0)
    condensate_nlo = np.array([
        condensate_2gev * evolution_coefficient_nlo(mu) / cs_2gev_nlo
        for mu in mu_values
    ])

    # LO running
    cs_2gev_lo = evolution_coefficient_lo(2.0)
    condensate_lo = np.array([
        condensate_2gev * evolution_coefficient_lo(mu) / cs_2gev_lo
        for mu in mu_values
    ])

    # Print at key scales
    print(f"\nChiral condensate -<psibar psi>^MSbar(mu) [GeV^3]:")
    print(f"  Input: -<psibar psi>(2 GeV) = {condensate_2gev} GeV^3")
    for mu_ref in [1.0, 2.0, 3.0, 4.0, 5.0]:
        cond = run_condensate(condensate_2gev, 2.0, mu_ref)
        cond_mev = (cond * 1e9)**(1.0/3.0)  # MeV
        print(f"  mu = {mu_ref:.0f} GeV: {cond:.6f} GeV^3  =  [{cond_mev:.1f} MeV]^3")

    # Paper result at 1 GeV
    cond_1gev_nnlo = run_condensate(condensate_2gev, 2.0, 1.0)
    print(f"\n  NNLO at 1 GeV: {cond_1gev_nnlo:.4f} GeV^3")
    print(f"  Paper at 1 GeV: 0.0124 GeV^3")

    # RGI condensate (Eq. 14)
    cond_rgi = condensate_2gev / cs_2gev
    print(f"\n  RGI condensate: {cond_rgi:.4f} GeV^3")
    print(f"  Paper RGI: 0.0088 GeV^3")

    # Cube root in MeV
    def cube_root_mev(x):
        return (x * 1e9)**(1.0/3.0) * 1e-3  # GeV^3 -> MeV

    # Plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # Left: Scale dependence
    ax1.plot(mu_values, condensate_nnlo * 1000, 'b-', linewidth=2, label='NNLO')
    ax1.plot(mu_values, condensate_nlo * 1000, 'r--', linewidth=1.5, label='NLO')
    ax1.plot(mu_values, condensate_lo * 1000, 'g:', linewidth=1.5, label='LO')

    # Mark paper values
    ax1.plot(2.0, condensate_2gev * 1000, 'ko', markersize=10, zorder=5,
             label=r'Paper: $\mu=2$ GeV')
    ax1.plot(1.0, 0.0124 * 1000, 'ks', markersize=8, zorder=5,
             label=r'Paper: $\mu=1$ GeV')

    ax1.set_xlabel(r'$\mu$ [GeV]')
    ax1.set_ylabel(r'$-\langle\bar{\psi}\psi\rangle^{\overline{MS}}(\mu)$ [$10^{-3}$ GeV$^3$]')
    ax1.set_title('Scale dependence of the chiral condensate')
    ax1.legend(loc='upper right')
    ax1.grid(True, linestyle='--', alpha=0.3)

    # Right: Cube root (in MeV) vs scale
    cond_mev_nnlo = np.array([(c * 1e9)**(1./3.) for c in condensate_nnlo])
    cond_mev_nlo = np.array([(c * 1e9)**(1./3.) for c in condensate_nlo])

    ax2.plot(mu_values, cond_mev_nnlo, 'b-', linewidth=2, label='NNLO')
    ax2.plot(mu_values, cond_mev_nlo, 'r--', linewidth=1.5, label='NLO')

    # Paper values
    ax2.plot(2.0, (0.0147e9)**(1./3.), 'ko', markersize=10, zorder=5,
             label=r'Paper: 245 MeV at $\mu=2$ GeV')
    ax2.plot(1.0, (0.0124e9)**(1./3.), 'ks', markersize=8, zorder=5,
             label=r'Paper: 231 MeV at $\mu=1$ GeV')

    # Sum rule comparison
    ax2.axhspan((0.012e9)**(1./3.), (0.016e9)**(1./3.), alpha=0.15, color='orange',
                label='Sum rule range')

    ax2.set_xlabel(r'$\mu$ [GeV]')
    ax2.set_ylabel(r'$[-\langle\bar{\psi}\psi\rangle^{\overline{MS}}(\mu)]^{1/3}$ [MeV]')
    ax2.set_title(r'Cube root of condensate vs. scale')
    ax2.legend(loc='upper right', fontsize=9)
    ax2.grid(True, linestyle='--', alpha=0.3)

    plt.tight_layout()
    plt.savefig('../plots/fig4_scale_dependence.png', bbox_inches='tight')
    print("\nFigure saved to plots/fig4_scale_dependence.png")
    plt.close()

    # Export data
    with open('../data/fig4_scale_dependence.csv', 'w') as f:
        f.write("# Figure 4: Scale dependence of the chiral condensate\n")
        f.write("# Input: -<psibar psi>^MSbar(2 GeV) = 0.0147 GeV^3\n")
        f.write("# Columns: mu_GeV, condensate_lo, condensate_nlo, condensate_nnlo\n")
        f.write("mu_GeV,condensate_lo_GeV3,condensate_nlo_GeV3,condensate_nnlo_GeV3\n")
        for i in range(len(mu_values)):
            f.write(f"{mu_values[i]:.8f},{condensate_lo[i]:.8f},"
                    f"{condensate_nlo[i]:.8f},{condensate_nnlo[i]:.8f}\n")
    print("Data exported to data/fig4_scale_dependence.csv")


main()
