"""
Figure 2: RI/MSbar matching coefficient and RG evolution coefficient.

Computes:
- Delta_Z^{RI/MSbar}(mu) — the matching coefficient (Eq. B.1)
- c_S^MSbar(mu) — the RG evolution coefficient (Eq. B.6)
Both at NNLO for quenched QCD.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, '.')
from qcd_perturbative import (
    delta_z_ri_msbar, evolution_coefficient,
    ri_msbar_matching_coefficients, alpha_s_value,
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


def main():
    print("=" * 60)
    print("Figure 2: RI/MSbar matching and RG evolution")
    print("=" * 60)

    c1, c2 = ri_msbar_matching_coefficients()
    print(f"\nRI/MSbar matching coefficients (N_c={N_C}, N_f={N_F}):")
    print(f"  C^(1) = {c1:.4f}")
    print(f"  C^(2) = {c2:.4f}")

    mu_values = np.linspace(1.0, 10.0, 400)

    dz_values = np.array([delta_z_ri_msbar(mu) for mu in mu_values])
    cs_values = np.array([evolution_coefficient(mu) for mu in mu_values])

    # Print at key scales
    print("\nDelta_Z^{RI/MSbar}(mu):")
    for mu_ref in [1.0, 2.0, 3.0, 4.0]:
        dz = delta_z_ri_msbar(mu_ref)
        cs = evolution_coefficient(mu_ref)
        print(f"  mu={mu_ref:.0f} GeV: Delta_Z = {dz:.6f}, c_S = {cs:.6f}")

    # Ratio c_S(2 GeV)/c_S(mu) used for RG running
    cs_2gev = evolution_coefficient(2.0)
    ratio_values = cs_2gev / cs_values

    # Plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # Left: Matching coefficient
    ax1.plot(mu_values, dz_values, 'b-', linewidth=2)
    for mu_lat, label in [(2.12, r'$a^{-1}_{6.0}$'), (2.7, r'$a^{-1}_{6.2}$'),
                           (4.0, r'$a^{-1}_{6.4}$')]:
        dz_lat = delta_z_ri_msbar(mu_lat)
        ax1.plot(mu_lat, dz_lat, 'ro', markersize=8)
        ax1.annotate(label, (mu_lat, dz_lat), textcoords="offset points",
                    xytext=(8, 5), fontsize=10)

    ax1.set_xlabel(r'$\mu$ [GeV]')
    ax1.set_ylabel(r'$\Delta Z^{RI/\overline{MS}}(\mu)$')
    ax1.set_title(r'RI/$\overline{MS}$ Matching Coefficient')
    ax1.grid(True, linestyle='--', alpha=0.3)

    # Right: Evolution ratio
    ax2.plot(mu_values, ratio_values, 'r-', linewidth=2,
             label=r'$c_S^{\overline{MS}}(2\,\mathrm{GeV}) / c_S^{\overline{MS}}(\mu)$')
    for mu_lat, label in [(2.12, r'$a^{-1}_{6.0}$'), (2.7, r'$a^{-1}_{6.2}$'),
                           (4.0, r'$a^{-1}_{6.4}$')]:
        ratio_lat = cs_2gev / evolution_coefficient(mu_lat)
        ax2.plot(mu_lat, ratio_lat, 'ro', markersize=8)
        ax2.annotate(label, (mu_lat, ratio_lat), textcoords="offset points",
                    xytext=(8, 5), fontsize=10)

    ax2.set_xlabel(r'$\mu$ [GeV]')
    ax2.set_ylabel(r'RG evolution ratio')
    ax2.set_title(r'RG Running to $\mu = 2$ GeV')
    ax2.legend(loc='best')
    ax2.grid(True, linestyle='--', alpha=0.3)

    plt.tight_layout()
    plt.savefig('../plots/fig2_matching_evolution.png', bbox_inches='tight')
    print("\nFigure saved to plots/fig2_matching_evolution.png")
    plt.close()

    # Export data
    with open('../data/fig2_matching_evolution.csv', 'w') as f:
        f.write("# Figure 2: RI/MSbar matching and RG evolution at NNLO\n")
        f.write("# Quenched QCD: N_c=3, N_f=0\n")
        f.write("# Columns: mu_GeV, delta_z_ri_msbar, c_s_msbar, rg_ratio_to_2GeV\n")
        f.write("mu_GeV,delta_z_ri_msbar,c_s_msbar,rg_ratio_to_2GeV\n")
        for i in range(len(mu_values)):
            f.write(f"{mu_values[i]:.8f},{dz_values[i]:.8f},"
                    f"{cs_values[i]:.8f},{ratio_values[i]:.8f}\n")
    print("Data exported to data/fig2_matching_evolution.csv")


main()
