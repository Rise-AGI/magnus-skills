"""
Reproduce Figure 1 from Martinez-Fuentes & Herrera-Velazquez (2017).

Log-log plot of the breakdown distance rho_d vs N_D for tolerances of 1%, 5%, and 10%.
Includes the empirical fit log10(rho_d) = 1.18828 - 1.07493*log10(N_D) and
the curve log10(rho) = -log10(N_D) (i.e., N_D * rho = 1).

Outputs: ../data/fig1_breakdown.csv, ../plots/fig1_breakdown.png
"""

import sys
import os
import csv
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(".")))
from debye_solver import find_breakdown_distance

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    ND_values = [0.1, 0.3, 1, 2, 5, 10, 40, 1e2, 5e2, 1e3, 1e4, 1e5, 1e6, 1e7]

    rho_1 = []
    rho_5 = []
    rho_10 = []

    for N_D in ND_values:
        print(f"Computing N_D = {N_D:.1e}...")
        rho_1.append(find_breakdown_distance(N_D, 1.0))
        rho_5.append(find_breakdown_distance(N_D, 5.0))
        rho_10.append(find_breakdown_distance(N_D, 10.0))

    # Save data
    outdir = os.path.join(os.path.dirname(os.path.abspath(".")), "data")
    os.makedirs(outdir, exist_ok=True)
    outpath = os.path.join(outdir, "fig1_breakdown.csv")

    with open(outpath, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["N_D", "rho_d_1pct", "rho_d_5pct", "rho_d_10pct"])
        for i, N_D in enumerate(ND_values):
            writer.writerow([f"{N_D:.1e}", f"{rho_1[i]:.8e}", f"{rho_5[i]:.8e}", f"{rho_10[i]:.8e}"])

    print(f"Data saved to {outpath}")

    # Plot
    plotdir = os.path.join(os.path.dirname(os.path.abspath(".")), "plots")
    os.makedirs(plotdir, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8, 6))

    ax.loglog(ND_values, rho_1, 'bo', markersize=8, label=r'$\tau = 1\%$')
    ax.loglog(ND_values, rho_5, 'gs', markersize=8, label=r'$\tau = 5\%$')
    ax.loglog(ND_values, rho_10, 'rD', markersize=8, label=r'$\tau = 10\%$')

    # Empirical fit: log10(rho_d) = 1.18828 - 1.07493 * log10(N_D)
    # valid for N_D from 40 to 1e7
    ND_fit = np.logspace(np.log10(40), 7, 200)
    rho_fit = 10**(1.18828 - 1.07493 * np.log10(ND_fit))
    ax.loglog(ND_fit, rho_fit, 'b-', linewidth=1.5, label=r'Fit: $\log_{10}\rho_d = 1.188 - 1.075 \log_{10} N_D$')

    # Curve N_D * rho = 1
    ND_line = np.logspace(-1, 7, 200)
    rho_line = 1.0 / ND_line
    ax.loglog(ND_line, rho_line, 'r-', linewidth=1.5, label=r'$N_D \rho = 1$')

    ax.set_xlabel(r'$N_D$', fontsize=14)
    ax.set_ylabel(r'$\rho_d$', fontsize=14)
    ax.set_title('Breakdown distances of the Debye linearization', fontsize=14)
    ax.legend(fontsize=10, loc='upper right')
    ax.grid(True, which='both', alpha=0.3)
    ax.set_xlim([0.05, 2e7])

    plotpath = os.path.join(plotdir, "fig1_breakdown.png")
    fig.savefig(plotpath, dpi=150, bbox_inches='tight')
    print(f"Plot saved to {plotpath}")
    plt.close(fig)


if __name__ == "__main__":
    main()
