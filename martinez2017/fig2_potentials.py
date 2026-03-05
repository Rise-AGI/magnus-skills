"""
Plot comparisons of exact vs approximate Debye potentials for selected N_D values.

Shows the normalized potential Psi(rho) for both exact (numerical) and approximate
(Yukawa) solutions side by side for N_D = 1, 5, 40, 1000.

Outputs: ../data/fig2_potentials.csv, ../plots/fig2_potentials.png
"""

import sys
import os
import csv
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(".")))
from debye_solver import solve_exact, psi_approx

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    ND_cases = [1, 5, 40, 1000]
    rho_start_map = {1: 11, 5: 7, 40: 5, 1000: 2}

    outdir = os.path.join(os.path.dirname(os.path.abspath(".")), "data")
    os.makedirs(outdir, exist_ok=True)
    plotdir = os.path.join(os.path.dirname(os.path.abspath(".")), "plots")
    os.makedirs(plotdir, exist_ok=True)

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.flatten()

    all_data = {}

    for idx, N_D in enumerate(ND_cases):
        print(f"Solving for N_D = {N_D}...")
        rho_start = rho_start_map[N_D]
        n_steps = 200000

        rho_arr, psi_exact = solve_exact(N_D, rho_start, 1e-4, n_steps)
        psi_a = psi_approx(rho_arr)

        # Subsample for CSV
        step = max(1, len(rho_arr) // 500)
        rho_sub = rho_arr[::step]
        psi_ex_sub = psi_exact[::step]
        psi_a_sub = psi_a[::step]

        all_data[N_D] = (rho_sub, psi_ex_sub, psi_a_sub)

        ax = axes[idx]
        mask = rho_arr > 0.01
        ax.semilogy(rho_arr[mask], psi_exact[mask], 'b-', linewidth=2, label='Exact')
        ax.semilogy(rho_arr[mask], psi_a[mask], 'r--', linewidth=2, label='Approximate (Yukawa)')
        ax.set_xlabel(r'$\rho = r/\lambda_D$', fontsize=12)
        ax.set_ylabel(r'$\Psi(\rho)$', fontsize=12)
        ax.set_title(f'$N_D = {N_D}$', fontsize=13)
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)
        ax.set_xlim([0, rho_start])

    fig.suptitle('Exact vs Approximate Debye Shielding Potential', fontsize=15)
    fig.tight_layout()
    plotpath = os.path.join(plotdir, "fig2_potentials.png")
    fig.savefig(plotpath, dpi=150, bbox_inches='tight')
    print(f"Plot saved to {plotpath}")
    plt.close(fig)

    # Save CSV for one representative case (N_D=40)
    csvpath = os.path.join(outdir, "fig2_potentials.csv")
    with open(csvpath, "w", newline="") as f:
        writer = csv.writer(f)
        header = []
        for N_D in ND_cases:
            header.extend([f"rho_ND{N_D}", f"psi_exact_ND{N_D}", f"psi_approx_ND{N_D}"])
        writer.writerow(header)
        max_len = max(len(all_data[N_D][0]) for N_D in ND_cases)
        for i in range(max_len):
            row = []
            for N_D in ND_cases:
                rho_s, psi_e, psi_a = all_data[N_D]
                if i < len(rho_s):
                    row.extend([f"{rho_s[i]:.8e}", f"{psi_e[i]:.8e}", f"{psi_a[i]:.8e}"])
                else:
                    row.extend(["", "", ""])
            writer.writerow(row)

    print(f"Data saved to {csvpath}")


if __name__ == "__main__":
    main()
