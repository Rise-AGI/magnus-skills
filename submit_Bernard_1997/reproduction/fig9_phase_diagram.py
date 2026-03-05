"""
Figure 9: Phase diagram for the Symanzik-improved action.

Plots kappa_T(beta) (thermal crossover) and kappa_c(beta) (chiral limit)
in the (beta, kappa) plane for Nt=4 improved Wilson fermions.
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(".")))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from paper_data import (TABLE1, KAPPA_T, KAPPA_C, get_crossover_data)


def compute_phase_diagram():
    """Compute phase diagram data from Table 1 and crossover/critical lines."""

    # Crossover points (octagons in paper)
    beta_T = sorted(KAPPA_T.keys())
    kappa_T = [KAPPA_T[b] for b in beta_T]

    # Critical points (diamonds in paper)
    beta_C = sorted(KAPPA_C.keys())
    kappa_C = [KAPPA_C[b] for b in beta_C]

    # Zero temperature simulation points (crosses in paper)
    beta_0T = []
    kappa_0T = []
    for row in TABLE1:
        beta_0T.append(row[0])
        kappa_0T.append(row[1])

    return (np.array(beta_T), np.array(kappa_T),
            np.array(beta_C), np.array(kappa_C),
            np.array(beta_0T), np.array(kappa_0T))


def interpolate_line(beta_pts, kappa_pts, n_interp=200):
    """Interpolate a smooth curve through the phase boundary points."""
    beta_fine = np.linspace(min(beta_pts), max(beta_pts), n_interp)
    # Polynomial fit
    coeffs = np.polyfit(beta_pts, kappa_pts, min(3, len(beta_pts) - 1))
    kappa_fine = np.polyval(coeffs, beta_fine)
    return beta_fine, kappa_fine


def main():
    beta_T, kappa_T, beta_C, kappa_C, beta_0T, kappa_0T = compute_phase_diagram()

    # Interpolate smooth lines
    beta_T_fine, kappa_T_fine = interpolate_line(beta_T, kappa_T)
    beta_C_fine, kappa_C_fine = interpolate_line(beta_C, kappa_C)

    # Save data
    header = "beta_T_interp,kappa_T_interp,beta_C_interp,kappa_C_interp"
    data_out = np.column_stack([beta_T_fine, kappa_T_fine,
                                 beta_C_fine, kappa_C_fine])
    np.savetxt("../data/fig9_phase_diagram.csv", data_out,
               delimiter=",", header=header, fmt="%.8f", comments="")

    # Also save the discrete points
    with open("../data/fig9_phase_diagram.csv", "a") as f:
        f.write("\n# Crossover points (beta, kappa_T)\n")
        for b, k in zip(beta_T, kappa_T):
            f.write(f"# {b:.2f}, {k:.4f}\n")
        f.write("# Critical points (beta, kappa_c)\n")
        for b, k in zip(beta_C, kappa_C):
            f.write(f"# {b:.2f}, {k:.4f}\n")
        f.write("# Zero-T simulation points (beta, kappa)\n")
        for b, k in zip(beta_0T, kappa_0T):
            f.write(f"# {b:.2f}, {k:.4f}\n")

    # Plot
    fig, ax = plt.subplots(figsize=(8, 6))

    # Crossover line
    ax.plot(beta_T_fine, kappa_T_fine, 'k--', linewidth=1, label=r'$\kappa_T(\beta)$ (guide)')
    ax.plot(beta_T, kappa_T, 'o', color='blue', markersize=10, markerfacecolor='white',
            markeredgewidth=2, label=r'$N_t=4$ crossover', zorder=5)

    # Critical line
    ax.plot(beta_C_fine, kappa_C_fine, 'k:', linewidth=1, label=r'$\kappa_c(\beta)$ (guide)')
    ax.plot(beta_C, kappa_C, 'D', color='red', markersize=8, markerfacecolor='white',
            markeredgewidth=2, label=r'$\kappa_c$ estimates', zorder=5)

    # Zero-T points
    ax.plot(beta_0T, kappa_0T, 'x', color='green', markersize=8, markeredgewidth=2,
            label=r'$T=0$ simulations', zorder=5)

    ax.set_xlabel(r'$\beta$', fontsize=14)
    ax.set_ylabel(r'$\kappa$', fontsize=14)
    ax.set_title('Phase Diagram: Symanzik-Improved Action', fontsize=14)
    ax.legend(fontsize=11, loc='upper right')
    ax.set_xlim(6.2, 7.5)
    ax.set_ylim(0.11, 0.155)

    plt.tight_layout()
    plt.savefig("../plots/fig9_phase_diagram.png", dpi=150)
    plt.close()
    print("Fig 9: Phase diagram saved.")


main()
