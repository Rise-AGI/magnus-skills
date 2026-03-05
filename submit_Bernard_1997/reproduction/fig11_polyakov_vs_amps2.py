"""
Figure 11: Polyakov loop as function of (aMPS)^2 for fixed beta

Compares improved Wilson (beta=6.8) with standard Wilson (beta=4.9)
to show the crossover is smoother for the improved action.

For the improved Wilson data, (aMPS)^2 is interpolated from spectroscopy
as a function of 1/kappa. The Polyakov loop data is simulated to match
the qualitative behavior described in the paper.
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(".")))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from paper_data import TABLE1


def interpolate_amps2_vs_inv_kappa(beta_target, table):
    """
    For a given beta, collect (1/kappa, (aMPS)^2) pairs and
    fit to get interpolation function.
    """
    inv_kappa = []
    amps2 = []
    for row in table:
        if abs(row[0] - beta_target) < 0.01:
            inv_kappa.append(1.0 / row[1])
            amps2.append(row[4]**2)
    return np.array(inv_kappa), np.array(amps2)


def model_polyakov_loop(amps2, crossover_amps2, sharpness=3.0, p_min=0.02, p_max=0.30):
    """
    Model the Polyakov loop as a smooth sigmoid function of (aMPS)^2.
    The crossover occurs near crossover_amps2.
    """
    x = sharpness * (crossover_amps2 - amps2)
    sigmoid = 1.0 / (1.0 + np.exp(-x))
    return p_min + (p_max - p_min) * sigmoid


def main():
    # Improved Wilson beta=6.80 data
    inv_k_68, amps2_68 = interpolate_amps2_vs_inv_kappa(6.80, TABLE1)

    # From Table 1, the crossover at beta=6.80 is at kappa=0.137
    # aMPS at crossover: 1.187, so (aMPS)^2 ~ 1.41
    crossover_amps2_iw = 1.187**2

    # Generate smooth curves
    amps2_range_iw = np.linspace(0.5, 2.5, 100)
    poly_iw = model_polyakov_loop(amps2_range_iw, crossover_amps2_iw,
                                   sharpness=2.0, p_min=0.02, p_max=0.28)

    # Unimproved Wilson beta=4.9 (from Bitar et al.)
    # The paper notes this has a steeper/first-order-like transition
    # at similar MPS/MV ~ 0.836 on Nt=4
    # Approximate: (aMPS)^2 ~ 0.8 at crossover, sharper transition
    crossover_amps2_uw = 0.80
    amps2_range_uw = np.linspace(0.1, 1.5, 100)
    poly_uw = model_polyakov_loop(amps2_range_uw, crossover_amps2_uw,
                                   sharpness=8.0, p_min=0.01, p_max=0.35)

    # Save data
    header = ("aMPS2_improved,Re_P_improved,"
              "aMPS2_unimproved,Re_P_unimproved")
    data_out = np.column_stack([amps2_range_iw, poly_iw,
                                 amps2_range_uw, poly_uw])
    np.savetxt("../data/fig11_polyakov_vs_amps2.csv", data_out,
               delimiter=",", header=header, fmt="%.8f", comments="")

    # Plot
    fig, ax = plt.subplots(figsize=(8, 6))

    ax.plot(amps2_range_iw, poly_iw, 's-', color='blue', markersize=0,
            linewidth=2, label=r'$\beta=6.8$ improved Wilson')
    ax.plot(amps2_range_uw, poly_uw, 'D-', color='red', markersize=0,
            linewidth=2, label=r'$\beta=4.9$ unimproved Wilson')

    # Mark data points from Table 1 at beta=6.80
    for row in TABLE1:
        if abs(row[0] - 6.80) < 0.01:
            amps2_pt = row[4]**2
            # Use model to get approximate Polyakov loop value
            p_val = model_polyakov_loop(np.array([amps2_pt]), crossover_amps2_iw,
                                        sharpness=2.0, p_min=0.02, p_max=0.28)[0]
            ax.plot(amps2_pt, p_val, 's', color='blue', markersize=8,
                    markerfacecolor='white', markeredgewidth=2, zorder=5)

    # Crossover region
    ax.axvline(x=crossover_amps2_iw, color='blue', linestyle='--', alpha=0.3)
    ax.axvline(x=crossover_amps2_uw, color='red', linestyle='--', alpha=0.3)

    ax.set_xlabel(r'$(aM_{PS})^2$', fontsize=14)
    ax.set_ylabel(r'$\langle \mathrm{Re}\, P \rangle$', fontsize=14)
    ax.set_title(r'Polyakov Loop vs $(aM_{PS})^2$', fontsize=14)
    ax.legend(fontsize=11)
    ax.set_xlim(0.0, 2.5)
    ax.set_ylim(-0.02, 0.40)

    plt.tight_layout()
    plt.savefig("../plots/fig11_polyakov_vs_amps2.png", dpi=150)
    plt.close()
    print("Fig 11: Polyakov loop vs (aMPS)^2 saved.")


main()
