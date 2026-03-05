"""
Figure 5: Density of states per unit cell.

Top: DOS for t'=0.2t
Bottom: DOS for t'=0
With zoom insets near neutrality point.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(os.path.realpath(sys.argv[0]))))
from graphene_params import A_CC, T_HOP, T_PRIME, dos_numerical

# Ensure output directories exist
os.makedirs("../data", exist_ok=True)
os.makedirs("../plots", exist_ok=True)


def compute_dos(tp_val, n_k=800, n_bins=600):
    a = A_CC
    t = T_HOP

    # Create k-point grid over the full Brillouin zone
    b1 = 2 * np.pi / (3 * a) * np.array([1.0, np.sqrt(3)])
    b2 = 2 * np.pi / (3 * a) * np.array([1.0, -np.sqrt(3)])

    n1 = np.linspace(0, 1, n_k, endpoint=False)
    n2 = np.linspace(0, 1, n_k, endpoint=False)
    N1, N2 = np.meshgrid(n1, n2)

    kx = N1 * b1[0] + N2 * b2[0]
    ky = N1 * b1[1] + N2 * b2[1]

    E_centers, dos_vals = dos_numerical(kx, ky, t, tp_val, n_bins, E_range=(-3.5, 3.5))

    return E_centers, dos_vals


def save_data(E_tp, dos_tp, E_0, dos_0):
    # Trim to same length
    n = min(len(E_tp), len(E_0))
    data = np.column_stack([E_tp[:n], dos_tp[:n], E_0[:n], dos_0[:n]])
    header = "E_over_t_tp02,DOS_tp02,E_over_t_tp0,DOS_tp0"
    np.savetxt("../data/fig5_dos.csv", data, delimiter=",",
               header=header, fmt="%.8e", comments="")


def plot_figure(E_tp, dos_tp, E_0, dos_0):
    fig, axes = plt.subplots(2, 2, figsize=(10, 8),
                              gridspec_kw={'width_ratios': [3, 1]})

    # Top left: DOS with t'=0.2t
    axes[0, 0].plot(E_tp, dos_tp, 'b-', linewidth=0.8)
    axes[0, 0].set_xlabel(r"$\epsilon/t$")
    axes[0, 0].set_ylabel(r"$\rho(\epsilon)$")
    axes[0, 0].set_title(r"$t' = 0.2t$")
    axes[0, 0].set_xlim(-4, 2)
    axes[0, 0].set_ylim(0, 5)

    # Top right: zoom near neutrality
    mask_tp = (E_tp > -1.0) & (E_tp < 1.0)
    axes[0, 1].plot(E_tp[mask_tp], dos_tp[mask_tp], 'b-', linewidth=0.8)
    axes[0, 1].set_xlabel(r"$\epsilon/t$")
    axes[0, 1].set_title("Zoom near E=0")

    # Bottom left: DOS with t'=0
    axes[1, 0].plot(E_0, dos_0, 'b-', linewidth=0.8)
    axes[1, 0].set_xlabel(r"$\epsilon/t$")
    axes[1, 0].set_ylabel(r"$\rho(\epsilon)$")
    axes[1, 0].set_title(r"$t' = 0$")
    axes[1, 0].set_xlim(-4, 4)
    axes[1, 0].set_ylim(0, 5)

    # Bottom right: zoom near neutrality
    mask_0 = (E_0 > -1.0) & (E_0 < 1.0)
    axes[1, 1].plot(E_0[mask_0], dos_0[mask_0], 'b-', linewidth=0.8)
    axes[1, 1].set_xlabel(r"$\epsilon/t$")
    axes[1, 1].set_title("Zoom near E=0")

    plt.tight_layout()
    plt.savefig("../plots/fig5_dos.png", dpi=150, bbox_inches='tight')
    plt.close()


if True:
    print("Computing DOS for t'=0.2t...")
    E_tp, dos_tp = compute_dos(T_PRIME)
    print("Computing DOS for t'=0...")
    E_0, dos_0 = compute_dos(0.0)
    save_data(E_tp, dos_tp, E_0, dos_0)
    plot_figure(E_tp, dos_tp, E_0, dos_0)
    print("Figure 5 done: DOS saved.")
