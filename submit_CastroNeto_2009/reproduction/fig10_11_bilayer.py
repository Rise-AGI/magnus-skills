"""
Figure 10/11: Bilayer graphene band structure.

Figure 10: V=0, gamma3=0
Figure 11: V!=0, gamma3=0
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(os.path.realpath(sys.argv[0]))))
from graphene_params import T_HOP, A_CC, GAMMA1

# Ensure output directories exist
os.makedirs("../data", exist_ok=True)
os.makedirs("../plots", exist_ok=True)


def bilayer_hamiltonian(kx, ky, V=0.0, gamma1=GAMMA1, gamma3=0.0, t=T_HOP):
    """Build 4x4 bilayer Hamiltonian at (kx, ky) near K point."""
    a = A_CC
    vF = 3 * t * a / 2  # eV*m
    k = kx + 1j * ky  # complex wavevector

    H = np.array([
        [-V,            vF * k,        0,                  3 * gamma3 * a * np.conj(k)],
        [vF * np.conj(k), -V,          gamma1,             0],
        [0,             gamma1,         V,                  vF * k],
        [3 * gamma3 * a * k, 0,        vF * np.conj(k),    V]
    ], dtype=complex)
    return H


def compute_bilayer_bands(k_array, V=0.0, gamma1=GAMMA1, gamma3=0.0):
    bands = np.zeros((len(k_array), 4))
    for i, kval in enumerate(k_array):
        H = bilayer_hamiltonian(kval, 0.0, V, gamma1, gamma3)
        evals = np.linalg.eigvalsh(H)
        bands[i, :] = np.sort(evals.real)
    return bands


def save_data(k_arr, bands_V0, bands_V):
    data_V0 = np.column_stack([k_arr, bands_V0])
    header_V0 = "k,E1,E2,E3,E4"
    np.savetxt("../data/fig10_bilayer_V0.csv", data_V0, delimiter=",",
               header=header_V0, fmt="%.8e", comments="")

    data_V = np.column_stack([k_arr, bands_V])
    header_V = "k,E1,E2,E3,E4"
    np.savetxt("../data/fig11_bilayer_V.csv", data_V, delimiter=",",
               header=header_V, fmt="%.8e", comments="")


def plot_figure(k_arr, bands_V0, bands_V):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    k_inv_nm = k_arr * 1e-9  # Convert to 1/nm for plot labeling (arbitrary)

    # Figure 10: V=0
    for i in range(4):
        ax1.plot(k_arr, bands_V0[:, i], 'b-', linewidth=1.0)
    ax1.set_xlabel(r'$k$ (1/m)')
    ax1.set_ylabel('E (eV)')
    ax1.set_title(r'Bilayer graphene: $V=0$, $\gamma_3=0$')
    ax1.axhline(0, color='gray', linewidth=0.5, linestyle='--')
    ax1.set_ylim(-1.0, 1.0)

    # Figure 11: V != 0
    for i in range(4):
        ax2.plot(k_arr, bands_V[:, i], 'b-', linewidth=1.0)
    ax2.set_xlabel(r'$k$ (1/m)')
    ax2.set_ylabel('E (eV)')
    ax2.set_title(r'Bilayer graphene: $V=0.1$ eV, $\gamma_3=0$')
    ax2.axhline(0, color='gray', linewidth=0.5, linestyle='--')
    ax2.set_ylim(-1.0, 1.0)

    plt.tight_layout()
    plt.savefig("../plots/fig10_11_bilayer.png", dpi=150, bbox_inches='tight')
    plt.close()


if True:
    # k range near K point
    k_range = np.linspace(-3e9, 3e9, 500)

    print("Computing bilayer bands (V=0)...")
    bands_V0 = compute_bilayer_bands(k_range, V=0.0)

    print("Computing bilayer bands (V=0.1 eV)...")
    bands_V = compute_bilayer_bands(k_range, V=0.1)

    save_data(k_range, bands_V0, bands_V)
    plot_figure(k_range, bands_V0, bands_V)
    print("Figure 10/11 done: bilayer band structure saved.")
