"""
Figure 17: Energy spectrum of graphene nanoribbons.

Top: armchair nanoribbon (N=200 unit cells)
Bottom: zigzag nanoribbon (N=200 unit cells)
With zoom panels on the right.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(os.path.realpath(sys.argv[0]))))
from graphene_params import T_HOP, A_CC

# Ensure output directories exist
os.makedirs("../data", exist_ok=True)
os.makedirs("../plots", exist_ok=True)


def zigzag_nanoribbon_spectrum(N, k_array):
    """Tight-binding spectrum for zigzag nanoribbon with N zigzag chains."""
    t = T_HOP
    a = A_CC
    dim = 2 * N
    energies = np.zeros((len(k_array), dim))

    for ik, k in enumerate(k_array):
        H = np.zeros((dim, dim), dtype=complex)

        for n in range(N):
            ai = 2 * n      # A site
            bi = 2 * n + 1  # B site

            # A-B coupling within the same unit cell
            phase = 1 + np.exp(1j * k)
            H[ai, bi] = -t * phase
            H[bi, ai] = -t * np.conj(phase)

            # B_n to A_{n+1} coupling
            if n < N - 1:
                H[bi, 2*(n+1)] = -t
                H[2*(n+1), bi] = -t

        evals = np.linalg.eigvalsh(H)
        energies[ik, :] = np.sort(evals)

    return energies


def armchair_nanoribbon_spectrum(N, k_array):
    """Tight-binding spectrum for armchair nanoribbon with N dimer lines."""
    t = T_HOP
    dim = 2 * N
    energies = np.zeros((len(k_array), dim))

    for ik, k in enumerate(k_array):
        H = np.zeros((dim, dim), dtype=complex)

        for n in range(N):
            ai = 2 * n
            bi = 2 * n + 1

            # Intra-dimer coupling
            H[ai, bi] = -t
            H[bi, ai] = -t

            # Inter-dimer coupling (along ribbon direction)
            if n < N - 1:
                # A_n to B_{n+1}
                H[ai, 2*(n+1)+1] = -t * np.exp(1j * k / 2)
                H[2*(n+1)+1, ai] = -t * np.exp(-1j * k / 2)
                # B_n to A_{n+1}
                H[bi, 2*(n+1)] = -t * np.exp(-1j * k / 2)
                H[2*(n+1), bi] = -t * np.exp(1j * k / 2)

        evals = np.linalg.eigvalsh(H)
        energies[ik, :] = np.sort(evals)

    return energies


def compute_spectra():
    N = 200  # unit cells wide (as in paper)
    n_k = 300

    # Zigzag: k in [-pi, pi]
    k_zz = np.linspace(-np.pi, np.pi, n_k)
    print("Computing zigzag spectrum (N=200)...")
    E_zz = zigzag_nanoribbon_spectrum(N, k_zz)

    # Armchair: k in [-pi, pi]
    k_ac = np.linspace(-np.pi, np.pi, n_k)
    print("Computing armchair spectrum (N=200)...")
    E_ac = armchair_nanoribbon_spectrum(N, k_ac)

    return k_zz, E_zz, k_ac, E_ac


def save_data(k_zz, E_zz, k_ac, E_ac):
    # Save only a subset of bands near zero energy (first 14 eigenvalues around Fermi level)
    N = E_zz.shape[1]
    mid = N // 2
    n_bands = 7  # 7 above, 7 below

    # Zigzag
    zz_bands = E_zz[:, mid-n_bands:mid+n_bands]
    data_zz = np.column_stack([k_zz, zz_bands])
    cols_zz = "k," + ",".join([f"E_zz_{i}" for i in range(2*n_bands)])
    np.savetxt("../data/fig17_zigzag.csv", data_zz, delimiter=",",
               header=cols_zz, fmt="%.8e", comments="")

    # Armchair
    ac_bands = E_ac[:, mid-n_bands:mid+n_bands]
    data_ac = np.column_stack([k_ac, ac_bands])
    cols_ac = "k," + ",".join([f"E_ac_{i}" for i in range(2*n_bands)])
    np.savetxt("../data/fig17_armchair.csv", data_ac, delimiter=",",
               header=cols_ac, fmt="%.8e", comments="")


def plot_figure(k_zz, E_zz, k_ac, E_ac):
    N = E_zz.shape[1]
    mid = N // 2
    n_bands = 7

    fig, axes = plt.subplots(2, 2, figsize=(12, 8),
                              gridspec_kw={'width_ratios': [2, 1]})

    # Top left: armchair full
    for i in range(mid-n_bands, mid+n_bands):
        axes[0, 0].plot(k_ac, E_ac[:, i], 'b-', linewidth=0.5)
    axes[0, 0].set_ylabel('E (eV)')
    axes[0, 0].set_title('Armchair nanoribbon (N=200)')
    axes[0, 0].set_xlim(-np.pi, np.pi)
    axes[0, 0].set_ylim(-1.5, 1.5)

    # Top right: armchair zoom
    for i in range(mid-n_bands, mid+n_bands):
        axes[0, 1].plot(k_ac, E_ac[:, i], 'b-', linewidth=0.5)
    axes[0, 1].set_title('Zoom')
    axes[0, 1].set_xlim(-np.pi, np.pi)
    axes[0, 1].set_ylim(-0.3, 0.3)

    # Bottom left: zigzag full
    for i in range(mid-n_bands, mid+n_bands):
        axes[1, 0].plot(k_zz, E_zz[:, i], 'b-', linewidth=0.5)
    axes[1, 0].set_xlabel(r'$k$')
    axes[1, 0].set_ylabel('E (eV)')
    axes[1, 0].set_title('Zigzag nanoribbon (N=200)')
    axes[1, 0].set_xlim(-np.pi, np.pi)
    axes[1, 0].set_ylim(-1.5, 1.5)

    # Bottom right: zigzag zoom
    for i in range(mid-n_bands, mid+n_bands):
        axes[1, 1].plot(k_zz, E_zz[:, i], 'b-', linewidth=0.5)
    axes[1, 1].set_xlabel(r'$k$')
    axes[1, 1].set_title('Zoom (edge states)')
    axes[1, 1].set_xlim(-np.pi, np.pi)
    axes[1, 1].set_ylim(-0.3, 0.3)

    plt.tight_layout()
    plt.savefig("../plots/fig17_nanoribbons.png", dpi=150, bbox_inches='tight')
    plt.close()


if True:
    k_zz, E_zz, k_ac, E_ac = compute_spectra()
    save_data(k_zz, E_zz, k_ac, E_ac)
    plot_figure(k_zz, E_zz, k_ac, E_ac)
    print("Figure 17 done: nanoribbon spectra saved.")
