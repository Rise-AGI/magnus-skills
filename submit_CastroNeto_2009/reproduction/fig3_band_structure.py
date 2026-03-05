"""
Figure 3: Energy band structure of graphene.

Left: Full energy spectrum E(k) for t=2.7 eV, t'=0.2t along high-symmetry path.
Right: Zoom near Dirac point.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(os.path.realpath(sys.argv[0]))))
from graphene_params import (A_CC, T_HOP, T_PRIME, energy_dispersion_full,

# Ensure output directories exist
os.makedirs("../data", exist_ok=True)
os.makedirs("../plots", exist_ok=True)
                              lattice_vectors, dirac_points)


def compute_band_structure():
    a = A_CC
    t = T_HOP
    tp = T_PRIME

    # High symmetry points: Gamma -> M -> K -> Gamma
    Gamma = np.array([0.0, 0.0])
    K_pt = np.array([2 * np.pi / (3 * a), 2 * np.pi / (3 * np.sqrt(3) * a)])
    M_pt = np.array([2 * np.pi / (3 * a), 0.0])

    n_points = 300

    # Build path segments
    segments = [
        (Gamma, M_pt, "Gamma-M"),
        (M_pt, K_pt, "M-K"),
        (K_pt, Gamma, "K-Gamma"),
    ]

    all_kx = []
    all_ky = []
    all_dist = []
    tick_positions = [0.0]
    tick_labels = [r'$\Gamma$']

    cumulative_dist = 0.0

    for start, end, label in segments:
        kx_seg = np.linspace(start[0], end[0], n_points, endpoint=False)
        ky_seg = np.linspace(start[1], end[1], n_points, endpoint=False)

        seg_len = np.linalg.norm(end - start)
        dist_seg = np.linspace(0, seg_len, n_points, endpoint=False) + cumulative_dist

        all_kx.extend(kx_seg)
        all_ky.extend(ky_seg)
        all_dist.extend(dist_seg)

        cumulative_dist += seg_len
        tick_positions.append(cumulative_dist)

    # Add final point (Gamma)
    all_kx.append(Gamma[0])
    all_ky.append(Gamma[1])
    all_dist.append(cumulative_dist)

    tick_labels.extend(['M', 'K', r'$\Gamma$'])

    kx_arr = np.array(all_kx)
    ky_arr = np.array(all_ky)
    dist_arr = np.array(all_dist)

    Ep, Em = energy_dispersion_full(kx_arr, ky_arr, t, tp, a)

    return dist_arr, Ep, Em, tick_positions, tick_labels


def save_data(dist, Ep, Em):
    data = np.column_stack([dist, Ep, Em])
    header = "k_path_dist,E_plus,E_minus"
    np.savetxt("../data/fig3_band_structure.csv", data, delimiter=",",
               header=header, fmt="%.8e", comments="")


def plot_figure(dist, Ep, Em, tick_pos, tick_labels):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 6))

    # Left: full spectrum
    ax1.plot(dist, Ep, 'b-', linewidth=0.8)
    ax1.plot(dist, Em, 'b-', linewidth=0.8)
    ax1.set_ylabel('E/t')
    ax1.set_xlim(dist[0], dist[-1])
    ax1.set_ylim(-3.5 * T_HOP, 3.5 * T_HOP)
    ax1.set_xticks(tick_pos)
    ax1.set_xticklabels(tick_labels)
    for tp in tick_pos:
        ax1.axvline(tp, color='gray', linewidth=0.5, linestyle='--')
    ax1.set_title("Energy spectrum (t=2.7 eV, t'=0.2t)")

    # Right: zoom near K point (Dirac point)
    K_idx = 2  # K is the 3rd tick
    K_pos = tick_pos[K_idx]
    zoom_width = (tick_pos[-1] - tick_pos[0]) * 0.08
    mask = (dist > K_pos - zoom_width) & (dist < K_pos + zoom_width)

    ax2.plot(dist[mask], Ep[mask], 'b-', linewidth=1.2)
    ax2.plot(dist[mask], Em[mask], 'b-', linewidth=1.2)
    ax2.set_ylabel('E (eV)')
    ax2.set_xlabel('k')
    ax2.set_title("Zoom near Dirac point")
    ax2.axhline(0, color='gray', linewidth=0.5, linestyle='--')

    plt.tight_layout()
    plt.savefig("../plots/fig3_band_structure.png", dpi=150, bbox_inches='tight')
    plt.close()


if True:
    dist, Ep, Em, tick_pos, tick_labels = compute_band_structure()
    save_data(dist, Ep, Em)
    plot_figure(dist, Ep, Em, tick_pos, tick_labels)
    print("Figure 3 done: band structure saved.")
