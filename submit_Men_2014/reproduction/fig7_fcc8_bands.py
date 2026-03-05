"""
Reproduce Figure 7: FCC8 photonic crystal band structure.

FCC8 = FCC lattice of hollow dielectric spheres connected by rods to 12 neighbors.
Parameters: r1=0.12a, r2=0.19a, r3=0.08a, eps=12.96
Expected gap: ~18.3% between bands 8-9.
"""
import sys
import os
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))

from pwe3d import PWESolver3D, reciprocal_vectors, make_k_path, compute_gap
from crystals import fcc_lattice, make_fcc8_eps, fcc_kpath

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    a1, a2, a3 = fcc_lattice()
    b1, b2, b3 = reciprocal_vectors(a1, a2, a3)

    eps_func = make_fcc8_eps(r1=0.12, r2=0.19, r3=0.08, eps_hi=12.96)

    n_max = 4
    grid_res = 32
    n_bands = 14

    print(f"FCC8: n_max={n_max}, grid_res={grid_res}")
    print(f"  N_G = {(2*n_max+1)**3}")

    solver = PWESolver3D(a1, a2, a3, eps_func, n_max=n_max, grid_res=grid_res)

    path_pts, labels = fcc_kpath(b1, b2, b3)
    k_list, k_dist, tick_pos = make_k_path(path_pts, n_per_segment=20)

    print(f"  Computing {len(k_list)} k-points...")
    all_freqs = solver.compute_bands(k_list, n_bands=n_bands, verbose=True)

    gap = compute_gap(all_freqs, 7, 8)  # bands 8-9 (0-indexed: 7,8)
    print(f"  Gap between bands 8-9: {gap:.2f}%")

    # Save data
    data_dir = os.path.join(os.path.dirname(__file__), '..', 'data')
    os.makedirs(data_dir, exist_ok=True)

    header = "k_dist," + ",".join([f"band_{i+1}" for i in range(n_bands)])
    data = np.column_stack([k_dist, all_freqs])
    np.savetxt(os.path.join(data_dir, 'fig7_fcc8_bands.csv'),
               data, delimiter=',', header=header, comments='',
               fmt='%.8f')

    with open(os.path.join(data_dir, 'fig7_fcc8_gap.csv'), 'w') as f:
        f.write("structure,gap_bands,gap_percent,n_max,grid_res\n")
        f.write(f"FCC8,8-9,{gap:.8f},{n_max},{grid_res}\n")

    # Plot
    plot_dir = os.path.join(os.path.dirname(__file__), '..', 'plots')
    os.makedirs(plot_dir, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8, 6))
    for i in range(n_bands):
        ax.plot(k_dist, all_freqs[:, i], 'b-', linewidth=0.8)

    max_below = np.max(all_freqs[:, 7])
    min_above = np.min(all_freqs[:, 8])
    if min_above > max_below:
        ax.axhspan(max_below, min_above, alpha=0.3, color='yellow',
                    label=f'Gap: {gap:.1f}%')
        ax.legend()

    for tp in tick_pos:
        ax.axvline(tp, color='gray', linewidth=0.5, linestyle='--')

    ax.set_xticks(tick_pos)
    ax.set_xticklabels(labels)
    ax.set_ylabel(r'Frequency $\omega a / 2\pi c$')
    ax.set_title(f'FCC8 Band Structure (n_max={n_max}, gap={gap:.1f}%)')
    ax.set_xlim(k_dist[0], k_dist[-1])
    ax.set_ylim(0, None)

    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'fig7_fcc8_bands.png'), dpi=150)
    print(f"  Plot saved to plots/fig7_fcc8_bands.png")


if __name__ == '__main__':
    main()
