"""
Reproduce Figure 6: Diamond2 photonic crystal band structure.

Diamond2 = FCC lattice with 2-atom basis, tetrahedral cylinder bonds.
Parameters: r=0.1a, eps=12.96
Expected gap: ~31.56% between bands 2-3.
"""
import sys
import os
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))

from pwe3d import PWESolver3D, reciprocal_vectors, make_k_path, compute_gap
from crystals import fcc_lattice, make_diamond2_eps, fcc_kpath

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    a1, a2, a3 = fcc_lattice()
    b1, b2, b3 = reciprocal_vectors(a1, a2, a3)

    eps_func = make_diamond2_eps(r=0.1, eps_hi=12.96)

    n_max = 4
    grid_res = 32
    n_bands = 8

    print(f"Diamond2: n_max={n_max}, grid_res={grid_res}")
    print(f"  N_G = {(2*n_max+1)**3}")

    solver = PWESolver3D(a1, a2, a3, eps_func, n_max=n_max, grid_res=grid_res)

    path_pts, labels = fcc_kpath(b1, b2, b3)
    k_list, k_dist, tick_pos = make_k_path(path_pts, n_per_segment=20)

    print(f"  Computing {len(k_list)} k-points...")
    all_freqs = solver.compute_bands(k_list, n_bands=n_bands, verbose=True)

    gap = compute_gap(all_freqs, 1, 2)  # bands 2-3 (0-indexed: 1,2)
    print(f"  Gap between bands 2-3: {gap:.2f}%")

    # Save data
    data_dir = os.path.join(os.path.dirname(__file__), '..', 'data')
    os.makedirs(data_dir, exist_ok=True)

    header = "k_dist," + ",".join([f"band_{i+1}" for i in range(n_bands)])
    data = np.column_stack([k_dist, all_freqs])
    np.savetxt(os.path.join(data_dir, 'fig6_diamond2_bands.csv'),
               data, delimiter=',', header=header, comments='',
               fmt='%.8f')

    with open(os.path.join(data_dir, 'fig6_diamond2_gap.csv'), 'w') as f:
        f.write("structure,gap_bands,gap_percent,n_max,grid_res\n")
        f.write(f"Diamond2,2-3,{gap:.8f},{n_max},{grid_res}\n")

    # Plot
    plot_dir = os.path.join(os.path.dirname(__file__), '..', 'plots')
    os.makedirs(plot_dir, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8, 6))
    for i in range(n_bands):
        ax.plot(k_dist, all_freqs[:, i], 'b-', linewidth=0.8)

    max_below = np.max(all_freqs[:, 1])
    min_above = np.min(all_freqs[:, 2])
    if min_above > max_below:
        ax.axhspan(max_below, min_above, alpha=0.3, color='yellow',
                    label=f'Gap: {gap:.1f}%')
        ax.legend()

    for tp in tick_pos:
        ax.axvline(tp, color='gray', linewidth=0.5, linestyle='--')

    ax.set_xticks(tick_pos)
    ax.set_xticklabels(labels)
    ax.set_ylabel(r'Frequency $\omega a / 2\pi c$')
    ax.set_title(f'Diamond2 Band Structure (n_max={n_max}, gap={gap:.1f}%)')
    ax.set_xlim(k_dist[0], k_dist[-1])
    ax.set_ylim(0, None)

    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'fig6_diamond2_bands.png'), dpi=150)
    print(f"  Plot saved to plots/fig6_diamond2_bands.png")


if __name__ == '__main__':
    main()
