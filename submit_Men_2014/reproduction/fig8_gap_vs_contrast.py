"""
Reproduce Figure 8: Gap size vs. dielectric contrast for Diamond2.

Computes the fractional band gap as a function of eps_hi/eps_lo
for the Diamond2 parametric structure with r=0.1a.

The paper shows a threshold around n~1.9:1 (eps~3.6:1) and saturation ~30%.
"""
import sys
import os
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))

from pwe3d import PWESolver3D, reciprocal_vectors, compute_gap
from crystals import fcc_lattice, make_diamond2_eps, fcc_ibz_kpoints

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    a1, a2, a3 = fcc_lattice()
    b1, b2, b3 = reciprocal_vectors(a1, a2, a3)

    # Sample IBZ k-points for gap computation
    n_max = 3
    grid_res = 24
    n_bands = 6

    # Dielectric contrast values (eps_hi with eps_lo=1)
    eps_values = [1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.96, 16.0]

    kpts = fcc_ibz_kpoints(b1, b2, b3, n_per_edge=6)

    results = []
    for eps_hi in eps_values:
        contrast = np.sqrt(eps_hi)  # n_hi/n_lo
        print(f"  eps_hi={eps_hi:.2f} (n={contrast:.2f})...")

        eps_func = make_diamond2_eps(r=0.1, eps_hi=eps_hi)
        solver = PWESolver3D(a1, a2, a3, eps_func, n_max=n_max, grid_res=grid_res)

        all_freqs = solver.compute_bands(kpts, n_bands=n_bands)
        gap = compute_gap(all_freqs, 1, 2)

        results.append((eps_hi, contrast, gap))
        print(f"    gap = {gap:.2f}%")

    # Save data
    data_dir = os.path.join(os.path.dirname(__file__), '..', 'data')
    os.makedirs(data_dir, exist_ok=True)

    with open(os.path.join(data_dir, 'fig8_gap_vs_contrast.csv'), 'w') as f:
        f.write("eps_hi,n_contrast,gap_percent\n")
        for eps_hi, contrast, gap in results:
            f.write(f"{eps_hi:.8f},{contrast:.8f},{gap:.8f}\n")

    # Plot
    plot_dir = os.path.join(os.path.dirname(__file__), '..', 'plots')
    os.makedirs(plot_dir, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8, 6))
    contrasts = [r[1] for r in results]
    gaps = [r[2] for r in results]

    ax.plot(contrasts, gaps, 'bo-', linewidth=1.5, markersize=6)
    ax.axhline(0, color='gray', linewidth=0.5)
    ax.set_xlabel(r'Index contrast $n_{hi}/n_{lo}$')
    ax.set_ylabel('Fractional band gap (%)')
    ax.set_title(f'Diamond2 Gap vs. Index Contrast (n_max={n_max})')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'fig8_gap_vs_contrast.png'), dpi=150)
    print(f"  Plot saved to plots/fig8_gap_vs_contrast.png")


if __name__ == '__main__':
    main()
