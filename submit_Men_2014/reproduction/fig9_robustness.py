"""
Reproduce Figure 9: Fabrication-Adaptive (FA) robustness analysis.

Shows how the band gap changes when the structure is dilated/eroded
by a parameter delta (modeling fabrication errors).

For Diamond2 with r=0.1a, we vary r by +/- delta and compute the gap.
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

    n_max = 3
    grid_res = 24
    n_bands = 6

    kpts = fcc_ibz_kpoints(b1, b2, b3, n_per_edge=6)

    # Delta values (fractional change in radius)
    r_nominal = 0.1
    deltas = np.linspace(-0.03, 0.03, 13)

    results = []
    for delta in deltas:
        r_eff = r_nominal + delta
        if r_eff <= 0.01:
            results.append((delta, 0.0))
            continue

        print(f"  delta={delta:.4f}, r_eff={r_eff:.4f}...")
        eps_func = make_diamond2_eps(r=r_eff, eps_hi=12.96)
        solver = PWESolver3D(a1, a2, a3, eps_func, n_max=n_max, grid_res=grid_res)

        all_freqs = solver.compute_bands(kpts, n_bands=n_bands)
        gap = compute_gap(all_freqs, 1, 2)

        results.append((delta, gap))
        print(f"    gap = {gap:.2f}%")

    # Save data
    data_dir = os.path.join(os.path.dirname(__file__), '..', 'data')
    os.makedirs(data_dir, exist_ok=True)

    with open(os.path.join(data_dir, 'fig9_robustness.csv'), 'w') as f:
        f.write("delta,gap_percent\n")
        for delta, gap in results:
            f.write(f"{delta:.8f},{gap:.8f}\n")

    # Plot
    plot_dir = os.path.join(os.path.dirname(__file__), '..', 'plots')
    os.makedirs(plot_dir, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8, 6))
    ds = [r[0] for r in results]
    gaps = [r[1] for r in results]

    ax.plot(ds, gaps, 'ro-', linewidth=1.5, markersize=6)
    ax.axhline(0, color='gray', linewidth=0.5)
    ax.axvline(0, color='gray', linewidth=0.5, linestyle='--')
    ax.set_xlabel(r'Fabrication error $\delta/a$')
    ax.set_ylabel('Fractional band gap (%)')
    ax.set_title(f'Diamond2 Gap Robustness to Fabrication Error (n_max={n_max})')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'fig9_robustness.png'), dpi=150)
    print(f"  Plot saved to plots/fig9_robustness.png")


if __name__ == '__main__':
    main()
