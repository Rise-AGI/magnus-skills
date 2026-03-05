"""
Fig 34: Spin correlation function of the Heisenberg chain.

Reproduces Figure 34 from Sandvik (2010).
Shows |C(r)| * r vs r for the ground state, where C(r) = <S_0 . S_r>.

The correlation function decays as (-1)^r / r with logarithmic corrections.
"""
import sys
import os
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) if os.path.dirname(os.path.abspath(__file__)) else '.')
from heisenberg_ed import compute_all_correlations

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    chain_lengths = [16, 20, 24]
    J = 1.0

    all_corr = {}
    for N in chain_lengths:
        print(f"Computing correlations for N={N}...")
        e0, corr = compute_all_correlations(N, J=J)
        all_corr[N] = corr
        print(f"  E0/N = {e0/N:.8f}")
        for r in range(len(corr)):
            print(f"  C({r}) = {corr[r]:.8f}")

    # Save data
    os.makedirs('../data', exist_ok=True)
    max_r = max(len(c) for c in all_corr.values())
    with open('../data/fig34_correlations.csv', 'w') as f:
        f.write('r,' + ','.join([f'C_r_N{N}' for N in chain_lengths]) + '\n')
        for r in range(max_r):
            vals = [str(r)]
            for N in chain_lengths:
                if r < len(all_corr[N]):
                    vals.append(f'{all_corr[N][r]:.8f}')
                else:
                    vals.append('')
            f.write(','.join(vals) + '\n')

    # Plot: |C(r)| * r vs r
    os.makedirs('../plots', exist_ok=True)
    fig, ax = plt.subplots(figsize=(8, 6))
    markers = ['o', 's', '^']
    for i, N in enumerate(chain_lengths):
        corr = all_corr[N]
        r_vals = np.arange(1, len(corr))
        c_vals = np.abs(corr[1:]) * r_vals
        ax.plot(r_vals, c_vals, markers[i] + '-', label=f'N={N}', markersize=5)

    ax.set_xlabel('$r$', fontsize=14)
    ax.set_ylabel('$|C(r)| \\times r$', fontsize=14)
    ax.set_title('Spin Correlations in the Heisenberg Chain', fontsize=14)
    ax.legend(fontsize=12)
    plt.tight_layout()
    plt.savefig('../plots/fig34_correlations.png', dpi=150)
    print("Saved fig34_correlations.png")


if __name__ == '__main__':
    main()
