"""
Figure 7: Self-consistency check — Fourth-order magnetization cumulant
U_4 vs lattice size L at three fixed K values near K_c.

Runs Wolff MC simulations at K = 0.221654622, 0.221654628, 0.221654634
for multiple lattice sizes and plots U_4(L).

Reproduces Fig. 7 of Ferrenberg, Xu & Landau (2018).
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(".")))
from ising3d import Ising3D, compute_u4


def run_u4_at_K(L, K, n_therm, n_measure, rng):
    """Run Wolff MC at coupling K and measure U_4."""
    model = Ising3D(L, K, rng=rng)
    model.initialize_hot()
    N = L ** 3

    for _ in range(n_therm):
        model.wolff_step()

    m2_samples = np.zeros(n_measure)
    m4_samples = np.zeros(n_measure)
    for i in range(n_measure):
        model.wolff_step()
        m = abs(model.magnetization()) / N
        m2_samples[i] = m ** 2
        m4_samples[i] = m ** 4

    return compute_u4(m2_samples, m4_samples)


def main():
    K_values = [0.221654622, 0.221654628, 0.221654634]
    K_labels = [r'$K = 0.221\,654\,622$', r'$K = 0.221\,654\,628$', r'$K = 0.221\,654\,634$']
    lattice_sizes = [8, 12, 16, 24, 32, 48]
    n_therm = 5000
    n_measure = 40000
    seed = 123

    print("Running self-consistency check (U_4 vs L)...")

    results = {}
    for K in K_values:
        results[K] = []
        for L in lattice_sizes:
            print(f"  K = {K:.9f}, L = {L}...", end=" ", flush=True)
            rng = np.random.default_rng(seed + L + int(K * 1e12) % 1000)
            u4 = run_u4_at_K(L, K, n_therm, n_measure, rng)
            results[K].append(u4)
            print(f"U_4 = {u4:.6f}")

    # --- CSV output ---
    with open("../data/fig7_u4_vs_L.csv", "w") as f:
        f.write("# Figure 7: Self-consistency check, U_4 vs L at fixed K values\n")
        f.write("# Columns: L, U4_K1 (K=0.221654622), U4_K2 (K=0.221654628), U4_K3 (K=0.221654634)\n")
        f.write("L,U4_K_0.221654622,U4_K_0.221654628,U4_K_0.221654634\n")
        for i, L in enumerate(lattice_sizes):
            f.write(f"{L}")
            for K in K_values:
                f.write(f",{results[K][i]:.8f}")
            f.write("\n")

    # --- Plot ---
    fig, ax = plt.subplots(figsize=(8, 6))
    colors = ['blue', 'red', 'green']
    markers = ['o', 's', '^']

    for idx, K in enumerate(K_values):
        ax.plot(lattice_sizes, results[K], f'{markers[idx]}-', color=colors[idx],
                label=K_labels[idx], markersize=6)

    ax.axhline(y=0.46548, color='black', linestyle='--', linewidth=0.8,
               label=r'$U^* = 0.46548$')

    ax.set_xlabel(r'$L$', fontsize=14)
    ax.set_ylabel(r'$U_4$', fontsize=14)
    ax.set_title(r'$U_4$ vs $L$ for fixed $K$ values', fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig("../plots/fig7_u4_vs_L.png", dpi=300)
    print("Saved fig7_u4_vs_L.png and fig7_u4_vs_L.csv")


if __name__ == "__main__":
    main()
