"""
Figure 5: Fourth-order magnetization cumulant U_4 vs inverse
temperature K for different lattice sizes.

This figure demonstrates the Binder cumulant crossing technique.
We run Wolff cluster MC simulations for small lattice sizes and
use histogram reweighting to obtain U_4(K) curves.

Reproduces Fig. 5 of Ferrenberg, Xu & Landau (2018).
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(".")))
from ising3d import Ising3D, compute_u4, KC_BEST


def run_simulation(L, K0, n_therm, n_measure, rng):
    """Run Wolff MC and collect energy/magnetization samples."""
    model = Ising3D(L, K0, rng=rng)
    model.initialize_hot()
    N = L ** 3

    # Thermalization
    for _ in range(n_therm):
        model.wolff_step()

    # Measurement
    energies = np.zeros(n_measure)
    mag_abs = np.zeros(n_measure)
    for i in range(n_measure):
        model.wolff_step()
        E = model.energy_fast()
        M = model.magnetization()
        energies[i] = E
        mag_abs[i] = abs(M) / N

    return energies, mag_abs


def compute_u4_reweighted(energies, mag_abs, K0, K_target):
    """Compute U_4 at K_target by reweighting from K0."""
    dK = K0 - K_target
    log_weights = dK * energies
    log_weights -= np.max(log_weights)  # numerical stability
    weights = np.exp(log_weights)
    W = np.sum(weights)

    m2 = mag_abs ** 2
    m4 = mag_abs ** 4

    m2_avg = np.sum(weights * m2) / W
    m4_avg = np.sum(weights * m4) / W

    if m2_avg == 0:
        return 0.0
    return 1.0 - m4_avg / (3.0 * m2_avg ** 2)


def main():
    K0 = 0.221654  # simulation point
    lattice_sizes = [8, 12, 16, 24, 32]
    n_therm = 5000
    n_measure = 50000
    seed = 42

    # K range for reweighting
    K_range = np.linspace(0.2210, 0.2225, 200)

    print("Running Wolff MC simulations for cumulant crossing...")

    results = {}
    for L in lattice_sizes:
        print(f"  L = {L}...", end=" ", flush=True)
        rng = np.random.default_rng(seed + L)
        energies, mag_abs = run_simulation(L, K0, n_therm, n_measure, rng)

        u4_curve = np.array([compute_u4_reweighted(energies, mag_abs, K0, K)
                             for K in K_range])
        results[L] = u4_curve
        print(f"done (U_4 at K0 = {u4_curve[np.argmin(abs(K_range - K0))]:.4f})")

    # --- CSV output ---
    with open("../data/fig5_cumulant_crossing.csv", "w") as f:
        f.write("# Figure 5: Fourth-order magnetization cumulant U_4 vs K\n")
        f.write("# Columns: K")
        for L in lattice_sizes:
            f.write(f",U4_L{L}")
        f.write("\n")
        f.write("K")
        for L in lattice_sizes:
            f.write(f",U4_L{L}")
        f.write("\n")
        for i, K in enumerate(K_range):
            f.write(f"{K:.10f}")
            for L in lattice_sizes:
                f.write(f",{results[L][i]:.8f}")
            f.write("\n")

    # --- Plot ---
    fig, ax = plt.subplots(figsize=(8, 6))
    colors = ['blue', 'red', 'green', 'purple', 'orange']

    for idx, L in enumerate(lattice_sizes):
        ax.plot(K_range, results[L], color=colors[idx], label=f'L = {L}')

    ax.axvline(x=KC_BEST, color='black', linestyle='--', linewidth=0.8,
               label=rf'$K_c = {KC_BEST}$')

    ax.set_xlabel(r'$K = J/k_BT$', fontsize=14)
    ax.set_ylabel(r'$U_4$', fontsize=14)
    ax.set_title(r'Fourth-order magnetization cumulant $U_4$ vs $K$', fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig("../plots/fig5_cumulant_crossing.png", dpi=300)
    print("Saved fig5_cumulant_crossing.png and fig5_cumulant_crossing.csv")


if __name__ == "__main__":
    main()
