"""
Figure 2: Critical exponent nu vs L_min with different numbers of
correction exponents.

Reproduces Table II / Fig. 2 of Ferrenberg, Xu & Landau (2018).
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(".")))
from ising3d import TABLE_II, NU_BEST

def main():
    data = TABLE_II
    L_min = np.array(data["L_min"])

    # 1 correction exponent
    nu1 = np.array(data["nu_1corr"])
    nu1_err = np.array(data["nu_1corr_err"])

    # 2 correction exponents
    nu2 = np.array(data["nu_2corr"])
    nu2_err = np.array(data["nu_2corr_err"])

    # 3 correction exponents (last entry is None)
    mask3 = np.array([x is not None for x in data["nu_3corr"]])
    L_min3 = L_min[mask3]
    nu3 = np.array([x for x in data["nu_3corr"] if x is not None])
    nu3_err = np.array([x for x in data["nu_3corr_err"] if x is not None])

    # --- CSV output ---
    with open("../data/fig2_nu_vs_lmin.csv", "w") as f:
        f.write("# Figure 2: Critical exponent nu vs L_min\n")
        f.write("# Columns: L_min, nu_1corr, nu_1corr_err, nu_2corr, nu_2corr_err, nu_3corr, nu_3corr_err\n")
        f.write("L_min,nu_1corr,nu_1corr_err,nu_2corr,nu_2corr_err,nu_3corr,nu_3corr_err\n")
        for i in range(len(L_min)):
            n3 = data["nu_3corr"][i]
            n3e = data["nu_3corr_err"][i]
            n3_str = f"{n3:.6f}" if n3 is not None else ""
            n3e_str = f"{n3e:.6f}" if n3e is not None else ""
            f.write(f"{L_min[i]},{nu1[i]:.6f},{nu1_err[i]:.6f},"
                    f"{nu2[i]:.6f},{nu2_err[i]:.6f},{n3_str},{n3e_str}\n")

    # --- Plot ---
    fig, ax = plt.subplots(figsize=(8, 6))

    ax.errorbar(L_min, nu1, yerr=nu1_err, fmt='o-', color='blue',
                label=r'$\omega_1 = 0.83$ fixed', capsize=3, markersize=5)
    ax.errorbar(L_min, nu2, yerr=nu2_err, fmt='s-', color='red',
                label=r'$\omega_1 = 0.83,\; \omega_2 = 4$ fixed', capsize=3, markersize=5)
    ax.errorbar(L_min3, nu3, yerr=nu3_err, fmt='^-', color='green',
                label=r'$\omega_1 = 0.83,\; \omega_2 = 4,\; \omega_\nu = 1.6$ fixed', capsize=3, markersize=5)

    ax.axhline(y=NU_BEST, color='black', linestyle='--', linewidth=0.8,
               label=rf'$\nu = {NU_BEST}$')

    ax.set_xlabel(r'$L_{\min}$', fontsize=14)
    ax.set_ylabel(r'$\nu$', fontsize=14)
    ax.set_title(r'Critical exponent $\nu$ vs $L_{\min}$', fontsize=14)
    ax.legend(fontsize=10)
    ax.set_xlim(10, 170)
    ax.set_ylim(0.6295, 0.6325)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig("../plots/fig2_nu_vs_lmin.png", dpi=300)
    print("Saved fig2_nu_vs_lmin.png and fig2_nu_vs_lmin.csv")


if __name__ == "__main__":
    main()
