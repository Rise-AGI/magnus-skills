"""
Figure 4: Critical coupling K_c vs L_min with 1, 2, 3 fixed correction
exponents.

Reproduces Table V / Fig. 4 of Ferrenberg, Xu & Landau (2018).
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(".")))
from ising3d import TABLE_V, KC_BEST

def main():
    data = TABLE_V
    L_min = np.array(data["L_min"])

    Kc1 = np.array(data["Kc_1corr"])
    Kc1_err = np.array(data["Kc_1corr_err"])

    Kc2 = np.array(data["Kc_2corr"])
    Kc2_err = np.array(data["Kc_2corr_err"])

    mask3 = np.array([x is not None for x in data["Kc_3corr"]])
    L_min3 = L_min[mask3]
    Kc3 = np.array([x for x in data["Kc_3corr"] if x is not None])
    Kc3_err = np.array([x for x in data["Kc_3corr_err"] if x is not None])

    # --- CSV output ---
    with open("../data/fig4_kc_vs_lmin.csv", "w") as f:
        f.write("# Figure 4: Critical coupling K_c vs L_min with fixed correction exponents\n")
        f.write("# Columns: L_min, Kc_1corr, Kc_1corr_err, Kc_2corr, Kc_2corr_err, Kc_3corr, Kc_3corr_err\n")
        f.write("L_min,Kc_1corr,Kc_1corr_err,Kc_2corr,Kc_2corr_err,Kc_3corr,Kc_3corr_err\n")
        for i in range(len(L_min)):
            k3 = data["Kc_3corr"][i]
            k3e = data["Kc_3corr_err"][i]
            k3_str = f"{k3:.10f}" if k3 is not None else ""
            k3e_str = f"{k3e:.1e}" if k3e is not None else ""
            f.write(f"{L_min[i]},{Kc1[i]:.10f},{Kc1_err[i]:.1e},"
                    f"{Kc2[i]:.10f},{Kc2_err[i]:.1e},{k3_str},{k3e_str}\n")

    # --- Plot ---
    fig, ax = plt.subplots(figsize=(8, 6))

    ax.errorbar(L_min, Kc1 * 1e9 - 221654620, yerr=Kc1_err * 1e9, fmt='o-', color='blue',
                label=r'$\omega_1 = 0.83$ fixed', capsize=3, markersize=5)
    ax.errorbar(L_min, Kc2 * 1e9 - 221654620, yerr=Kc2_err * 1e9, fmt='s-', color='red',
                label=r'$\omega_1 = 0.83,\; \omega_2 = 4$ fixed', capsize=3, markersize=5)
    ax.errorbar(L_min3, Kc3 * 1e9 - 221654620, yerr=Kc3_err * 1e9, fmt='^-', color='green',
                label=r'$\omega_1 = 0.83,\; \omega_2 = 4,\; \omega_\nu = 1.6$ fixed', capsize=3, markersize=5)

    ax.axhline(y=KC_BEST * 1e9 - 221654620, color='black', linestyle='--', linewidth=0.8,
               label=rf'$K_c = {KC_BEST}$')

    ax.set_xlabel(r'$L_{\min}$', fontsize=14)
    ax.set_ylabel(r'$K_c - 0.221654620 \;\; (\times 10^{-9})$', fontsize=14)
    ax.set_title(r'Critical coupling $K_c$ vs $L_{\min}$ (fixed correction exponents)', fontsize=14)
    ax.legend(fontsize=9)
    ax.set_xlim(10, 170)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig("../plots/fig4_kc_vs_lmin.png", dpi=300)
    print("Saved fig4_kc_vs_lmin.png and fig4_kc_vs_lmin.csv")


if __name__ == "__main__":
    main()
