"""
Figure 6: Critical coupling K_c from cumulant crossing technique
vs L_min, with one and two correction terms.

Reproduces Table VII, VIII / Fig. 6 of Ferrenberg, Xu & Landau (2018).
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(".")))
from ising3d import TABLE_VII, TABLE_VIII, KC_BEST

def main():
    # One correction term (Table VII)
    L1 = np.array(TABLE_VII["L_min"])
    Kc1 = np.array(TABLE_VII["Kc"])
    Kc1_err = np.array(TABLE_VII["Kc_err"])

    # Two correction terms (Table VIII)
    L2 = np.array(TABLE_VIII["L_min"])
    Kc2 = np.array(TABLE_VIII["Kc"])
    Kc2_err = np.array(TABLE_VIII["Kc_err"])

    # --- CSV output ---
    with open("../data/fig6_kc_crossing.csv", "w") as f:
        f.write("# Figure 6: K_c from cumulant crossing technique vs L_min\n")
        f.write("# One correction term data from Table VII\n")
        f.write("# Columns: L_min, Kc_1corr, Kc_1corr_err\n")
        f.write("L_min,Kc_1corr,Kc_1corr_err\n")
        for i in range(len(L1)):
            f.write(f"{L1[i]},{Kc1[i]:.11f},{Kc1_err[i]:.2e}\n")

    with open("../data/fig6_kc_crossing_2corr.csv", "w") as f:
        f.write("# Figure 6: K_c from cumulant crossing technique vs L_min (2 correction terms)\n")
        f.write("# Two correction terms data from Table VIII\n")
        f.write("# Columns: L_min, Kc_2corr, Kc_2corr_err\n")
        f.write("L_min,Kc_2corr,Kc_2corr_err\n")
        for i in range(len(L2)):
            f.write(f"{L2[i]},{Kc2[i]:.11f},{Kc2_err[i]:.2e}\n")

    # --- Plot ---
    # Offset for visualization: show (Kc - 0.22165462) * 1e8
    offset = 0.22165462
    scale = 1e8

    fig, ax = plt.subplots(figsize=(8, 6))

    ax.errorbar(L1, (Kc1 - offset) * scale, yerr=Kc1_err * scale,
                fmt='o-', color='blue', label='1 correction term', capsize=3, markersize=5)
    ax.errorbar(L2, (Kc2 - offset) * scale, yerr=Kc2_err * scale,
                fmt='s-', color='red', label='2 correction terms', capsize=3, markersize=5)

    ax.axhline(y=(KC_BEST - offset) * scale, color='black', linestyle='--',
               linewidth=0.8, label=rf'$K_c = {KC_BEST}$')

    ax.set_xlabel(r'$L_{\min}$', fontsize=14)
    ax.set_ylabel(r'$(K_c - 0.22165462) \times 10^{8}$', fontsize=14)
    ax.set_title(r'$K_c$ from cumulant crossing technique', fontsize=14)
    ax.legend(fontsize=10)
    ax.set_xlim(10, 200)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig("../plots/fig6_kc_crossing.png", dpi=300)
    print("Saved fig6_kc_crossing.png and fig6_kc_crossing.csv")


if __name__ == "__main__":
    main()
