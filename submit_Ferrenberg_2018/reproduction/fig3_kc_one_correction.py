"""
Figure 3: Critical coupling K_c vs L_min with one correction term
(unfixed omega).

Reproduces Table IV / Fig. 3 of Ferrenberg, Xu & Landau (2018).
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(".")))
from ising3d import TABLE_IV, KC_BEST

def main():
    data = TABLE_IV
    L_min = np.array(data["L_min"])
    Kc = np.array(data["Kc"])
    Kc_err = np.array(data["Kc_err"])

    # --- CSV output ---
    with open("../data/fig3_kc_one_correction.csv", "w") as f:
        f.write("# Figure 3: K_c vs L_min with one correction term (unfixed omega)\n")
        f.write("# From Table IV of Ferrenberg, Xu & Landau (2018)\n")
        f.write("# Columns: L_min, Kc, Kc_err\n")
        f.write("L_min,Kc,Kc_err\n")
        for i in range(len(L_min)):
            f.write(f"{L_min[i]},{Kc[i]:.10f},{Kc_err[i]:.2e}\n")

    # --- Plot ---
    offset = 0.221654620
    scale = 1e9

    fig, ax = plt.subplots(figsize=(8, 6))

    ax.errorbar(L_min, (Kc - offset) * scale, yerr=Kc_err * scale,
                fmt='o-', color='blue', capsize=3, markersize=6,
                label=r'$K_c(L_{\min})$, one correction term')

    ax.axhline(y=(KC_BEST - offset) * scale, color='black', linestyle='--',
               linewidth=0.8, label=rf'$K_c = {KC_BEST}$')

    ax.set_xlabel(r'$L_{\min}$', fontsize=14)
    ax.set_ylabel(r'$(K_c - 0.221654620) \times 10^{9}$', fontsize=14)
    ax.set_title(r'$K_c$ vs $L_{\min}$ (one correction term, unfixed $\omega$)', fontsize=14)
    ax.legend(fontsize=11)
    ax.set_xlim(10, 170)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig("../plots/fig3_kc_one_correction.png", dpi=300)
    print("Saved fig3_kc_one_correction.png and fig3_kc_one_correction.csv")


if __name__ == "__main__":
    main()
