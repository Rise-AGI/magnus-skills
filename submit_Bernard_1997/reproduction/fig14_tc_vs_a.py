"""
Figure 14: Tc/sqrt(sigma) vs a*sqrt(sigma)

Scaling plot showing critical temperature scaled by string tension
versus lattice spacing in units of 1/sqrt(sigma).
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(".")))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from paper_data import (TABLE2, TABLE3,
                        compute_tc_over_sqrtsigma, compute_a_sqrtsigma,
                        QUENCHED_TC_SQRTSIGMA_NT4)


def main():
    # Improved Wilson Nt=4
    iw_a_sqrts = []
    iw_a_sqrts_err = []
    iw_tc_sqrts = []
    iw_tc_sqrts_err = []

    for row in TABLE2:
        a2s, a2s_err = row[7], row[8]
        tc_s, tc_s_err = compute_tc_over_sqrtsigma(a2s, a2s_err, Nt=4)
        a_s, a_s_err = compute_a_sqrtsigma(a2s, a2s_err)

        iw_a_sqrts.append(a_s)
        iw_a_sqrts_err.append(a_s_err)
        iw_tc_sqrts.append(tc_s)
        iw_tc_sqrts_err.append(tc_s_err)

    # KS data
    ks_a_sqrts = []
    ks_a_sqrts_err = []
    ks_tc_sqrts = []
    ks_tc_sqrts_err = []
    ks_Nt = []

    for row in TABLE3:
        a2s, a2s_err = row[7], row[8]
        Nt = row[15]
        tc_s, tc_s_err = compute_tc_over_sqrtsigma(a2s, a2s_err, Nt=Nt)
        a_s, a_s_err = compute_a_sqrtsigma(a2s, a2s_err)

        ks_a_sqrts.append(a_s)
        ks_a_sqrts_err.append(a_s_err)
        ks_tc_sqrts.append(tc_s)
        ks_tc_sqrts_err.append(tc_s_err)
        ks_Nt.append(Nt)

    # Save data
    header = "a_sqrtsigma,a_sqrtsigma_err,Tc_sqrtsigma,Tc_sqrtsigma_err,action"
    lines = [header]
    for i in range(len(iw_a_sqrts)):
        lines.append(f"{iw_a_sqrts[i]:.8f},{iw_a_sqrts_err[i]:.8f},"
                     f"{iw_tc_sqrts[i]:.8f},{iw_tc_sqrts_err[i]:.8f},"
                     f"improved_wilson_Nt4")
    for i in range(len(ks_a_sqrts)):
        lines.append(f"{ks_a_sqrts[i]:.8f},{ks_a_sqrts_err[i]:.8f},"
                     f"{ks_tc_sqrts[i]:.8f},{ks_tc_sqrts_err[i]:.8f},"
                     f"KS_Nt{ks_Nt[i]}")
    with open("../data/fig14_tc_vs_a.csv", "w") as f:
        f.write("\n".join(lines) + "\n")

    # Plot
    fig, ax = plt.subplots(figsize=(8, 6))

    ax.errorbar(iw_a_sqrts, iw_tc_sqrts, xerr=iw_a_sqrts_err, yerr=iw_tc_sqrts_err,
                fmt='o', color='blue', markersize=10, markerfacecolor='white',
                markeredgewidth=2, capsize=3,
                label=r'$N_t=4$ improved Wilson', zorder=5)

    ks4 = [i for i in range(len(ks_Nt)) if ks_Nt[i] == 4]
    ks6 = [i for i in range(len(ks_Nt)) if ks_Nt[i] == 6]

    if ks4:
        ax.errorbar([ks_a_sqrts[i] for i in ks4], [ks_tc_sqrts[i] for i in ks4],
                    xerr=[ks_a_sqrts_err[i] for i in ks4],
                    yerr=[ks_tc_sqrts_err[i] for i in ks4],
                    fmt='D', color='red', markersize=8, capsize=3,
                    label=r'$N_t=4$ KS', zorder=4)
    if ks6:
        ax.errorbar([ks_a_sqrts[i] for i in ks6], [ks_tc_sqrts[i] for i in ks6],
                    xerr=[ks_a_sqrts_err[i] for i in ks6],
                    yerr=[ks_tc_sqrts_err[i] for i in ks6],
                    fmt='s', color='green', markersize=8, capsize=3,
                    label=r'$N_t=6$ KS', zorder=4)

    ax.set_xlabel(r'$a\sqrt{\sigma}$', fontsize=14)
    ax.set_ylabel(r'$T_c/\sqrt{\sigma}$', fontsize=14)
    ax.set_title(r'$T_c/\sqrt{\sigma}$ vs $a\sqrt{\sigma}$', fontsize=14)
    ax.legend(fontsize=11)
    ax.set_xlim(0.0, 0.8)
    ax.set_ylim(0.3, 0.6)

    plt.tight_layout()
    plt.savefig("../plots/fig14_tc_vs_a.png", dpi=150)
    plt.close()
    print("Fig 14: Tc vs a scaling saved.")


main()
