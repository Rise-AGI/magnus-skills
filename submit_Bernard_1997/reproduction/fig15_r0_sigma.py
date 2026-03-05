"""
Figure 15: r0*sqrt(sigma) vs a/r0

Scaling plot showing the dimensionless quantity r0*sqrt(sigma) as a
function of lattice spacing in units of r0.
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(".")))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from paper_data import (TABLE2, TABLE3,
                        compute_r0_sqrtsigma, compute_a_over_r0)


def main():
    # Improved Wilson Nt=4
    iw_a_r0 = []
    iw_a_r0_err = []
    iw_r0_sqrts = []
    iw_r0_sqrts_err = []

    for row in TABLE2:
        a2s, a2s_err = row[7], row[8]
        r0a, r0a_err = row[13], row[14]

        ar, ar_err = compute_a_over_r0(r0a, r0a_err)
        rs, rs_err = compute_r0_sqrtsigma(r0a, r0a_err, a2s, a2s_err)

        iw_a_r0.append(ar)
        iw_a_r0_err.append(ar_err)
        iw_r0_sqrts.append(rs)
        iw_r0_sqrts_err.append(rs_err)

    # KS data
    ks_a_r0 = []
    ks_a_r0_err = []
    ks_r0_sqrts = []
    ks_r0_sqrts_err = []
    ks_Nt = []

    for row in TABLE3:
        a2s, a2s_err = row[7], row[8]
        r0a, r0a_err = row[13], row[14]
        Nt = row[15]

        ar, ar_err = compute_a_over_r0(r0a, r0a_err)
        rs, rs_err = compute_r0_sqrtsigma(r0a, r0a_err, a2s, a2s_err)

        ks_a_r0.append(ar)
        ks_a_r0_err.append(ar_err)
        ks_r0_sqrts.append(rs)
        ks_r0_sqrts_err.append(rs_err)
        ks_Nt.append(Nt)

    # Save data
    header = "a_over_r0,a_over_r0_err,r0_sqrtsigma,r0_sqrtsigma_err,action"
    lines = [header]
    for i in range(len(iw_a_r0)):
        lines.append(f"{iw_a_r0[i]:.8f},{iw_a_r0_err[i]:.8f},"
                     f"{iw_r0_sqrts[i]:.8f},{iw_r0_sqrts_err[i]:.8f},"
                     f"improved_wilson_Nt4")
    for i in range(len(ks_a_r0)):
        lines.append(f"{ks_a_r0[i]:.8f},{ks_a_r0_err[i]:.8f},"
                     f"{ks_r0_sqrts[i]:.8f},{ks_r0_sqrts_err[i]:.8f},"
                     f"KS_Nt{ks_Nt[i]}")
    with open("../data/fig15_r0_sigma.csv", "w") as f:
        f.write("\n".join(lines) + "\n")

    # Plot
    fig, ax = plt.subplots(figsize=(8, 6))

    ax.errorbar(iw_a_r0, iw_r0_sqrts, xerr=iw_a_r0_err, yerr=iw_r0_sqrts_err,
                fmt='o', color='blue', markersize=10, markerfacecolor='white',
                markeredgewidth=2, capsize=3,
                label=r'$N_t=4$ improved Wilson', zorder=5)

    ks4 = [i for i in range(len(ks_Nt)) if ks_Nt[i] == 4]
    ks6 = [i for i in range(len(ks_Nt)) if ks_Nt[i] == 6]

    if ks4:
        ax.errorbar([ks_a_r0[i] for i in ks4], [ks_r0_sqrts[i] for i in ks4],
                    xerr=[ks_a_r0_err[i] for i in ks4],
                    yerr=[ks_r0_sqrts_err[i] for i in ks4],
                    fmt='D', color='red', markersize=8, capsize=3,
                    label=r'$N_t=4$ KS', zorder=4)
    if ks6:
        ax.errorbar([ks_a_r0[i] for i in ks6], [ks_r0_sqrts[i] for i in ks6],
                    xerr=[ks_a_r0_err[i] for i in ks6],
                    yerr=[ks_r0_sqrts_err[i] for i in ks6],
                    fmt='s', color='green', markersize=8, capsize=3,
                    label=r'$N_t=6$ KS', zorder=4)

    # Continuum limit reference line
    ax.axhline(y=1.18, color='gray', linestyle='--', alpha=0.5,
               label=r'Continuum $r_0\sqrt{\sigma} \approx 1.18$')

    ax.set_xlabel(r'$a/r_0$', fontsize=14)
    ax.set_ylabel(r'$r_0\sqrt{\sigma}$', fontsize=14)
    ax.set_title(r'$r_0\sqrt{\sigma}$ vs $a/r_0$', fontsize=14)
    ax.legend(fontsize=11)
    ax.set_xlim(0.0, 0.8)
    ax.set_ylim(0.6, 1.6)

    plt.tight_layout()
    plt.savefig("../plots/fig15_r0_sigma.png", dpi=150)
    plt.close()
    print("Fig 15: r0*sqrt(sigma) saved.")


main()
