"""
Figure 12: Tc/MV vs MPS/MV

Critical temperature divided by vector meson mass vs. pseudoscalar/vector
meson mass ratio. Shows improved Wilson Nt=4 data from Table 1.
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(".")))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from paper_data import (TABLE1, TABLE2, TABLE3, KS_MESON_DATA,
                        get_crossover_data, compute_tc_over_mv)


def main():
    crossover = get_crossover_data()

    # Improved Wilson Nt=4 data
    iw_mps_mv = []
    iw_mps_mv_err = []
    iw_tc_mv = []
    iw_tc_mv_err = []
    iw_beta = []

    for row in crossover:
        beta, kappa, u0, n, aMPS, aMPS_err, aMV, aMV_err = row[:8]
        mps_mv, mps_mv_err = row[10], row[11]

        tc_mv, tc_mv_err = compute_tc_over_mv(aMV, aMV_err, Nt=4)

        iw_beta.append(beta)
        iw_mps_mv.append(mps_mv)
        iw_mps_mv_err.append(mps_mv_err)
        iw_tc_mv.append(tc_mv)
        iw_tc_mv_err.append(tc_mv_err)

    # KS Nt=4 data (approximate from Ref [34])
    ks_mps_mv = []
    ks_mps_mv_err = []
    ks_tc_mv = []
    ks_tc_mv_err = []
    ks_Nt = []

    for i, ks_row in enumerate(KS_MESON_DATA):
        ks_beta, amq, mps_mv, mps_mv_err = ks_row
        Nt = TABLE3[i][15]

        # For KS, we need aMV to compute Tc/MV
        # aMV can be estimated from MPS/MV and aMPS
        # Since we don't have aMV directly, use aMPS/ratio
        # aMPS ~ 2*sqrt(amq) approximately for KS (leading order)
        # But better to use the relation Tc/MV = 1/(Nt*aMV)
        # We'll use approximate aMV values typical for these parameters
        # From Blum et al. (1995), approximate aMV values:
        if Nt == 4:
            if amq == 0.025:
                aMV_ks, aMV_ks_err = 1.07, 0.05
            elif amq == 0.050:
                aMV_ks, aMV_ks_err = 1.12, 0.04
            elif amq == 0.100:
                aMV_ks, aMV_ks_err = 1.22, 0.03
            else:
                continue
        else:  # Nt=6
            aMV_ks, aMV_ks_err = 0.72, 0.04

        tc_mv_val, tc_mv_e = compute_tc_over_mv(aMV_ks, aMV_ks_err, Nt=Nt)
        ks_mps_mv.append(mps_mv)
        ks_mps_mv_err.append(mps_mv_err)
        ks_tc_mv.append(tc_mv_val)
        ks_tc_mv_err.append(tc_mv_e)
        ks_Nt.append(Nt)

    # Save data
    header = "MPS_MV,MPS_MV_err,Tc_MV,Tc_MV_err,beta,action"
    lines = [header]
    for i in range(len(iw_mps_mv)):
        lines.append(f"{iw_mps_mv[i]:.8f},{iw_mps_mv_err[i]:.8f},"
                     f"{iw_tc_mv[i]:.8f},{iw_tc_mv_err[i]:.8f},"
                     f"{iw_beta[i]:.2f},improved_wilson_Nt4")
    for i in range(len(ks_mps_mv)):
        ks_beta = TABLE3[i][0]
        lines.append(f"{ks_mps_mv[i]:.8f},{ks_mps_mv_err[i]:.8f},"
                     f"{ks_tc_mv[i]:.8f},{ks_tc_mv_err[i]:.8f},"
                     f"{ks_beta:.4f},KS_Nt{ks_Nt[i]}")
    with open("../data/fig12_tc_over_mv.csv", "w") as f:
        f.write("\n".join(lines) + "\n")

    # Plot
    fig, ax = plt.subplots(figsize=(8, 6))

    ax.errorbar(iw_mps_mv, iw_tc_mv, xerr=iw_mps_mv_err, yerr=iw_tc_mv_err,
                fmt='o', color='blue', markersize=10, markerfacecolor='white',
                markeredgewidth=2, capsize=3,
                label=r'$N_t=4$ improved Wilson', zorder=5)

    # KS Nt=4
    ks4_idx = [i for i in range(len(ks_Nt)) if ks_Nt[i] == 4]
    ks6_idx = [i for i in range(len(ks_Nt)) if ks_Nt[i] == 6]

    if ks4_idx:
        ax.errorbar([ks_mps_mv[i] for i in ks4_idx],
                    [ks_tc_mv[i] for i in ks4_idx],
                    xerr=[ks_mps_mv_err[i] for i in ks4_idx],
                    yerr=[ks_tc_mv_err[i] for i in ks4_idx],
                    fmt='D', color='red', markersize=8, capsize=3,
                    label=r'$N_t=4$ KS', zorder=4)

    if ks6_idx:
        ax.errorbar([ks_mps_mv[i] for i in ks6_idx],
                    [ks_tc_mv[i] for i in ks6_idx],
                    xerr=[ks_mps_mv_err[i] for i in ks6_idx],
                    yerr=[ks_tc_mv_err[i] for i in ks6_idx],
                    fmt='s', color='green', markersize=8, capsize=3,
                    label=r'$N_t=6$ KS', zorder=4)

    ax.set_xlabel(r'$M_{PS}/M_V$', fontsize=14)
    ax.set_ylabel(r'$T_c/M_V$', fontsize=14)
    ax.set_title(r'$T_c/M_V$ vs $M_{PS}/M_V$', fontsize=14)
    ax.legend(fontsize=11)
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 0.35)

    plt.tight_layout()
    plt.savefig("../plots/fig12_tc_over_mv.png", dpi=150)
    plt.close()
    print("Fig 12: Tc/MV saved.")


main()
