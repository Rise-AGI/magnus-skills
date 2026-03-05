"""
Figure 13: Tc/sqrt(sigma) and r0*Tc vs MPS/MV

Critical temperature scaled by (a) sqrt of string tension and
(b) inverse Sommer parameter vs. pseudoscalar/vector meson mass ratio.
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(".")))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from paper_data import (TABLE1, TABLE2, TABLE3, KS_MESON_DATA,
                        get_crossover_data, compute_tc_over_sqrtsigma,
                        compute_r0_tc, QUENCHED_TC_SQRTSIGMA_NT4)


def get_mps_mv_for_potential(table_potential, table_spectrum, is_ks=False):
    """Match potential fit data with spectroscopy to get MPS/MV."""
    results = []
    if is_ks:
        for i, row in enumerate(table_potential):
            if i < len(KS_MESON_DATA):
                mps_mv = KS_MESON_DATA[i][2]
                mps_mv_err = KS_MESON_DATA[i][3]
                results.append((mps_mv, mps_mv_err))
            else:
                results.append((None, None))
    else:
        for row in table_potential:
            beta_pot, kappa_pot = row[0], row[1]
            matched = False
            for spec in table_spectrum:
                if abs(spec[0] - beta_pot) < 0.01 and abs(spec[1] - kappa_pot) < 0.001:
                    results.append((spec[10], spec[11]))  # MPS/MV, err
                    matched = True
                    break
            if not matched:
                results.append((None, None))
    return results


def main():
    # Improved Wilson data
    mps_mv_iw = get_mps_mv_for_potential(TABLE2, TABLE1, is_ks=False)

    iw_tc_sqrts = []
    iw_tc_sqrts_err = []
    iw_r0tc = []
    iw_r0tc_err = []
    iw_mps_mv = []
    iw_mps_mv_err = []

    for i, row in enumerate(TABLE2):
        mm = mps_mv_iw[i]
        if mm[0] is None:
            continue
        a2s, a2s_err = row[7], row[8]
        r0a, r0a_err = row[13], row[14]

        tc_s, tc_s_err = compute_tc_over_sqrtsigma(a2s, a2s_err, Nt=4)
        r0t, r0t_err = compute_r0_tc(r0a, r0a_err, Nt=4)

        iw_mps_mv.append(mm[0])
        iw_mps_mv_err.append(mm[1])
        iw_tc_sqrts.append(tc_s)
        iw_tc_sqrts_err.append(tc_s_err)
        iw_r0tc.append(r0t)
        iw_r0tc_err.append(r0t_err)

    # KS data
    ks_tc_sqrts = []
    ks_tc_sqrts_err = []
    ks_r0tc = []
    ks_r0tc_err = []
    ks_mps_mv = []
    ks_mps_mv_err = []
    ks_Nt = []

    for i, row in enumerate(TABLE3):
        a2s, a2s_err = row[7], row[8]
        r0a, r0a_err = row[13], row[14]
        Nt = row[15]

        tc_s, tc_s_err = compute_tc_over_sqrtsigma(a2s, a2s_err, Nt=Nt)
        r0t, r0t_err = compute_r0_tc(r0a, r0a_err, Nt=Nt)

        if i < len(KS_MESON_DATA):
            mm, mm_err = KS_MESON_DATA[i][2], KS_MESON_DATA[i][3]
        else:
            continue

        ks_mps_mv.append(mm)
        ks_mps_mv_err.append(mm_err)
        ks_tc_sqrts.append(tc_s)
        ks_tc_sqrts_err.append(tc_s_err)
        ks_r0tc.append(r0t)
        ks_r0tc_err.append(r0t_err)
        ks_Nt.append(Nt)

    # Save data
    header = "MPS_MV,MPS_MV_err,Tc_sqrtsigma,Tc_sqrtsigma_err,r0_Tc,r0_Tc_err,action"
    lines = [header]
    for i in range(len(iw_mps_mv)):
        lines.append(f"{iw_mps_mv[i]:.8f},{iw_mps_mv_err[i]:.8f},"
                     f"{iw_tc_sqrts[i]:.8f},{iw_tc_sqrts_err[i]:.8f},"
                     f"{iw_r0tc[i]:.8f},{iw_r0tc_err[i]:.8f},"
                     f"improved_wilson_Nt4")
    for i in range(len(ks_mps_mv)):
        lines.append(f"{ks_mps_mv[i]:.8f},{ks_mps_mv_err[i]:.8f},"
                     f"{ks_tc_sqrts[i]:.8f},{ks_tc_sqrts_err[i]:.8f},"
                     f"{ks_r0tc[i]:.8f},{ks_r0tc_err[i]:.8f},"
                     f"KS_Nt{ks_Nt[i]}")
    with open("../data/fig13_tc_scaling.csv", "w") as f:
        f.write("\n".join(lines) + "\n")

    # Plot: two panels
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # Panel (a): Tc/sqrt(sigma)
    ax1.errorbar(iw_mps_mv, iw_tc_sqrts, xerr=iw_mps_mv_err, yerr=iw_tc_sqrts_err,
                 fmt='o', color='blue', markersize=10, markerfacecolor='white',
                 markeredgewidth=2, capsize=3,
                 label=r'$N_t=4$ improved Wilson', zorder=5)

    ks4 = [i for i in range(len(ks_Nt)) if ks_Nt[i] == 4]
    ks6 = [i for i in range(len(ks_Nt)) if ks_Nt[i] == 6]

    if ks4:
        ax1.errorbar([ks_mps_mv[i] for i in ks4], [ks_tc_sqrts[i] for i in ks4],
                     xerr=[ks_mps_mv_err[i] for i in ks4],
                     yerr=[ks_tc_sqrts_err[i] for i in ks4],
                     fmt='D', color='red', markersize=8, capsize=3,
                     label=r'$N_t=4$ KS', zorder=4)
    if ks6:
        ax1.errorbar([ks_mps_mv[i] for i in ks6], [ks_tc_sqrts[i] for i in ks6],
                     xerr=[ks_mps_mv_err[i] for i in ks6],
                     yerr=[ks_tc_sqrts_err[i] for i in ks6],
                     fmt='s', color='green', markersize=8, capsize=3,
                     label=r'$N_t=6$ KS', zorder=4)

    # Quenched arrow
    ax1.annotate('', xy=(1.02, QUENCHED_TC_SQRTSIGMA_NT4),
                 xytext=(1.08, QUENCHED_TC_SQRTSIGMA_NT4),
                 arrowprops=dict(arrowstyle='->', color='black', lw=1.5))
    ax1.text(1.09, QUENCHED_TC_SQRTSIGMA_NT4, r'quenched $N_t=4$', fontsize=9,
             va='center')

    ax1.set_xlabel(r'$M_{PS}/M_V$', fontsize=14)
    ax1.set_ylabel(r'$T_c/\sqrt{\sigma}$', fontsize=14)
    ax1.set_title('(a)', fontsize=14)
    ax1.legend(fontsize=10)
    ax1.set_xlim(0.0, 1.2)
    ax1.set_ylim(0.3, 0.7)

    # Panel (b): r0*Tc
    ax2.errorbar(iw_mps_mv, iw_r0tc, xerr=iw_mps_mv_err, yerr=iw_r0tc_err,
                 fmt='o', color='blue', markersize=10, markerfacecolor='white',
                 markeredgewidth=2, capsize=3,
                 label=r'$N_t=4$ improved Wilson', zorder=5)

    if ks4:
        ax2.errorbar([ks_mps_mv[i] for i in ks4], [ks_r0tc[i] for i in ks4],
                     xerr=[ks_mps_mv_err[i] for i in ks4],
                     yerr=[ks_r0tc_err[i] for i in ks4],
                     fmt='D', color='red', markersize=8, capsize=3,
                     label=r'$N_t=4$ KS', zorder=4)
    if ks6:
        ax2.errorbar([ks_mps_mv[i] for i in ks6], [ks_r0tc[i] for i in ks6],
                     xerr=[ks_mps_mv_err[i] for i in ks6],
                     yerr=[ks_r0tc_err[i] for i in ks6],
                     fmt='s', color='green', markersize=8, capsize=3,
                     label=r'$N_t=6$ KS', zorder=4)

    ax2.set_xlabel(r'$M_{PS}/M_V$', fontsize=14)
    ax2.set_ylabel(r'$r_0 T_c$', fontsize=14)
    ax2.set_title('(b)', fontsize=14)
    ax2.legend(fontsize=10)
    ax2.set_xlim(0.0, 1.0)
    ax2.set_ylim(0.3, 0.65)

    plt.tight_layout()
    plt.savefig("../plots/fig13_tc_scaling.png", dpi=150)
    plt.close()
    print("Fig 13: Tc scaling saved.")


main()
