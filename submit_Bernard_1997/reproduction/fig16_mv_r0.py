"""
Figure 16: MV*r0 vs MPS/MV

Vector meson mass times Sommer parameter vs. pseudoscalar/vector meson
mass ratio for improved Wilson and KS fermions.
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(".")))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from paper_data import (TABLE1, TABLE2, TABLE3, KS_MESON_DATA,
                        get_crossover_data, compute_mv_r0)


def main():
    # Improved Wilson: match Table 1 (crossover points) with Table 2
    iw_mps_mv = []
    iw_mps_mv_err = []
    iw_mv_r0 = []
    iw_mv_r0_err = []

    for pot_row in TABLE2:
        beta_pot, kappa_pot = pot_row[0], pot_row[1]
        r0a, r0a_err = pot_row[13], pot_row[14]

        # Find matching spectroscopy entry
        for spec_row in TABLE1:
            if abs(spec_row[0] - beta_pot) < 0.01 and abs(spec_row[1] - kappa_pot) < 0.001:
                aMV = spec_row[6]
                aMV_err = spec_row[7]
                mps_mv = spec_row[10]
                mps_mv_err = spec_row[11]

                mv_r0, mv_r0_err = compute_mv_r0(aMV, aMV_err, r0a, r0a_err)

                iw_mps_mv.append(mps_mv)
                iw_mps_mv_err.append(mps_mv_err)
                iw_mv_r0.append(mv_r0)
                iw_mv_r0_err.append(mv_r0_err)
                break

    # KS data
    ks_mps_mv = []
    ks_mps_mv_err = []
    ks_mv_r0 = []
    ks_mv_r0_err = []
    ks_Nt = []

    for i, row in enumerate(TABLE3):
        r0a, r0a_err = row[13], row[14]
        Nt = row[15]

        if i >= len(KS_MESON_DATA):
            continue
        mm, mm_err = KS_MESON_DATA[i][2], KS_MESON_DATA[i][3]
        amq = row[1]

        # Approximate aMV for KS
        if Nt == 4:
            if amq == 0.025:
                aMV_ks, aMV_ks_err = 1.07, 0.05
            elif amq == 0.050:
                aMV_ks, aMV_ks_err = 1.12, 0.04
            elif amq == 0.100:
                aMV_ks, aMV_ks_err = 1.22, 0.03
            else:
                continue
        else:
            aMV_ks, aMV_ks_err = 0.72, 0.04

        mv_r0, mv_r0_err = compute_mv_r0(aMV_ks, aMV_ks_err, r0a, r0a_err)
        ks_mps_mv.append(mm)
        ks_mps_mv_err.append(mm_err)
        ks_mv_r0.append(mv_r0)
        ks_mv_r0_err.append(mv_r0_err)
        ks_Nt.append(Nt)

    # Save data
    header = "MPS_MV,MPS_MV_err,MV_r0,MV_r0_err,action"
    lines = [header]
    for i in range(len(iw_mps_mv)):
        lines.append(f"{iw_mps_mv[i]:.8f},{iw_mps_mv_err[i]:.8f},"
                     f"{iw_mv_r0[i]:.8f},{iw_mv_r0_err[i]:.8f},"
                     f"improved_wilson_Nt4")
    for i in range(len(ks_mps_mv)):
        lines.append(f"{ks_mps_mv[i]:.8f},{ks_mps_mv_err[i]:.8f},"
                     f"{ks_mv_r0[i]:.8f},{ks_mv_r0_err[i]:.8f},"
                     f"KS_Nt{ks_Nt[i]}")
    with open("../data/fig16_mv_r0.csv", "w") as f:
        f.write("\n".join(lines) + "\n")

    # Plot
    fig, ax = plt.subplots(figsize=(8, 6))

    ax.errorbar(iw_mps_mv, iw_mv_r0, xerr=iw_mps_mv_err, yerr=iw_mv_r0_err,
                fmt='o', color='blue', markersize=10, markerfacecolor='white',
                markeredgewidth=2, capsize=3,
                label=r'$N_t=4$ improved Wilson', zorder=5)

    ks4 = [i for i in range(len(ks_Nt)) if ks_Nt[i] == 4]
    ks6 = [i for i in range(len(ks_Nt)) if ks_Nt[i] == 6]

    if ks4:
        ax.errorbar([ks_mps_mv[i] for i in ks4], [ks_mv_r0[i] for i in ks4],
                    xerr=[ks_mps_mv_err[i] for i in ks4],
                    yerr=[ks_mv_r0_err[i] for i in ks4],
                    fmt='D', color='red', markersize=8, capsize=3,
                    label=r'$N_t=4$ KS', zorder=4)
    if ks6:
        ax.errorbar([ks_mps_mv[i] for i in ks6], [ks_mv_r0[i] for i in ks6],
                    xerr=[ks_mps_mv_err[i] for i in ks6],
                    yerr=[ks_mv_r0_err[i] for i in ks6],
                    fmt='s', color='green', markersize=8, capsize=3,
                    label=r'$N_t=6$ KS', zorder=4)

    # Physical rho mass times r0 for reference
    # r0 = 0.49 fm, M_rho = 770 MeV => M_rho * r0 ~ 1.91
    ax.axhline(y=1.91, color='gray', linestyle='--', alpha=0.5,
               label=r'Physical $M_\rho r_0 \approx 1.91$')

    ax.set_xlabel(r'$M_{PS}/M_V$', fontsize=14)
    ax.set_ylabel(r'$M_V r_0$', fontsize=14)
    ax.set_title(r'$M_V r_0$ vs $M_{PS}/M_V$', fontsize=14)
    ax.legend(fontsize=11)
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 6.0)

    plt.tight_layout()
    plt.savefig("../plots/fig16_mv_r0.png", dpi=150)
    plt.close()
    print("Fig 16: MV*r0 saved.")


main()
