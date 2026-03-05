"""
Fig 6a: Evolution of plasma currents during the axisymmetric VDE simulation.

Reproduces the ITER VDE current evolution from Section 4.
I_p = 14.5 MA, B = 4.8 T, n_e = 5e19 m^-3, T ~ 2.35 eV post-TQ.
Shows I_total and I_r for cases with and without RE generation.
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(os.path.join(os.getcwd(), __file__))) if '__file__' not in dir() else os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from re_fluid import VDE1D

def main():
    data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data')
    os.makedirs(data_dir, exist_ok=True)
    plots_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'plots')
    os.makedirs(plots_dir, exist_ok=True)

    # Case with REs
    print("Running VDE simulation with REs...")
    vde_re = VDE1D(
        I_p_MA=14.5, B_phi=4.8, n_e=5e19,
        eta_axis=1.24e-4, a=2.0, R=6.2,
        T_eV=2.35, Z_eff=1.0, lnLambda=15.0,
        I_r_seed_fraction=1e-3
    )
    t_re, I_total_re, I_r_re = vde_re.run(t_final_ms=12.0, dt_ms=0.001)
    print(f"  With REs: I_total final = {I_total_re[-1]:.2f} MA, I_r final = {I_r_re[-1]:.2f} MA")

    # Baseline without REs (just resistive decay)
    print("Running VDE simulation without REs (baseline)...")
    vde_base = VDE1D(
        I_p_MA=14.5, B_phi=4.8, n_e=5e19,
        eta_axis=1.24e-4, a=2.0, R=6.2,
        T_eV=2.35, Z_eff=1.0, lnLambda=15.0,
        I_r_seed_fraction=0.0  # No REs
    )
    t_base, I_total_base, I_r_base = vde_base.run(t_final_ms=12.0, dt_ms=0.001)
    print(f"  Baseline: I_total final = {I_total_base[-1]:.2f} MA")

    # Save data
    # Interpolate to common time grid
    t_common = np.linspace(0, 12, 500)
    I_total_re_interp = np.interp(t_common, t_re, I_total_re)
    I_r_re_interp = np.interp(t_common, t_re, I_r_re)
    I_total_base_interp = np.interp(t_common, t_base, I_total_base)
    I_th_re_interp = I_total_re_interp - I_r_re_interp

    header = "time_ms,I_total_with_RE_MA,I_r_MA,I_th_with_RE_MA,I_total_baseline_MA"
    data = np.column_stack([t_common, I_total_re_interp, I_r_re_interp,
                           I_th_re_interp, I_total_base_interp])
    np.savetxt(os.path.join(data_dir, 'fig6a_vde_currents.csv'),
               data, delimiter=',', header=header, comments='', fmt='%.8e')

    # Plot
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(t_common, I_total_re_interp, 'b-', linewidth=2,
            label=r'$I_{\mathrm{total}}$ (with REs)')
    ax.plot(t_common, I_r_re_interp, 'r-', linewidth=2,
            label=r'$I_r$ (RE current)')
    ax.plot(t_common, I_th_re_interp, 'g--', linewidth=2,
            label=r'$I_{\mathrm{th}}$ (thermal, with REs)')
    ax.plot(t_common, I_total_base_interp, 'k:', linewidth=2,
            label=r'$I_{\mathrm{total}}$ (baseline, no REs)')
    ax.set_xlabel('Time [ms]', fontsize=12)
    ax.set_ylabel('Current [MA]', fontsize=12)
    ax.set_title('Fig 6a: ITER VDE current evolution', fontsize=13)
    ax.legend(fontsize=10)
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 16)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(plots_dir, 'fig6a_vde_currents.png'), dpi=150)
    plt.close(fig)
    print("Saved fig6a_vde_currents.csv and fig6a_vde_currents.png")


if __name__ == '__main__':
    main()
