"""
Fig 3b: Current density profiles before and after the current quench.

Generates idealized pre-quench and post-quench profiles based on the
0D current evolution and known physics of peaked RE current formation.
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from re_fluid import PlasmaParams, CurrentQuench0D
from scipy.constants import pi

def main():
    params = PlasmaParams(
        R=10.0, a=1.0, B_phi0=1.0,
        T0_keV=1.7, n0=1e20,
        Z_eff=1.0, eta0=1.1e-7,
        lnLambda=15.0, I_p_MA=0.67
    )

    # Run 0D model to get final currents
    cq = CurrentQuench0D(params, T_final_eV=5.0, tau_TQ_ms=5.0, I_r_seed_frac=1e-3)
    times, I_total, I_re, I_th = cq.run(t_final_ms=70.0, n_points=500)

    I_total_final = I_total[-1] * 1e6  # back to Amps
    I_re_final = I_re[-1] * 1e6
    I_th_final = I_th[-1] * 1e6

    # Radial coordinate
    r_over_a = np.linspace(0, 1, 200)

    # Pre-quench: parabolic j profile
    j0_pre = 2 * params.I_p / (pi * params.a**2)
    j_pre = j0_pre * (1 - r_over_a**2)

    # Post-quench profiles
    # RE current: centrally peaked (resistive diffusion of E toward center)
    # Use a peaked profile: j_r ~ (1 - (r/a)^2)^2
    j_r_profile = (1 - r_over_a**2)**2
    j_r_norm = 2 * pi * np.trapezoid(j_r_profile * r_over_a * params.a**2, r_over_a)
    j_r_post = j_r_profile * I_re_final / j_r_norm

    # Thermal current: broader, reduced
    j_th_profile = (1 - r_over_a**2)
    j_th_norm = 2 * pi * np.trapezoid(j_th_profile * r_over_a * params.a**2, r_over_a)
    j_th_post = j_th_profile * I_th_final / j_th_norm

    j_total_post = j_r_post + j_th_post

    data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data')
    os.makedirs(data_dir, exist_ok=True)
    plots_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'plots')
    os.makedirs(plots_dir, exist_ok=True)

    header = "r_over_a,j_prequench_Am2,j_total_postquench_Am2,j_re_postquench_Am2,j_thermal_postquench_Am2"
    data = np.column_stack([r_over_a, j_pre, j_total_post, j_r_post, j_th_post])
    np.savetxt(os.path.join(data_dir, 'fig3b_profiles.csv'),
               data, delimiter=',', header=header, comments='', fmt='%.8e')

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(r_over_a, j_pre * 1e-6, 'k-', linewidth=2, label='Pre-quench (total)')
    ax.plot(r_over_a, j_total_post * 1e-6, 'b-', linewidth=2, label='Post-quench (total)')
    ax.plot(r_over_a, j_r_post * 1e-6, 'r--', linewidth=2, label='Post-quench (RE)')
    ax.plot(r_over_a, j_th_post * 1e-6, 'g:', linewidth=2, label='Post-quench (thermal)')
    ax.set_xlabel('r/a', fontsize=12)
    ax.set_ylabel(r'$j$ [MA/m$^2$]', fontsize=12)
    ax.set_title('Fig 3b: Current density profiles', fontsize=13)
    ax.legend(fontsize=10)
    ax.set_xlim(0, 1)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(plots_dir, 'fig3b_profiles.png'), dpi=150)
    plt.close(fig)
    print("Saved fig3b_profiles.csv and fig3b_profiles.png")


if __name__ == '__main__':
    main()
