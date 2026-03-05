"""
Figure 7: Klein tunneling transmission T(phi).

Top panel: D = 110 nm, V0 = 200 meV (dashed) and 285 meV (solid).
Bottom panel: D = 50 nm, same V0 values.
E = 80 meV, lambda = 50 nm (kF = 2*pi/lambda).
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(os.path.realpath(sys.argv[0]))))
from graphene_params import klein_tunneling_T


def compute_klein_tunneling():
    # Parameters from paper
    E = 0.080     # eV
    V0_1 = 0.200  # eV
    V0_2 = 0.285  # eV
    D1 = 110e-9   # meters
    D2 = 50e-9    # meters
    vF = 1e6      # m/s

    phi_deg = np.linspace(-89.5, 89.5, 500)
    phi_rad = np.radians(phi_deg)

    # Top panel: D = 110 nm
    T_D1_V1 = klein_tunneling_T(phi_rad, E, V0_1, D1, vF)
    T_D1_V2 = klein_tunneling_T(phi_rad, E, V0_2, D1, vF)

    # Bottom panel: D = 50 nm
    T_D2_V1 = klein_tunneling_T(phi_rad, E, V0_1, D2, vF)
    T_D2_V2 = klein_tunneling_T(phi_rad, E, V0_2, D2, vF)

    return phi_deg, T_D1_V1, T_D1_V2, T_D2_V1, T_D2_V2


def save_data(phi_deg, T_D1_V1, T_D1_V2, T_D2_V1, T_D2_V2):
    data = np.column_stack([phi_deg, T_D1_V1, T_D1_V2, T_D2_V1, T_D2_V2])
    header = "phi_deg,T_D110nm_V200meV,T_D110nm_V285meV,T_D50nm_V200meV,T_D50nm_V285meV"
    np.savetxt("../data/fig7_klein_tunneling.csv", data, delimiter=",",
               header=header, fmt="%.8e", comments="")


def plot_figure(phi_deg, T_D1_V1, T_D1_V2, T_D2_V1, T_D2_V2):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8))

    # Top: D = 110 nm
    ax1.plot(phi_deg, T_D1_V1, 'b--', linewidth=1.2, label=r'$V_0 = 200$ meV')
    ax1.plot(phi_deg, T_D1_V2, 'r-', linewidth=1.2, label=r'$V_0 = 285$ meV')
    ax1.set_ylabel(r'$T(\phi)$')
    ax1.set_ylim(0, 1.05)
    ax1.set_xlim(-90, 90)
    ax1.set_title('D = 110 nm')
    ax1.legend()
    ax1.set_xticks(np.arange(-90, 91, 15))

    # Bottom: D = 50 nm
    ax2.plot(phi_deg, T_D2_V1, 'b--', linewidth=1.2, label=r'$V_0 = 200$ meV')
    ax2.plot(phi_deg, T_D2_V2, 'r-', linewidth=1.2, label=r'$V_0 = 285$ meV')
    ax2.set_xlabel(r'angle $\phi$ (degrees)')
    ax2.set_ylabel(r'$T(\phi)$')
    ax2.set_ylim(0, 1.05)
    ax2.set_xlim(-90, 90)
    ax2.set_title('D = 50 nm')
    ax2.legend()
    ax2.set_xticks(np.arange(-90, 91, 15))

    plt.tight_layout()
    plt.savefig("../plots/fig7_klein_tunneling.png", dpi=150, bbox_inches='tight')
    plt.close()


if True:
    phi_deg, T_D1_V1, T_D1_V2, T_D2_V1, T_D2_V2 = compute_klein_tunneling()
    save_data(phi_deg, T_D1_V1, T_D1_V2, T_D2_V1, T_D2_V2)
    plot_figure(phi_deg, T_D1_V1, T_D1_V2, T_D2_V1, T_D2_V2)
    print("Figure 7 done: Klein tunneling saved.")
