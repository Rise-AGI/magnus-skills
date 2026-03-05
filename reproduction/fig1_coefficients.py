"""
Figure 1 data and scattering coefficients.

Computes:
- T, R, P for specific parameters from paper text (wR/c = 800, theta = 90 deg)
- Field profiles in flat and cylindrical geometries
- Coupling efficiency Delta^2 vs wR/c
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__) if '__file__' in dir() else '.')

import numpy as np
import mpmath
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from spp_core import (solve_mode_index, transmittance, reflectance,
                       upper_bound_transmittance, flat_surface_loss,
                       decay_constants, spp_wavevector, c)

eps_o = 1.0
theta = np.pi / 2

data_dir = os.path.join(os.path.dirname(__file__) if '__file__' in dir() else '.', '..', 'data')
plots_dir = os.path.join(os.path.dirname(__file__) if '__file__' in dir() else '.', '..', 'plots')
os.makedirs(data_dir, exist_ok=True)
os.makedirs(plots_dir, exist_ok=True)

# ---- Scattering coefficients for paper's specific parameters ----
omega_R_over_c = 800
omega = 2 * np.pi * c / (600e-9)
R = omega_R_over_c * c / omega

results = []
for label, eps_i in [("lossless", -15.0 + 0j), ("lossy", -15.0 + 0.5j)]:
    k = spp_wavevector(omega, eps_i, eps_o)
    m = solve_mode_index(omega, R, eps_i, eps_o)
    T = transmittance(m, k, R, theta)
    R_c = reflectance(m, k, R, theta)
    Tu = upper_bound_transmittance(m, theta)
    flat = flat_surface_loss(omega, eps_i, eps_o, R * theta)

    print(f"\n=== {label} (eps_i = {eps_i}) ===")
    print(f"m = {m:.6f}")
    print(f"T = {T:.6e}, R = {R_c:.6e}, A = {1-T-R_c:.6e}")
    print(f"Tu = {Tu:.6e}, Flat = {flat:.6e}")

    results.append({'label': label, 'eps_i': eps_i, 'T': T, 'R_coeff': R_c,
                     'A': 1-T-R_c, 'Tu': Tu, 'flat': flat, 'm': m})

csv_path = os.path.join(data_dir, 'scattering_coefficients.csv')
with open(csv_path, 'w') as f:
    f.write('case,eps_i_real,eps_i_imag,T,R_coeff,A,Tu,flat_surface,m_real,m_imag\n')
    for r in results:
        f.write(f"{r['label']},{r['eps_i'].real:.1f},{r['eps_i'].imag:.1f},"
                f"{r['T']:.8e},{r['R_coeff']:.8e},{r['A']:.8e},{r['Tu']:.8e},"
                f"{r['flat']:.8e},{r['m'].real:.8e},{r['m'].imag:.8e}\n")
print(f"\nCoefficients saved to {csv_path}")

# ---- Coupling efficiency Delta^2 ----
print("\n=== Coupling efficiency ===")
mpmath.mp.dps = 15

orc_values = [50, 100, 150, 200, 300, 400, 600, 800, 1000]
csv_path2 = os.path.join(data_dir, 'coupling_efficiency.csv')
delta2_list = []

for orc in orc_values:
    R_val = orc * c / omega
    eps_i_val = -15.0 + 0j
    k_val = spp_wavevector(omega, eps_i_val, eps_o)
    m_val = solve_mode_index(omega, R_val, eps_i_val, eps_o)

    ko = omega * np.sqrt(eps_o + 0j) / c
    gamma_i, gamma_o = decay_constants(omega, eps_i_val, eps_o)
    eta = 3.0

    # Compute Delta^2 (Eq. 5) using mpmath for Hankel functions
    Nr = 200
    r_max_val = R_val + eta / gamma_o.real
    r_arr = np.linspace(R_val, r_max_val, Nr)

    m_mp = mpmath.mpc(m_val.real, m_val.imag)
    ko_R_mp = mpmath.mpc(float((ko * R_val).real), float((ko * R_val).imag))

    H_m_R = mpmath.besselj(m_mp, ko_R_mp) + 1j * mpmath.bessely(m_mp, ko_R_mp)

    cyl_field = np.zeros(Nr, dtype=complex)
    for ir, r in enumerate(r_arr):
        ko_r_mp = mpmath.mpc(float((ko * r).real), float((ko * r).imag))
        H_m_r = mpmath.besselj(m_mp, ko_r_mp) + 1j * mpmath.bessely(m_mp, ko_r_mp)
        cyl_field[ir] = complex(H_m_r / H_m_R)

    flat_field = np.exp(-gamma_o * (r_arr - R_val))
    diff = cyl_field - flat_field
    num = np.trapz(np.abs(diff)**2, r_arr)
    den = np.trapz(np.abs(flat_field)**2, r_arr)
    D2 = (num / den).real

    delta2_list.append(D2)
    print(f"  wR/c = {orc}: Delta^2 = {D2:.6f}")

with open(csv_path2, 'w') as f:
    f.write('omega_R_over_c,Delta2\n')
    for orc, d2 in zip(orc_values, delta2_list):
        f.write(f"{orc},{d2:.8e}\n")
print(f"Coupling efficiency saved to {csv_path2}")

# ---- Field profiles ----
print("\n=== Field profiles ===")
eps_i_fp = -15.0 + 0.5j
m_fp = solve_mode_index(omega, R, eps_i_fp, eps_o)
ki = omega * np.sqrt(eps_i_fp + 0j) / c
ko = omega * np.sqrt(eps_o + 0j) / c
gamma_i, gamma_o = decay_constants(omega, eps_i_fp, eps_o)

r_inner = np.linspace(R * 0.92, R, 100)
r_outer = np.linspace(R, R * 1.08, 100)

flat_inner = np.exp(gamma_i.real * (r_inner - R))
flat_outer = np.exp(-gamma_o.real * (r_outer - R))

m_mp = mpmath.mpc(m_fp.real, m_fp.imag)

# Cylindrical field inside (Bessel J)
ki_R_mp = mpmath.mpc(float((ki * R).real), float((ki * R).imag))
J_m_R = mpmath.besselj(m_mp, ki_R_mp)
cyl_inner = np.zeros(len(r_inner))
for ir, r in enumerate(r_inner):
    ki_r_mp = mpmath.mpc(float((ki * r).real), float((ki * r).imag))
    J_val = mpmath.besselj(m_mp, ki_r_mp)
    cyl_inner[ir] = abs(complex(J_val / J_m_R))

# Cylindrical field outside (Hankel H1)
ko_R_mp = mpmath.mpc(float((ko * R).real), float((ko * R).imag))
H_m_R = mpmath.besselj(m_mp, ko_R_mp) + 1j * mpmath.bessely(m_mp, ko_R_mp)
cyl_outer = np.zeros(len(r_outer))
for ir, r in enumerate(r_outer):
    ko_r_mp = mpmath.mpc(float((ko * r).real), float((ko * r).imag))
    H_val = mpmath.besselj(m_mp, ko_r_mp) + 1j * mpmath.bessely(m_mp, ko_r_mp)
    cyl_outer[ir] = abs(complex(H_val / H_m_R))

r_all = np.concatenate([r_inner, r_outer])
r_norm = (r_all - R) / R
flat_all = np.concatenate([np.abs(flat_inner), np.abs(flat_outer)])
cyl_all = np.concatenate([cyl_inner, cyl_outer])

csv_path3 = os.path.join(data_dir, 'fig1_field_profiles.csv')
data_mat = np.column_stack([r_norm, flat_all, cyl_all])
np.savetxt(csv_path3, data_mat, delimiter=',',
           header='r_minus_R_over_R,flat_spp_field,cylindrical_mode_field',
           comments='', fmt='%.8e')
print(f"Field profiles saved to {csv_path3}")

fig, ax = plt.subplots(figsize=(8, 5))
n = len(r_inner)
ax.plot(r_norm[:n], flat_all[:n], 'b-', label='Flat SPP (Region I)', linewidth=1.5)
ax.plot(r_norm[:n], cyl_all[:n], 'r--', label='Cylindrical mode (Region II)', linewidth=1.5)
ax.plot(r_norm[n:], flat_all[n:], 'b-', linewidth=1.5)
ax.plot(r_norm[n:], cyl_all[n:], 'r--', linewidth=1.5)
ax.axvline(0, color='k', linewidth=0.5, linestyle=':')
ax.set_xlabel(r'$(r - R)/R$', fontsize=12)
ax.set_ylabel('Field amplitude (normalized)', fontsize=12)
ax.set_title(r'SPP Field Profiles ($\omega R/c = 800$, Silver-Air)')
ax.legend(fontsize=11)
ax.grid(True, alpha=0.3)
plt.savefig(os.path.join(plots_dir, 'fig1_field_profiles.png'), dpi=150, bbox_inches='tight')
print("Field profile plot saved")
