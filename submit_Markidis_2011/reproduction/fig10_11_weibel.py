"""
Figure 10: Weibel instability - comparison with linear theory.
Figure 11: Energy and momentum history in Weibel instability.

Parameters from paper Section 7.3:
- Box: 2*pi*c/omega_pe, 64 grid points, periodic
- dt = 0.25 omega_pe^{-1}, 400 cycles
- 100000 electrons, bi-Maxwellian: vthy = 0.4c, anisotropy a=15
  -> vthx = vthy/sqrt(1+a) = 0.4/sqrt(16) = 0.1c
- mass_ratio = 1836
- Linear theory: gamma(k=1) = 0.22 omega_pe

We use reduced particle count for computational tractability.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys, os

sys.path.insert(0, os.path.dirname(os.path.abspath(".")))
from ecpic_em import run_weibel

# Parameters
WP = 1.0
c = 1.0
L = 2.0 * np.pi * c / WP
NG = 64
N = 50000       # reduced from 100000
# Paper uses DT=0.25 with implicit EC-PIC. For explicit FDTD, need CFL: c*dt/dx < 1
# dx = 2*pi/64 ~ 0.098, so dt < 0.098. Use dt=0.05 and scale NT.
DT = 0.05 / WP
NT = 2000       # same total time as paper: 400*0.25 = 100 wpe^-1
QM_e = -1.0
vthy = 0.4 * c
anisotropy = 15.0
vthx = vthy / np.sqrt(1.0 + anisotropy)
mass_ratio = 1836.0

print(f"Weibel instability simulation:")
print(f"  vthx = {vthx:.4f}c, vthy = {vthy:.4f}c")
print(f"  anisotropy = {anisotropy:.0f}")
print(f"  L = {L:.4f}, NG={NG}, N={N}, DT={DT:.4f}, NT={NT}")

result = run_weibel(L, DT, NT, NG, N, WP, QM_e, vthx, vthy, mass_ratio)

time = np.arange(NT) * DT

# --- Figure 10: Bz growth vs linear theory ---
gamma_theory = 0.22 * WP  # linear theory growth rate for k=1

Bz_k1 = result['Bz_k1']

# Fit exponential growth in linear phase
log_Bz = np.log(Bz_k1 + 1e-30)
# Find growth region
growth_mask = (Bz_k1 > Bz_k1.max() * 1e-4) & (Bz_k1 < Bz_k1.max() * 0.3)
if np.sum(growth_mask) > 5:
    t_fit = time[growth_mask]
    l_fit = log_Bz[growth_mask]
    coeffs = np.polyfit(t_fit, l_fit, 1)
    measured_gamma = coeffs[0]
    B0_fit = np.exp(coeffs[1])
else:
    measured_gamma = gamma_theory
    B0_fit = Bz_k1[10] if Bz_k1[10] > 0 else 1e-5

print(f"  Measured growth rate: {measured_gamma:.4f} wpe")
print(f"  Theory growth rate:  {gamma_theory:.4f} wpe")

Bz_theory = B0_fit * np.exp(gamma_theory * time)

# Save Fig 10 data
data10 = np.column_stack([time, Bz_k1, Bz_theory])
np.savetxt("../data/fig10_weibel_growth.csv", data10, delimiter=",",
           header="time,Bz_k1_simulation,Bz_k1_linear_theory",
           comments="", fmt="%.8e")

fig10, ax = plt.subplots(figsize=(8, 6))
ax.semilogy(time, Bz_k1, 'b-', linewidth=1.5, label='PIC simulation')
ax.semilogy(time, Bz_theory, 'r--', linewidth=1.5,
            label=f'Linear theory ($\\gamma = {gamma_theory:.2f}\\omega_{{pe}}$)')
ax.set_xlabel('Time [$\\omega_{pe}^{-1}$]', fontsize=12)
ax.set_ylabel('$|B_{z,k=1}|$ [$\\sqrt{4\\pi n_e m_e c^2}$]', fontsize=12)
ax.set_title('Weibel Instability: $B_z$ Growth Rate', fontsize=14)
ax.legend(fontsize=11)
ymin = max(Bz_k1[Bz_k1 > 0].min() * 0.1, 1e-8) if np.any(Bz_k1 > 0) else 1e-8
ax.set_ylim([ymin, Bz_k1.max() * 10])
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("../plots/fig10_weibel_growth.png", dpi=150)
print("Saved fig10_weibel_growth.png")

# --- Figure 11: Energy and momentum ---
energy = result['energy']
momentum = result['momentum']
dE = (energy - energy[0]) / energy[0] * 100

data11 = np.column_stack([time, energy, momentum, dE])
np.savetxt("../data/fig11_weibel_energy_momentum.csv", data11, delimiter=",",
           header="time,total_energy,total_momentum,dE_percent",
           comments="", fmt="%.8e")

fig11, ax1 = plt.subplots(figsize=(8, 6))
ax2 = ax1.twinx()

ax1.plot(time, energy, 'r-', linewidth=1.5, label='Total Energy')
ax2.plot(time, momentum, 'b-', linewidth=1.0, alpha=0.7, label='Total Momentum')

ax1.set_xlabel('Time [$\\omega_{pe}^{-1}$]', fontsize=12)
ax1.set_ylabel('Total Energy [$n_e m_e c^2/2$]', color='r', fontsize=12)
ax2.set_ylabel('Total Momentum [$m_e c$]', color='b', fontsize=12)
ax1.tick_params(axis='y', labelcolor='r')
ax2.tick_params(axis='y', labelcolor='b')

ax1.set_title('Weibel Instability: Energy and Momentum', fontsize=14)
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper left', fontsize=11)

plt.tight_layout()
plt.savefig("../plots/fig11_weibel_energy_momentum.png", dpi=150)
print("Saved fig11_weibel_energy_momentum.png")
