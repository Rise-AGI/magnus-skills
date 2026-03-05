"""
Figure 7: Total energy history for Maxwellian plasma (finite grid instability test).

Compares energy conservation between:
- Energy-conserving PIC (red): ~1e-7% variation
- Explicit momentum-conserving PIC (blue): ~2% increase (finite grid instability)

Parameters from paper Section 7.1:
- Thermal velocity: v_the = 0.2c
- Box length: 50*pi*c/omega_pe
- 64 grid cells, 50000 particles (reduced here)
- dt = 0.5 omega_pe^{-1}
- 200 cycles
- Debye length ~ 10x smaller than grid spacing -> finite grid instability in explicit PIC
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys, os

sys.path.insert(0, os.path.dirname(os.path.abspath(".")))
from ecpic_es import run_electrostatic_pic, run_explicit_pic

# Parameters (from paper Section 7.1)
WP = 1.0
c = 1.0
VTH = 0.2 * c     # thermal velocity
L = 50.0 * np.pi * c / WP  # box length
NG = 64
N = 2000          # reduced from 50000 for speed
DT = 0.5 / WP     # time step
NT = 200           # cycles
QM = -1.0

dx = L / NG
print(f"Grid spacing: {dx:.4f} c/wpe")
print(f"Debye length: {VTH/WP:.4f} c/wpe")
print(f"Ratio dx/lambda_D: {dx/(VTH/WP):.1f}")

np.random.seed(42)

# Maxwellian plasma: uniform positions, Gaussian velocities
x0 = np.random.uniform(0, L, N)
v0 = VTH * np.random.randn(N)

print("Running EC-PIC for Maxwellian plasma (finite grid test)...")
result_ec = run_electrostatic_pic(L, DT, NT, NG, N, WP, QM,
                                   x0.copy(), v0.copy(), tol=1e-7)

print("Running explicit PIC for Maxwellian plasma...")
result_ex = run_explicit_pic(L, DT, NT, NG, N, WP, QM,
                              x0.copy(), v0.copy())

time = np.arange(NT + 1) * DT

E_ec = result_ec['energy']
E_ex = result_ex['energy']
dE_ec = (E_ec - E_ec[0]) / E_ec[0] * 100
dE_ex = (E_ex - E_ex[0]) / E_ex[0] * 100

print(f"  EC-PIC max energy variation: {np.max(np.abs(dE_ec)):.2e}%")
print(f"  Explicit max energy variation: {np.max(np.abs(dE_ex)):.2e}%")

# Save data
data = np.column_stack([time, E_ec, E_ex, dE_ec, dE_ex])
header = "time,energy_ec,energy_explicit,dE_ec_percent,dE_explicit_percent"
np.savetxt("../data/fig7_finite_grid_energy.csv", data, delimiter=",",
           header=header, comments="", fmt="%.8e")

# Plot (dual y-axis)
fig, ax1 = plt.subplots(figsize=(8, 6))
ax2 = ax1.twinx()

ax1.plot(time, E_ec, 'r-', linewidth=1.5, label='EC-PIC')
ax2.plot(time, E_ex, 'b-', linewidth=1.5, label='Explicit PIC')

ax1.set_xlabel('Time [$\\omega_{pe}^{-1}$]', fontsize=12)
ax1.set_ylabel('Total Energy (EC-PIC) [$n_e m_e c^2/2$]', color='r', fontsize=12)
ax2.set_ylabel('Total Energy (Explicit) [$n_e m_e c^2/2$]', color='b', fontsize=12)
ax1.tick_params(axis='y', labelcolor='r')
ax2.tick_params(axis='y', labelcolor='b')

ax1.set_title('Maxwellian Plasma: Energy Conservation\\n(Finite Grid Instability Test)', fontsize=14)

lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper left', fontsize=11)

plt.tight_layout()
plt.savefig("../plots/fig7_finite_grid_energy.png", dpi=150)
print("Saved fig7_finite_grid_energy.png and fig7_finite_grid_energy.csv")
