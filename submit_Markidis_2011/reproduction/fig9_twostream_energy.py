"""
Figure 9: Energy conservation comparison for the two-stream instability.

Compares total energy history between:
- Energy-conserving PIC (red line, left y axis): ~1e-4% variation
- Explicit momentum-conserving PIC (blue line, right y axis): ~5% variation

Parameters (same as Figure 8):
- Drift velocity: +/- 0.2c
- Box length: 2.053 c/omega_pe
- 64 grid points, reduced particles for speed
- dt = 0.1 omega_pe^{-1}
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys, os

sys.path.insert(0, os.path.dirname(os.path.abspath(".")))
from ecpic_es import run_electrostatic_pic, run_explicit_pic

# Parameters
WP = 1.0
c = 1.0
V0 = 0.2 * c
L = 2.053 * c / WP
NG = 64
N = 3000
DT = 0.1 / WP
NT = 300
QM = -1.0
VT = 0.0

np.random.seed(42)
x0 = np.linspace(0, L - L / N, N)
v0 = VT * np.random.randn(N)
pm = np.ones(N)
pm[1::2] = -1.0
v0 = v0 + pm * V0

XP1 = 1.0
mode = 1
x0 = x0 + XP1 * (L / N) * np.sin(2.0 * np.pi * x0 / L * mode)
x0 = x0 % L

print("Running EC-PIC for two-stream instability (energy comparison)...")
result_ec = run_electrostatic_pic(L, DT, NT, NG, N, WP, QM,
                                   x0.copy(), v0.copy(), tol=1e-7)

print("Running explicit PIC for two-stream instability...")
result_ex = run_explicit_pic(L, DT, NT, NG, N, WP, QM,
                              x0.copy(), v0.copy())

time = np.arange(NT + 1) * DT

# Energy relative variation
E_ec = result_ec['energy']
E_ex = result_ex['energy']
dE_ec = (E_ec - E_ec[0]) / E_ec[0] * 100  # in percent
dE_ex = (E_ex - E_ex[0]) / E_ex[0] * 100

print(f"  EC-PIC max energy variation: {np.max(np.abs(dE_ec)):.2e}%")
print(f"  Explicit max energy variation: {np.max(np.abs(dE_ex)):.2e}%")

# Save data
data = np.column_stack([time, E_ec, E_ex, dE_ec, dE_ex])
header = "time,energy_ec,energy_explicit,dE_ec_percent,dE_explicit_percent"
np.savetxt("../data/fig9_twostream_energy.csv", data, delimiter=",",
           header=header, comments="", fmt="%.8e")

# Plot (dual y-axis like in the paper)
fig, ax1 = plt.subplots(figsize=(8, 6))
ax2 = ax1.twinx()

ax1.plot(time, E_ec, 'r-', linewidth=1.5, label='EC-PIC')
ax2.plot(time, E_ex, 'b-', linewidth=1.5, label='Explicit PIC')

ax1.set_xlabel('Time [$\\omega_{pe}^{-1}$]', fontsize=12)
ax1.set_ylabel('Total Energy (EC-PIC) [$n_e m_e c^2/2$]', color='r', fontsize=12)
ax2.set_ylabel('Total Energy (Explicit) [$n_e m_e c^2/2$]', color='b', fontsize=12)
ax1.tick_params(axis='y', labelcolor='r')
ax2.tick_params(axis='y', labelcolor='b')

ax1.set_title('Two-Stream Instability: Energy Conservation', fontsize=14)

lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper left', fontsize=11)

plt.tight_layout()
plt.savefig("../plots/fig9_twostream_energy.png", dpi=150)
print("Saved fig9_twostream_energy.png and fig9_twostream_energy.csv")
