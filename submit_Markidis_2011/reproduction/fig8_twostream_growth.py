"""
Figure 8: Two-stream instability - comparison with linear theory.

The k=1 spectral component of E grows exponentially with growth rate
gamma = 0.35355 * omega_pe as predicted by linear theory.

Parameters from the paper:
- Drift velocity: +/- 0.2c
- Box length: 2.053 c/omega_pe
- 64 grid points, 200000 particles
- dt = 0.1 omega_pe^{-1}
- Tolerance: 1e-8

We use reduced particle count for computational feasibility.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(".")))
from ecpic_es import run_electrostatic_pic

# Parameters (from paper Section 7.2)
WP = 1.0          # plasma frequency
c = 1.0           # speed of light (normalized)
V0 = 0.2 * c      # beam drift velocity
L = 2.053 * c / WP  # box length
NG = 64            # grid cells
N = 3000          # particles (reduced from 200000 for speed)
DT = 0.1 / WP     # time step
NT = 300           # time steps
QM = -1.0          # charge-to-mass ratio for electrons
VT = 0.0           # no thermal spread (cold beams)

# Seed for reproducibility
np.random.seed(42)

# Initial positions: uniform
x0 = np.linspace(0, L - L / N, N)

# Initial velocities: two counter-streaming beams
v0 = VT * np.random.randn(N)
pm = np.ones(N)
pm[1::2] = -1.0
v0 = v0 + pm * V0

# Perturbation (as in appendix code)
XP1 = 1.0
V1 = 0.0
mode = 1
v0 = v0 + V1 * np.sin(2.0 * np.pi * x0 / L * mode)
x0 = x0 + XP1 * (L / N) * np.sin(2.0 * np.pi * x0 / L * mode)
x0 = x0 % L

print("Running energy-conserving PIC for two-stream instability...")
print(f"  N={N}, NG={NG}, NT={NT}, DT={DT:.3f}, L={L:.4f}")

# Save field at every step to extract spectral components
save_steps = list(range(NT))
result = run_electrostatic_pic(L, DT, NT, NG, N, WP, QM, x0.copy(), v0.copy(),
                                tol=1e-7, save_field_at=save_steps)

# Extract k=1 spectral component of E field
dx = L / NG
xg = np.arange(NG) * dx
k1 = 2.0 * np.pi / L  # k = 1 * omega_pe / c

E_k1 = np.zeros(NT)
for it in range(NT):
    if it in result['field_data']:
        E = result['field_data'][it]
        # Fourier component at k=1
        E_hat = np.fft.fft(E)
        E_k1[it] = np.abs(E_hat[1]) * 2.0 / NG

# Time array
time = np.arange(NT) * DT

# Linear theory: growth rate for k=1 mode
gamma_theory = 0.35355 * WP

# Normalize to find initial amplitude
# Find a good fit region in the linear growth phase
log_Ek1 = np.log(E_k1 + 1e-30)
# Find where linear growth starts and ends
growth_start = 20
growth_end = min(NT - 1, 150)
# Fit the growth region
valid = (E_k1[growth_start:growth_end] > 0)
if np.sum(valid) > 5:
    t_fit = time[growth_start:growth_end][valid]
    log_fit = log_Ek1[growth_start:growth_end][valid]
    from numpy.polynomial import polynomial as P
    coeffs = np.polyfit(t_fit, log_fit, 1)
    measured_gamma = coeffs[0]
    print(f"  Measured growth rate: {measured_gamma:.4f}")
    print(f"  Theory growth rate:  {gamma_theory:.4f}")
    # Linear theory line
    E_theory = np.exp(coeffs[1]) * np.exp(gamma_theory * time)
else:
    measured_gamma = gamma_theory
    E_theory = E_k1[growth_start] * np.exp(gamma_theory * (time - time[growth_start]))

# Save data
data = np.column_stack([time, E_k1, E_theory[:NT]])
header = "time,E_k1_simulation,E_k1_linear_theory"
np.savetxt("../data/fig8_twostream_growth.csv", data, delimiter=",",
           header=header, comments="", fmt="%.8e")

# Plot
fig, ax = plt.subplots(figsize=(8, 6))
ax.semilogy(time, E_k1, 'b-', linewidth=1.5, label='EC-PIC simulation')
ax.semilogy(time, E_theory[:NT], 'r--', linewidth=1.5,
            label=f'Linear theory ($\\gamma = {gamma_theory:.4f}\\omega_{{pe}}$)')
ax.set_xlabel('Time [$\\omega_{pe}^{-1}$]', fontsize=12)
ax.set_ylabel('$|E_{k=1}|$ [$\\sqrt{4\\pi n_e m_e c^2}$]', fontsize=12)
ax.set_title('Two-Stream Instability: Growth Rate Comparison', fontsize=14)
ax.legend(fontsize=11)
ax.set_xlim([0, time[-1]])
ymin = max(E_k1[E_k1 > 0].min() * 0.1, 1e-8)
ax.set_ylim([ymin, E_k1.max() * 10])
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("../plots/fig8_twostream_growth.png", dpi=150)
print("Saved fig8_twostream_growth.png and fig8_twostream_growth.csv")
