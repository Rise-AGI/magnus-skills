"""
Figure 5: Energy histories with different solver tolerance values.

Shows that smaller tolerance -> better energy conservation in the EC-PIC.
Tests with tolerances: 1e-5, 1e-7, 1e-9

Parameters: two-stream instability setup from paper.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys, os

sys.path.insert(0, os.path.dirname(os.path.abspath(".")))
from ecpic_es import run_electrostatic_pic

# Parameters
WP = 1.0
c = 1.0
V0 = 0.2 * c
L = 2.053 * c / WP
NG = 64
N = 2000
DT = 0.1 / WP
NT = 150
QM = -1.0

np.random.seed(42)
x0 = np.linspace(0, L - L / N, N)
v0 = np.zeros(N)
pm = np.ones(N)
pm[1::2] = -1.0
v0 = v0 + pm * V0
XP1 = 1.0
mode = 1
x0 = x0 + XP1 * (L / N) * np.sin(2.0 * np.pi * x0 / L * mode)
x0 = x0 % L

tolerances = [1e-5, 1e-7, 1e-9]
results = {}

for tol in tolerances:
    print(f"Running EC-PIC with tolerance = {tol:.0e}...")
    result = run_electrostatic_pic(L, DT, NT, NG, N, WP, QM,
                                    x0.copy(), v0.copy(), tol=tol)
    results[tol] = result

time = np.arange(NT + 1) * DT

# Save data
cols = [time]
headers = ["time"]
for tol in tolerances:
    E = results[tol]['energy']
    dE = (E - E[0]) / E[0] * 100
    cols.append(dE)
    headers.append(f"dE_tol_{tol:.0e}_percent")

data = np.column_stack(cols)
np.savetxt("../data/fig5_tolerance_energy.csv", data, delimiter=",",
           header=",".join(headers), comments="", fmt="%.8e")

# Plot
fig, ax = plt.subplots(figsize=(8, 6))
colors = ['b', 'r', 'g']
for i, tol in enumerate(tolerances):
    E = results[tol]['energy']
    dE = (E - E[0]) / E[0] * 100
    ax.plot(time, dE, color=colors[i], linewidth=1.5,
            label=f'$\\epsilon_a = \\epsilon_r = {tol:.0e}$')

ax.set_xlabel('Time [$\\omega_{pe}^{-1}$]', fontsize=12)
ax.set_ylabel('$\\Delta E / E_0$ [%]', fontsize=12)
ax.set_title('Energy Conservation vs Solver Tolerance\\n(Two-Stream Instability)', fontsize=14)
ax.legend(fontsize=11)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("../plots/fig5_tolerance_energy.png", dpi=150)
print("Saved fig5_tolerance_energy.png and fig5_tolerance_energy.csv")
