"""
Figure 5: Anomalous commutator W in the continuum limit for different derivatives.

Computes the dimensionless prefactor W/(k/pi) using Eq. (83).
For the continuum Schwinger model, this should equal -1.
Results (Eqs. 84-88):
  naive -> -2  (fermion doubling)
  Wilson (any r) -> -2
  SLAC -> -1 (correct)
  modified SLAC (mu) -> -(1 + mu/(1-mu)) = -1/(1-mu)

Reference: Section III of Melnikov & Weinstein (2000)
"""
import sys
import os
import numpy as np

os.environ["MPLBACKEND"] = "Agg"
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(__file__))
from lattice_schwinger import (
    naive_derivative, wilson_derivative,
    anomalous_commutator_continuum_limit,
    anomalous_commutator_analytical,
)

data_dir = os.path.join(os.path.dirname(__file__), "..", "data")
plot_dir = os.path.join(os.path.dirname(__file__), "..", "plots")
os.makedirs(data_dir, exist_ok=True)
os.makedirs(plot_dir, exist_ok=True)

# Numerical computation for smooth derivatives (naive, Wilson)
print("=== Numerical computations (smooth derivatives) ===")
W_naive = anomalous_commutator_continuum_limit(naive_derivative, npts=50000)
print(f"Naive: W/(k/pi) = {W_naive:.6f} (expected: -2)")

results_num = {"Naive": W_naive}
for r in [0.5, 1.0, 2.0]:
    W = anomalous_commutator_continuum_limit(wilson_derivative, r, npts=50000)
    results_num[f"Wilson r={r}"] = W
    print(f"Wilson r={r}: W/(k/pi) = {W:.6f} (expected: -2)")

# Analytical results for modified SLAC (Eq. 88)
print("\n=== Analytical: Modified SLAC (Eq. 88) ===")
mu_scan = np.linspace(0.05, 0.95, 100)
W_analytical = np.array([
    anomalous_commutator_analytical("modified_slac", mu=mu)
    for mu in mu_scan
])

for mu in [0.3, 0.5, 0.7]:
    W = anomalous_commutator_analytical("modified_slac", mu=mu)
    print(f"mu={mu}: W/(k/pi) = {W:.4f}")

# Save data
data = np.column_stack([mu_scan, W_analytical])
np.savetxt(
    os.path.join(data_dir, "fig5_anomalous_commutator.csv"),
    data,
    delimiter=",",
    header="mu,W_analytical",
    comments="",
    fmt="%.8f"
)

with open(os.path.join(data_dir, "fig5_discrete_results.csv"), "w") as f:
    f.write("derivative,W_numerical,W_analytical\n")
    f.write(f"Naive,{results_num['Naive']:.8f},-2.0\n")
    for r in [0.5, 1.0, 2.0]:
        f.write(f"Wilson_r{r},{results_num[f'Wilson r={r}']:.8f},-2.0\n")
    f.write(f"SLAC,,-1.0\n")

print("Saved fig5 data")

# Plot
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

# Left: W vs mu for modified SLAC (analytical, Eq. 88)
ax1.plot(mu_scan, W_analytical, 'b-', linewidth=2,
         label=r'$W = -(1 + \mu/(1-\mu))$')
ax1.axhline(y=-1, color='green', linestyle=':', linewidth=1.5,
            label='Continuum ($-1$)')
ax1.axhline(y=-2, color='orange', linestyle=':', linewidth=1.5,
            label='Doubled ($-2$)')
ax1.set_xlabel(r"$\mu$", fontsize=14)
ax1.set_ylabel(r"$W / (k/\pi)$", fontsize=14)
ax1.set_title("Anomalous commutator: Modified SLAC (Eq. 88)", fontsize=13)
ax1.legend(fontsize=11)
ax1.grid(True, alpha=0.3)
ax1.set_xlim(0, 1)
ax1.set_ylim(-12, 0)

# Right: Summary bar chart
names = ["Naive", "Wilson\n($r=1$)", "SLAC\n(exact)", "Mod. SLAC\n($\\mu=0.3$)",
         "Mod. SLAC\n($\\mu=0.7$)"]
computed = [
    results_num["Naive"],
    results_num["Wilson r=1.0"],
    -1.0,  # SLAC exact
    anomalous_commutator_analytical("modified_slac", mu=0.3),
    anomalous_commutator_analytical("modified_slac", mu=0.7),
]
expected = [-2.0, -2.0, -1.0,
            -(1 + 0.3/0.7),
            -(1 + 0.7/0.3)]

x = np.arange(len(names))
bars = ax2.bar(x, computed, 0.5, color='steelblue', label='Result')
ax2.axhline(y=-1, color='green', linestyle=':', linewidth=1.5,
            label='Continuum value ($-1$)')
ax2.set_xticks(x)
ax2.set_xticklabels(names, fontsize=10)
ax2.set_ylabel(r"$W / (k/\pi)$", fontsize=14)
ax2.set_title("Anomalous commutator: derivative comparison", fontsize=13)
ax2.legend(fontsize=10)
ax2.grid(True, alpha=0.3, axis='y')

fig.tight_layout()
fig.savefig(os.path.join(plot_dir, "fig5_anomalous_commutator.png"), dpi=150)
plt.close(fig)
print("Saved fig5 plot")
