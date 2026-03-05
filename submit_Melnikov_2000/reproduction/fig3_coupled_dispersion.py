"""
Figure 3: Coupled dispersion relation eigenvalues for the modified SLAC derivative.

Plots the eigenvalues of the matrix M (Eq. 93) as a function of k
for different values of c = mu/(1-mu), showing the transition between
the small-k (Goldstone + massive) and large-k (uncoupled) regimes.

Reference: Fig. 3 of Melnikov & Weinstein (2000)
"""
import sys
import os
import numpy as np

os.environ["MPLBACKEND"] = "Agg"
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(__file__))
from lattice_schwinger import coupled_dispersion_eigenvalues

# Parameters
e2_over_pi = 1.0  # e^2/pi = 1 (natural units, sets mass scale)
c_values = [2.0, 5.0, 20.0]
npts = 300
k_max = 3.0

data_dir = os.path.join(os.path.dirname(__file__), "..", "data")
plot_dir = os.path.join(os.path.dirname(__file__), "..", "plots")
os.makedirs(data_dir, exist_ok=True)
os.makedirs(plot_dir, exist_ok=True)

k_vals = np.linspace(0.001, k_max, npts)

# Compute and save
header = "k"
columns = [k_vals]
for c in c_values:
    E1, E2 = coupled_dispersion_eigenvalues(k_vals, c, e2_over_pi)
    header += f",E1_c_{c},E2_c_{c}"
    columns.extend([E1, E2])

data = np.column_stack(columns)
np.savetxt(
    os.path.join(data_dir, "fig3_coupled_dispersion.csv"),
    data,
    delimiter=",",
    header=header,
    comments="",
    fmt="%.8f"
)
print(f"Saved fig3 data: {data.shape}")

# Plot
fig, axes = plt.subplots(1, len(c_values), figsize=(5 * len(c_values), 5), sharey=False)
if len(c_values) == 1:
    axes = [axes]

for idx, c in enumerate(c_values):
    ax = axes[idx]
    E1 = columns[2 * idx + 1]
    E2 = columns[2 * idx + 2]

    ax.plot(k_vals, E1, 'b-', linewidth=2, label=r'$E_1$ (lower)')
    ax.plot(k_vals, E2, 'r-', linewidth=2, label=r'$E_2$ (upper)')

    # Show small-k asymptotes
    E1_small = np.sqrt(c) * k_vals
    E2_small = np.sqrt((1 + c) * e2_over_pi) * np.ones_like(k_vals)
    ax.plot(k_vals, E1_small, 'b--', alpha=0.4, label=r'$\sqrt{c}\,|k|$')
    ax.plot(k_vals, E2_small, 'r--', alpha=0.4, label=r'$\sqrt{(1+c)e^2/\pi}$')

    ax.set_xlabel(r"$k$", fontsize=12)
    ax.set_ylabel(r"$E$", fontsize=12)
    ax.set_title(f"$c = {c}$", fontsize=13)
    ax.legend(fontsize=9)
    ax.set_xlim(0, k_max)
    ax.grid(True, alpha=0.3)

fig.suptitle("Coupled dispersion relation (modified SLAC)", fontsize=14, y=1.02)
fig.tight_layout()
fig.savefig(os.path.join(plot_dir, "fig3_coupled_dispersion.png"), dpi=150, bbox_inches='tight')
plt.close(fig)
print("Saved fig3 plot")
