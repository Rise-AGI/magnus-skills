"""
Figure 2: One-particle energy spectrum for modified SLAC derivative.

Plots E_k (in units of 1/a) vs xi = ka for different values of mu.
Reference: Fig. 2 of Melnikov & Weinstein (2000)
"""
import sys
import os
import numpy as np

MPLBACKEND = os.environ.get("MPLBACKEND", "Agg")
os.environ["MPLBACKEND"] = "Agg"
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(os.getcwd())))
sys.path.insert(0, os.path.dirname(__file__))
from lattice_schwinger import modified_slac_derivative, energy_spectrum

# Parameters matching the paper
mu_values = [0.3, 0.7]
npts = 500

# Output directories
data_dir = os.path.join(os.path.dirname(__file__), "..", "data")
plot_dir = os.path.join(os.path.dirname(__file__), "..", "plots")
os.makedirs(data_dir, exist_ok=True)
os.makedirs(plot_dir, exist_ok=True)

xi = np.linspace(0, np.pi, npts)

# Compute and save data
header = "xi"
columns = [xi]
for mu in mu_values:
    Z, X = modified_slac_derivative(xi, mu)
    E = energy_spectrum(Z, X)
    header += f",E_mu_{mu}"
    columns.append(E)

data = np.column_stack(columns)
np.savetxt(
    os.path.join(data_dir, "fig2_modified_slac_spectrum.csv"),
    data,
    delimiter=",",
    header=header,
    comments="",
    fmt="%.8f"
)
print(f"Saved fig2 data: {data.shape}")

# Plot
fig, ax = plt.subplots(figsize=(8, 6))
for i, mu in enumerate(mu_values):
    ax.plot(xi / np.pi, columns[i + 1], label=f"$\\mu = {mu}$", linewidth=2)

ax.set_xlabel(r"$ka / \pi$", fontsize=14)
ax.set_ylabel(r"$E_k$ (units of $1/a$)", fontsize=14)
ax.set_title("One-particle energy: Modified SLAC derivative", fontsize=14)
ax.legend(fontsize=12)
ax.set_xlim(0, 1)
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(plot_dir, "fig2_modified_slac_spectrum.png"), dpi=150)
plt.close(fig)
print("Saved fig2 plot")
