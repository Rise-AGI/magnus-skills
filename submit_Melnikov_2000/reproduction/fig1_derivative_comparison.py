"""
Figure 1: Comparison of one-particle energy spectra for different fermion derivatives.

Compares naive, Wilson (r=1), SLAC, and modified SLAC derivatives.
This supplements the paper's discussion in Section III.
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
    naive_derivative, wilson_derivative, slac_derivative,
    modified_slac_derivative, energy_spectrum
)

npts = 500

data_dir = os.path.join(os.path.dirname(__file__), "..", "data")
plot_dir = os.path.join(os.path.dirname(__file__), "..", "plots")
os.makedirs(data_dir, exist_ok=True)
os.makedirs(plot_dir, exist_ok=True)

xi = np.linspace(0, np.pi, npts)

derivatives = {
    "naive": naive_derivative(xi),
    "wilson_r1": wilson_derivative(xi, r=1.0),
    "slac": slac_derivative(xi),
    "mod_slac_0.5": modified_slac_derivative(xi, mu=0.5),
}

header = "xi"
columns = [xi]
for name, (Z, X) in derivatives.items():
    E = energy_spectrum(Z, X)
    header += f",E_{name}"
    columns.append(E)

data = np.column_stack(columns)
np.savetxt(
    os.path.join(data_dir, "fig1_derivative_comparison.csv"),
    data,
    delimiter=",",
    header=header,
    comments="",
    fmt="%.8f"
)
print(f"Saved fig1 data: {data.shape}")

# Plot
fig, ax = plt.subplots(figsize=(8, 6))
labels = {
    "naive": "Naive",
    "wilson_r1": "Wilson ($r=1$)",
    "slac": "SLAC",
    "mod_slac_0.5": "Modified SLAC ($\\mu=0.5$)",
}
for i, (name, _) in enumerate(derivatives.items()):
    ax.plot(xi / np.pi, columns[i + 1], linewidth=2, label=labels[name])

ax.set_xlabel(r"$ka / \pi$", fontsize=14)
ax.set_ylabel(r"$E_k$ (units of $1/a$)", fontsize=14)
ax.set_title("One-particle energy spectra: derivative comparison", fontsize=14)
ax.legend(fontsize=11)
ax.set_xlim(0, 1)
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(plot_dir, "fig1_derivative_comparison.png"), dpi=150)
plt.close(fig)
print("Saved fig1 plot")
