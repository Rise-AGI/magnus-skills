"""
Figure 4: One-particle energy spectrum for perfect Wilson derivative.

Plots E_k (in units of 1/a) vs xi = ka.
Reference: Fig. 4 of Melnikov & Weinstein (2000)
"""
import sys
import os
import numpy as np

os.environ["MPLBACKEND"] = "Agg"
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(__file__))
from lattice_schwinger import perfect_wilson_derivative, energy_spectrum

npts = 500

data_dir = os.path.join(os.path.dirname(__file__), "..", "data")
plot_dir = os.path.join(os.path.dirname(__file__), "..", "plots")
os.makedirs(data_dir, exist_ok=True)
os.makedirs(plot_dir, exist_ok=True)

xi = np.linspace(0, np.pi, npts)
Z, X = perfect_wilson_derivative(xi)
E = energy_spectrum(Z, X)

# Save data
data = np.column_stack([xi, Z, X, E])
np.savetxt(
    os.path.join(data_dir, "fig4_perfect_wilson_spectrum.csv"),
    data,
    delimiter=",",
    header="xi,Z,X,E",
    comments="",
    fmt="%.8f"
)
print(f"Saved fig4 data: {data.shape}")

# Plot
fig, ax = plt.subplots(figsize=(8, 6))
ax.plot(xi / np.pi, E, 'b-', linewidth=2, label=r'$E_k = \sqrt{Z_k^2 + X_k^2}$')
ax.plot(xi / np.pi, np.abs(Z), 'g--', linewidth=1.5, alpha=0.6, label=r'$|Z_k|$')
ax.plot(xi / np.pi, np.abs(X), 'r--', linewidth=1.5, alpha=0.6, label=r'$|X_k|$')

ax.set_xlabel(r"$ka / \pi$", fontsize=14)
ax.set_ylabel(r"$E_k$ (units of $1/a$)", fontsize=14)
ax.set_title("One-particle energy: Perfect Wilson derivative", fontsize=14)
ax.legend(fontsize=12)
ax.set_xlim(0, 1)
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(plot_dir, "fig4_perfect_wilson_spectrum.png"), dpi=150)
plt.close(fig)
print("Saved fig4 plot")
