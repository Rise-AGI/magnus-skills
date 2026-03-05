"""
Fig 20: Probability distribution p(x) for the 1D disordered chain.
Parameters: Lz=300, W=1, E=1.
Expected: nearly Gaussian with <x>=10.41, var(x)=14.45.
Statistical ensemble: N_stat = 10^5 (reduced from paper's 10^9 for feasibility).
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from anderson_model import transfer_matrix_1d

Lz = 300
W = 1.0
E = 1.0
n_samples = 200000

rng = np.random.default_rng(42)
x_vals = np.zeros(n_samples)

for i in range(n_samples):
    if i % 50000 == 0:
        print(f"Sample {i}/{n_samples}")
    x, g = transfer_matrix_1d(Lz, W, E, disorder="box", seed=rng.integers(0, 2**31))
    x_vals[i] = x

mean_x = np.mean(x_vals)
var_x = np.var(x_vals)
print(f"<x> = {mean_x:.4f}, var(x) = {var_x:.4f}")
print(f"Paper values: <x> = 10.41, var(x) = 14.45")

# Histogram
bins = np.linspace(0, max(30, mean_x + 4*np.sqrt(var_x)), 150)
hist, edges = np.histogram(x_vals, bins=bins, density=True)
centers = 0.5 * (edges[:-1] + edges[1:])

# Gaussian comparison
gaussian = np.exp(-(centers - mean_x)**2 / (2 * var_x)) / np.sqrt(2 * np.pi * var_x)

# Save data
data = np.column_stack([centers, hist])
np.savetxt("../data/fig20_px_distribution.csv", data,
           delimiter=",", header="x,p_x", fmt="%.8e", comments="")
np.savetxt("../data/fig20_px_gaussian.csv",
           np.column_stack([centers, gaussian]),
           delimiter=",", header="x,gaussian", fmt="%.8e", comments="")

# Plot
fig, ax = plt.subplots(figsize=(8, 6))
ax.plot(centers, hist, 'b-', lw=1.5, label=f"p(x), Lz={Lz}")
ax.plot(centers, gaussian, 'r--', lw=1.5,
        label=f"Gaussian: <x>={mean_x:.2f}, var={var_x:.2f}")
ax.set_xlabel("x")
ax.set_ylabel("p(x)")
ax.set_title(f"1D Anderson: Lz={Lz}, W={W}, E={E}")
ax.legend()
ax.set_xlim(0, 30)
plt.tight_layout()
plt.savefig("../plots/fig20_px_distribution.png", dpi=150)
print("Fig 20 saved.")
