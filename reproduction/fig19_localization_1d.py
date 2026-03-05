"""
Fig 19: Mean <x> vs system length Lz for 1D Anderson chain.
Parameters: E=1, W=1 (box disorder).
Expected: linear x(Lz) = 1.94 + 0.0282*Lz, giving lambda = 70.92.
Analytical: lambda = 1/Re(gamma) = 72 for E=1, W=1.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from anderson_model import transfer_matrix_1d

W = 1.0
E = 1.0
n_samples = 5000

# System lengths
lengths = np.arange(10, 510, 10)

mean_x = np.zeros(len(lengths))
std_x = np.zeros(len(lengths))

rng = np.random.default_rng(42)

for il, Lz in enumerate(lengths):
    if il % 10 == 0:
        print(f"Length {Lz}/{lengths[-1]}")
    x_vals = np.zeros(n_samples)
    for i in range(n_samples):
        x, g = transfer_matrix_1d(Lz, W, E, disorder="box", seed=rng.integers(0, 2**31))
        x_vals[i] = x
    mean_x[il] = np.mean(x_vals)
    std_x[il] = np.std(x_vals) / np.sqrt(n_samples)

# Linear fit
from numpy.polynomial import polynomial as P
coeffs = np.polyfit(lengths, mean_x, 1)
slope = coeffs[0]
intercept = coeffs[1]
lambda_est = 2.0 / slope  # x = 2*Lz/lambda => slope = 2/lambda

print(f"Linear fit: x = {intercept:.4f} + {slope:.6f} * Lz")
print(f"Estimated localization length: lambda = {lambda_est:.2f}")
print(f"Analytical prediction: lambda = 72")

# Save data
data = np.column_stack([lengths, mean_x, std_x])
np.savetxt("../data/fig19_localization_1d.csv", data,
           delimiter=",", header="Lz,mean_x,std_x", fmt="%.8e", comments="")

# Plot
fig, ax = plt.subplots(figsize=(8, 6))
ax.errorbar(lengths, mean_x, yerr=std_x, fmt='o', markersize=3, label="Numerical")
fit_line = intercept + slope * lengths
ax.plot(lengths, fit_line, 'r-', lw=2,
        label=f"Fit: x = {intercept:.2f} + {slope:.4f}*Lz\n$\\lambda$ = {lambda_est:.1f}")
ax.set_xlabel(r"$L_z$")
ax.set_ylabel(r"$\langle x \rangle$")
ax.set_title(f"1D Anderson: E={E}, W={W}")
ax.legend()
plt.tight_layout()
plt.savefig("../plots/fig19_localization_1d.png", dpi=150)
print("Fig 19 saved.")
