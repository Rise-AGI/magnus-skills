"""
Figure 5 (paper's Figure 7): Filtering of the transition correlator ln(C21).
Demonstrates the inverse-variance weighted filtering (Eqs. 6-7) applied to
synthetic ARA data for ln(C21) at T=5.
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(".")))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Generate synthetic ln(C21) data at T=5
# The signal is a smooth, roughly linear function of R
np.random.seed(123)

R_values = np.linspace(8, 20, 175)  # dense ARA R values

# True signal: roughly linear decrease
B_true = -0.5 * R_values - 2.0

# Errors grow with R (larger separations = more noise)
sigma = 0.05 * np.exp(0.15 * (R_values - 8))

# Noisy measurements
B_noisy = B_true + np.random.normal(0, 1, len(R_values)) * sigma

# Apply inverse-variance weighted filter (Eqs. 6-7)
R_f = 0.5  # filter radius from the paper

B_filtered = np.zeros_like(B_noisy)
for i in range(len(R_values)):
    mask = np.abs(R_values - R_values[i]) <= R_f
    if np.sum(mask) > 0:
        weights = 1.0 / sigma[mask]**2
        B_filtered[i] = np.sum(weights * B_noisy[mask]) / np.sum(weights)
    else:
        B_filtered[i] = B_noisy[i]

# Save data
data_dir = os.path.join(os.path.dirname(os.path.abspath(".")), "data")
os.makedirs(data_dir, exist_ok=True)

with open(os.path.join(data_dir, "fig5_filtering.csv"), "w") as f:
    f.write("R,B_true,B_noisy,sigma,B_filtered\n")
    for i in range(len(R_values)):
        f.write(f"{R_values[i]:.8f},{B_true[i]:.8f},{B_noisy[i]:.8f},{sigma[i]:.8f},{B_filtered[i]:.8f}\n")

# Plot (mimicking Figure 7 layout)
fig, ax = plt.subplots(1, 1, figsize=(8, 6))

# Unfiltered (upper, with error bars) - subsample for clarity
step = 5
ax.errorbar(R_values[::step], B_noisy[::step], yerr=sigma[::step],
            fmt='o', markersize=4, color='blue', alpha=0.6, capsize=2,
            label='ARA + smearing (unfiltered)')

# Filtered (shifted down by 8, as in paper)
shift = -8
ax.plot(R_values[::step], B_filtered[::step] + shift, 's', markersize=4,
        color='red', label=f'Filtered (shifted by {shift})')

# True signal for reference
ax.plot(R_values, B_true, 'k--', linewidth=1, alpha=0.5, label='True signal')
ax.plot(R_values, B_true + shift, 'k--', linewidth=1, alpha=0.5)

ax.axvline(13, color='gray', linestyle=':', alpha=0.5)
ax.annotate('String breaking', (13, B_true[len(B_true)//2]),
            textcoords="offset points", xytext=(5, 10), fontsize=9)

ax.set_xlabel('$R$ (lattice units)', fontsize=14)
ax.set_ylabel('$\\ln C_{21}(R, T=5)$', fontsize=14)
ax.set_title('Filtering of transition correlator (Bolder et al. 2000)', fontsize=13)
ax.legend(fontsize=10, loc='lower left')
ax.set_xlim(7.5, 20.5)
ax.grid(True, alpha=0.3)

plots_dir = os.path.join(os.path.dirname(os.path.abspath(".")), "plots")
os.makedirs(plots_dir, exist_ok=True)
plt.tight_layout()
plt.savefig(os.path.join(plots_dir, "fig5_filtering.png"), dpi=150)
print("Figure 5 saved.")
