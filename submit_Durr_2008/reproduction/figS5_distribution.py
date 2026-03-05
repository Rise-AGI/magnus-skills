"""
Figure S5: Systematic error distribution for the nucleon mass.

Reproduces Figure S5 from Durr et al., Science 322, 1224 (2008)
Supplementary Online Material.

Shows the distribution of results from 432 different fitting procedures:
    2 normalization methods x 2 chiral strategies x 3 mass ranges x
    2 continuum extrapolations x 18 time intervals = 432

Each result is weighted by fit quality. The median gives the central
value and the 68% confidence interval gives the systematic error.

Output:
    ../plots/figS5_distribution.png
    ../data/figS5_distribution.csv
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, ".")
from lattice_qcd import (
    systematic_error_analysis, EXPERIMENTAL_MASSES, LATTICE_XI_SET,
)

print("=" * 60)
print("Figure S5: Systematic Error Distribution (Nucleon)")
print("=" * 60)

# ---- Run analysis ----
result = systematic_error_analysis(particle="N", n_bootstrap=500, seed=42)

print(f"Number of procedures: {len(result['results'])}")
print(f"Median: {result['median']:.4f} GeV")
print(f"Statistical error (SEM): {result['stat_err']:.4f} GeV")
print(f"Systematic error: {result['sys_err']:.4f} GeV")
print(f"Total error: {np.sqrt(result['stat_err']**2 + result['sys_err']**2):.4f} GeV")
print(f"Experimental value: {result['experimental']:.3f} GeV")

# ---- Plot ----
fig, ax = plt.subplots(1, 1, figsize=(8, 5))

# Weighted histogram
n_bins = 50
hist_range = (result["median"] - 0.1, result["median"] + 0.1)
bin_edges = np.linspace(hist_range[0], hist_range[1], n_bins + 1)
bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

# Compute weighted histogram
weights_hist = np.zeros(n_bins)
for val, w in zip(result["results"], result["weights"]):
    idx = np.searchsorted(bin_edges, val) - 1
    if 0 <= idx < n_bins:
        weights_hist[idx] += w

# Normalize
weights_hist /= (weights_hist.sum() * (bin_edges[1] - bin_edges[0]))

ax.bar(bin_centers, weights_hist, width=bin_edges[1] - bin_edges[0],
       color="steelblue", alpha=0.7, edgecolor="navy", linewidth=0.5)

# Mark median
ax.axvline(result["median"], color="red", linewidth=2, linestyle="-",
           label=f'Median = {result["median"]:.3f} GeV')

# Mark experimental value
ax.axvline(result["experimental"], color="black", linewidth=2, linestyle="--",
           label=f'Experiment = {result["experimental"]:.3f} GeV')

# Mark 68% CI
lo = result["median"] - result["sys_err"]
hi = result["median"] + result["sys_err"]
ax.axvspan(lo, hi, alpha=0.15, color="red", label=f"68% CI: [{lo:.3f}, {hi:.3f}]")

ax.set_xlabel("$M_N$ [GeV]", fontsize=13)
ax.set_ylabel("Weighted probability density", fontsize=13)
ax.set_title("Distribution of nucleon mass from 432 fitting procedures", fontsize=13)
ax.legend(fontsize=10)
ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig("../plots/figS5_distribution.png", dpi=300)
print("\nSaved: ../plots/figS5_distribution.png")

# ---- Export CSV ----
with open("../data/figS5_distribution.csv", "w") as f:
    f.write("# Figure S5: Nucleon mass distribution from 432 fitting procedures\n")
    f.write("# Columns: bin_center_GeV, weighted_density\n")
    f.write("bin_center_GeV,weighted_density\n")
    for bc, wh in zip(bin_centers, weights_hist):
        f.write(f"{bc:.8f},{wh:.8f}\n")

print("Saved: ../data/figS5_distribution.csv")

# ---- Also compute for all hadrons ----
print("\n" + "=" * 60)
print("Summary: All hadron masses from systematic analysis")
print("=" * 60)
print(f"{'Particle':>10s} {'Median':>8s} {'Stat':>7s} {'Sys':>7s} {'Exp':>8s} {'sigma':>7s}")
print("-" * 55)

all_results = []
for particle in ["rho", "K*", "N", "Lambda", "Sigma", "Delta", "Sigma*", "Xi*", "Omega"]:
    r = systematic_error_analysis(particle=particle, n_bootstrap=200)
    total = np.sqrt(r["stat_err"]**2 + r["sys_err"]**2)
    diff = abs(r["median"] - r["experimental"]) if r["experimental"] else 0
    sigma = diff / total if total > 0 else 0
    print(f"{particle:>10s} {r['median']:8.3f} {r['stat_err']:7.3f} {r['sys_err']:7.3f} "
          f"{r['experimental']:8.3f} {sigma:7.1f}")
    all_results.append(r)

print("Done.")
