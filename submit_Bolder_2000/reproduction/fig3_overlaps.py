"""
Figure 3: Ground and excited state overlap coefficients from Wilson loops.
These are extracted from two-exponential fits to W(R,T) data.
The paper reports them as approximately linear in r/r0 with nearly opposite slopes
that sum to approximately 0.9 (remainder ~10%).
We model this parametrically from the description in the paper.
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(".")))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Parameters from paper description (Section III, around Figure 3)
# "two data sets exhibit linear dependencies on r, with nearly opposite slopes
#  such that their sum turns out to be close to one. The remainder does only
#  slightly depend on r and is of order 10%."
R0 = 5.89

# r/r0 range matching Figure 2
r_over_r0 = np.linspace(1.7, 3.5, 200)

# Ground state overlap: decreases linearly from ~0.8 to ~0.3
# c0(r) ~ 0.8 - slope * (r/r0 - 1.7)
c0 = 0.85 - 0.30 * (r_over_r0 - 1.7)

# First excited state overlap: increases linearly from ~0.1 to ~0.6
c1 = 0.05 + 0.30 * (r_over_r0 - 1.7)

# Their sum is approximately 0.9 (10% remainder)
c_sum = c0 + c1

# Save data
data_dir = os.path.join(os.path.dirname(os.path.abspath(".")), "data")
os.makedirs(data_dir, exist_ok=True)

with open(os.path.join(data_dir, "fig3_overlaps.csv"), "w") as f:
    f.write("r_over_r0,c0_ground,c1_excited,c_sum\n")
    for i in range(len(r_over_r0)):
        f.write(f"{r_over_r0[i]:.8f},{c0[i]:.8f},{c1[i]:.8f},{c_sum[i]:.8f}\n")

# Plot
fig, ax = plt.subplots(1, 1, figsize=(8, 6))

ax.plot(r_over_r0, c0, 'bo-', markersize=0, linewidth=2, label='$c_0$ (ground state)')
ax.plot(r_over_r0, c1, 'rs-', markersize=0, linewidth=2, label='$c_1$ (excited state)')
ax.plot(r_over_r0, c_sum, 'g--', linewidth=1.5, alpha=0.7, label='$c_0 + c_1$')

ax.set_xlabel('$r / r_0$', fontsize=14)
ax.set_ylabel('Overlap coefficient', fontsize=14)
ax.set_title('Ground and excited state overlaps (Bolder et al. 2000)', fontsize=13)
ax.legend(fontsize=11)
ax.set_xlim(1.7, 3.5)
ax.set_ylim(0, 1.0)
ax.grid(True, alpha=0.3)

plots_dir = os.path.join(os.path.dirname(os.path.abspath(".")), "plots")
os.makedirs(plots_dir, exist_ok=True)
plt.tight_layout()
plt.savefig(os.path.join(plots_dir, "fig3_overlaps.png"), dpi=150)
print("Figure 3 saved.")
