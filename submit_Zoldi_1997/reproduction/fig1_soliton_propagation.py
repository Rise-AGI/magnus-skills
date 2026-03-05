"""
Figure 1: Soliton Propagation via Split-Step Fourier Method

Demonstrates that the SSF method correctly preserves the bright soliton
solution of NLSE: A(x,t) = sech(x) * exp(i*t/2).
The soliton should maintain its shape and amplitude over propagation.
"""

import sys
sys.path.insert(0, ".")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from ssf_core import run_ssf_simulation, soliton_initial

# Parameters
N = 1024
L = 20.0
T_total = 10.0
n_steps = 5000
sigma = 1.0
d = 1.0

print("Running SSF soliton propagation (N=1024, T=10)...")
x, t_array, A_history = run_ssf_simulation(
    N=N, L=L, T_total=T_total, n_steps=n_steps,
    sigma=sigma, d=d,
    initial_func=lambda xx: soliton_initial(xx, amplitude=1.0, velocity=0.0)
)

# Save data
data = np.column_stack([t_array] + [A_history[:, i] for i in range(0, N, N // 200)])
header = "t," + ",".join([f"x={x[i]:.4f}" for i in range(0, N, N // 200)])
np.savetxt("../data/fig1_soliton_propagation.csv", data, delimiter=",",
           header=header, fmt="%.8e")

# Also save the 2D snapshot data for contour plot
# Save selected spatial points at all times
x_indices = np.arange(0, N, 1)
with open("../data/fig1_soliton_2d.csv", "w") as f:
    f.write("# Rows: time snapshots, Columns: spatial positions\n")
    f.write("# x values: " + ",".join([f"{x[i]:.6f}" for i in x_indices]) + "\n")
    f.write("# t values: " + ",".join([f"{t:.6f}" for t in t_array]) + "\n")
    for row in A_history:
        f.write(",".join([f"{row[i]:.8e}" for i in x_indices]) + "\n")

# Plot
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Left: Waterfall of |A(x,t)| at selected times
plot_times = [0, len(t_array)//4, len(t_array)//2, 3*len(t_array)//4, -1]
for idx in plot_times:
    ax1.plot(x, A_history[idx], label=f"t={t_array[idx]:.1f}")
ax1.set_xlabel("x")
ax1.set_ylabel("|A(x,t)|")
ax1.set_title("Soliton Profile at Different Times")
ax1.legend()
ax1.set_xlim(-10, 10)

# Right: Contour plot
T_mesh, X_mesh = np.meshgrid(x, t_array)
ax2.contourf(T_mesh, X_mesh, A_history, levels=50, cmap="inferno")
ax2.set_xlabel("x")
ax2.set_ylabel("t")
ax2.set_title("|A(x,t)| Evolution")
ax2.set_xlim(-10, 10)

plt.tight_layout()
plt.savefig("../plots/fig1_soliton_propagation.png", dpi=150)
print("Saved fig1_soliton_propagation.png")

# Verify soliton preservation
initial_peak = np.max(A_history[0])
final_peak = np.max(A_history[-1])
print(f"Initial peak: {initial_peak:.6f}")
print(f"Final peak:   {final_peak:.6f}")
print(f"Relative error: {abs(final_peak - initial_peak) / initial_peak:.2e}")
