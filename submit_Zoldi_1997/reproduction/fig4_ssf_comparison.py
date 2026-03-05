"""
Figure 4: SSF Method Comparison - Serial vs Parallel Implementation

Compares the serial SSF (using standard 1D FFT) with the parallel SSF
(using 2D matrix decomposed FFT) for soliton propagation.
Verifies both methods produce identical physical results.
Also times the 2D-decomposition approach vs standard FFT.
"""

import sys
sys.path.insert(0, ".")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from ssf_core import run_ssf_simulation, soliton_initial, time_fft_methods

# Part 1: Compare serial vs parallel SSF
N = 256  # Must be perfect square
L = 20.0
T_total = 5.0
n_steps = 2000

print("Running serial SSF...")
x_s, t_s, A_s = run_ssf_simulation(
    N=N, L=L, T_total=T_total, n_steps=n_steps,
    initial_func=lambda xx: soliton_initial(xx, amplitude=1.0),
    use_parallel=False
)

print("Running parallel SSF (2D decomposed FFT)...")
x_p, t_p, A_p = run_ssf_simulation(
    N=N, L=L, T_total=T_total, n_steps=n_steps,
    initial_func=lambda xx: soliton_initial(xx, amplitude=1.0),
    use_parallel=True
)

# Compute differences
max_diff = np.max(np.abs(A_s - A_p))
print(f"Max difference between serial and parallel: {max_diff:.2e}")

# Save comparison data
comparison_rows = []
for i in range(len(t_s)):
    rms_diff = np.sqrt(np.mean((A_s[i] - A_p[i])**2))
    max_d = np.max(np.abs(A_s[i] - A_p[i]))
    comparison_rows.append([t_s[i], rms_diff, max_d])

np.savetxt("../data/fig4_ssf_comparison.csv",
           np.array(comparison_rows),
           delimiter=",",
           header="t,rms_difference,max_difference",
           fmt="%.8e")

# Part 2: FFT timing comparison
print("\nTiming FFT methods...")
N_timing = [16, 64, 256, 1024, 4096, 16384]
timing_results = time_fft_methods(N_timing, n_repeats=20)

timing_rows = []
for N_t, (t_serial, t_parallel) in sorted(timing_results.items()):
    ratio = t_parallel / t_serial
    print(f"  N={N_t:6d}: serial={t_serial:.6f}s, parallel={t_parallel:.6f}s, ratio={ratio:.2f}")
    timing_rows.append([N_t, t_serial, t_parallel, ratio])

np.savetxt("../data/fig4_fft_timing.csv",
           np.array(timing_rows),
           delimiter=",",
           header="N,serial_time_sec,parallel_2d_time_sec,parallel_serial_ratio",
           fmt="%.8e")

# Plot
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Top left: Serial final profile
ax = axes[0, 0]
ax.plot(x_s, A_s[0], "b-", label="t=0", alpha=0.5)
ax.plot(x_s, A_s[-1], "r-", label=f"t={t_s[-1]:.1f}")
ax.set_xlabel("x")
ax.set_ylabel("|A(x,t)|")
ax.set_title("Serial SSF")
ax.legend()
ax.set_xlim(-10, 10)

# Top right: Parallel final profile
ax = axes[0, 1]
ax.plot(x_p, A_p[0], "b-", label="t=0", alpha=0.5)
ax.plot(x_p, A_p[-1], "r-", label=f"t={t_p[-1]:.1f}")
ax.set_xlabel("x")
ax.set_ylabel("|A(x,t)|")
ax.set_title("Parallel SSF (2D Decomposition)")
ax.legend()
ax.set_xlim(-10, 10)

# Bottom left: Difference over time
ax = axes[1, 0]
times = [row[0] for row in comparison_rows]
max_diffs = [row[2] for row in comparison_rows]
ax.semilogy(times, max_diffs, "k-")
ax.set_xlabel("Time")
ax.set_ylabel("Max |Serial - Parallel|")
ax.set_title("Serial vs Parallel Difference")
ax.grid(True, alpha=0.3)

# Bottom right: FFT timing
ax = axes[1, 1]
Ns = [row[0] for row in timing_rows]
t_ser = [row[1] for row in timing_rows]
t_par = [row[2] for row in timing_rows]
ax.loglog(Ns, t_ser, "bo-", label="Standard 1D FFT")
ax.loglog(Ns, t_par, "rs-", label="2D Decomposition FFT")
ax.set_xlabel("Array Size N")
ax.set_ylabel("Time (seconds)")
ax.set_title("FFT Timing Comparison")
ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("../plots/fig4_ssf_comparison.png", dpi=150)
print("Saved fig4_ssf_comparison.png")
