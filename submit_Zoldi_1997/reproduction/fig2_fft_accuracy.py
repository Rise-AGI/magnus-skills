"""
Figure 2: Parallel FFT Accuracy - 2D Matrix Decomposition vs Direct FFT

Validates that the 2D matrix FFT decomposition (Eq. 2-4 of the paper)
produces identical results to the standard 1D FFT, up to numerical precision.
"""

import sys
sys.path.insert(0, ".")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from ssf_core import parallel_fft_2d

N_values = [16, 64, 256, 1024, 4096]
errors = []
data_rows = []

print("Comparing parallel FFT (2D decomposition) vs direct FFT...")
for N in N_values:
    M = int(np.sqrt(N))
    assert M * M == N

    # Random complex signal
    np.random.seed(42)
    A = np.random.randn(N) + 1j * np.random.randn(N)

    # Direct FFT
    F_direct = np.fft.fft(A)

    # Parallel FFT via 2D decomposition
    F_parallel = parallel_fft_2d(A, M, M)

    # Compute relative error
    rel_err = np.max(np.abs(F_direct - F_parallel)) / np.max(np.abs(F_direct))
    errors.append(rel_err)
    data_rows.append([N, M, rel_err])
    print(f"  N={N:5d} (M0=M1={M:3d}): max relative error = {rel_err:.2e}")

# Save data
np.savetxt("../data/fig2_fft_accuracy.csv",
           np.array(data_rows),
           delimiter=",",
           header="N,M,max_relative_error",
           fmt=["%.0f", "%.0f", "%.8e"])

# Also do a detailed comparison for N=256
N = 256
M = 16
np.random.seed(42)
A = np.random.randn(N) + 1j * np.random.randn(N)
F_direct = np.fft.fft(A)
F_parallel = parallel_fft_2d(A, M, M)

detail_data = np.column_stack([
    np.arange(N),
    np.abs(F_direct),
    np.abs(F_parallel),
    np.abs(F_direct - F_parallel)
])
np.savetxt("../data/fig2_fft_detail_N256.csv",
           detail_data,
           delimiter=",",
           header="k,abs_F_direct,abs_F_parallel,abs_error",
           fmt="%.8e")

# Plot
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Left: Error vs N
ax1.semilogy(N_values, errors, "bo-", markersize=8)
ax1.set_xlabel("Array Size N")
ax1.set_ylabel("Max Relative Error")
ax1.set_title("2D Decomposition FFT vs Direct FFT")
ax1.grid(True, alpha=0.3)
ax1.set_xscale("log", base=2)

# Right: Detailed comparison for N=256
ax2.plot(np.arange(N), np.abs(F_direct), "b-", label="Direct FFT", alpha=0.7)
ax2.plot(np.arange(N), np.abs(F_parallel), "r--", label="2D Decomposition", alpha=0.7)
ax2.set_xlabel("Frequency index k")
ax2.set_ylabel("|F(k)|")
ax2.set_title(f"Spectrum Comparison (N={N})")
ax2.legend()

plt.tight_layout()
plt.savefig("../plots/fig2_fft_accuracy.png", dpi=150)
print("Saved fig2_fft_accuracy.png")
