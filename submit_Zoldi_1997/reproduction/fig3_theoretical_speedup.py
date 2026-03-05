"""
Figure 3: Speedup Analysis - Empirical Data and Theoretical Model

Plots the empirical speedup data from Tables 1 and 2 of the paper,
alongside the theoretical speedup formula (Eq. 6):
    SU = 2P / (1 + xi/K + f * 2^K / (P*K))

The theoretical model captures qualitative trends but uses fitted parameters
for each (P, paradigm) combination to match the empirical data.
"""

import sys
sys.path.insert(0, ".")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from ssf_core import theoretical_speedup

K_values = np.array([6, 7, 8, 9])  # N = 2^12 to 2^18

# Paper's empirical data (Tables 1 and 2)
# Shared memory (Table 1)
sm_data = {
    6: {"T1": 49.5, "T2": 29.5, "T4": 19.5},   # N=2^12, T=8000
    7: {"T1": 51.5, "T2": 30.5, "T4": 18.5},   # N=2^14, T=2000
    8: {"T1": 65.5, "T2": 33.5, "T4": 19.5},   # N=2^16, T=500
    9: {"T1": 97.5, "T2": 61.0, "T4": 34.5},   # N=2^18, T=125
}

# Distributed MPI (Table 2)
mpi_data = {
    6: {"T1": 37.9, "T2": 24.7, "T4": 18.8},   # N=2^12, S=8000
    7: {"T1": 44.5, "T2": 25.4, "T4": 16.3},   # N=2^14, S=2000
    8: {"T1": 59.4, "T2": 34.9, "T4": 20.1},   # N=2^16, S=500
    9: {"T1": 92.4, "T2": 65.9, "T4": 26.8},   # N=2^18, S=125
}

# Compute empirical speedups
sm_rows = []
for K in K_values:
    v = sm_data[K]
    N = 2**(2*K)
    su2 = v["T1"] / v["T2"]
    su4 = v["T1"] / v["T4"]
    sm_rows.append([K, N, v["T1"], v["T2"], v["T4"], su2, su4])

mpi_rows = []
for K in K_values:
    v = mpi_data[K]
    N = 2**(2*K)
    su2 = v["T1"] / v["T2"]
    su4 = v["T1"] / v["T4"]
    mpi_rows.append([K, N, v["T1"], v["T2"], v["T4"], su2, su4])

print("Empirical Speedup Data:")
print("\nShared Memory:")
print(f"  {'K':>3} {'N':>8} {'SU(P=2)':>8} {'SU(P=4)':>8}")
for row in sm_rows:
    print(f"  {int(row[0]):3d} {int(row[1]):8d} {row[5]:8.2f} {row[6]:8.2f}")

print("\nDistributed Memory (MPI):")
print(f"  {'K':>3} {'N':>8} {'SU(P=2)':>8} {'SU(P=4)':>8}")
for row in mpi_rows:
    print(f"  {int(row[0]):3d} {int(row[1]):8d} {row[5]:8.2f} {row[6]:8.2f}")

# Save empirical data
np.savetxt("../data/fig3_empirical_shared.csv",
           np.array(sm_rows), delimiter=",",
           header="K,N,T1_sec,T2_sec,T4_sec,SU_P2,SU_P4",
           fmt="%.8e")

np.savetxt("../data/fig3_empirical_mpi.csv",
           np.array(mpi_rows), delimiter=",",
           header="K,N,T1_sec,T2_sec,T4_sec,SU_P2,SU_P4",
           fmt="%.8e")

# Theoretical curves using Eq. (6) with fitted parameters
# Shared memory: non-monotonic speedup, peak at K=8
# Fit: xi_sm=4.5, f_sm=0.035 gives reasonable qualitative match
xi_sm, f_sm = 4.5, 0.035
# MPI: monotonically increasing speedup for P=4
xi_mpi, f_mpi = 5.0, 0.005

K_fine = np.linspace(5.5, 9.5, 200)
theory_rows = []
for K in K_fine:
    su_sm_2 = theoretical_speedup(K, 2, xi_sm, f_sm)
    su_sm_4 = theoretical_speedup(K, 4, xi_sm, f_sm)
    su_mpi_2 = theoretical_speedup(K, 2, xi_mpi, f_mpi)
    su_mpi_4 = theoretical_speedup(K, 4, xi_mpi, f_mpi)
    theory_rows.append([K, 2**(2*K), su_sm_2, su_sm_4, su_mpi_2, su_mpi_4])

np.savetxt("../data/fig3_theoretical_speedup.csv",
           np.array(theory_rows), delimiter=",",
           header="K,N,SU_shared_P2,SU_shared_P4,SU_MPI_P2,SU_MPI_P4",
           fmt="%.8e")

# Plot
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Left: Shared Memory
ax = axes[0]
K_emp = [row[0] for row in sm_rows]
su2_emp = [row[5] for row in sm_rows]
su4_emp = [row[6] for row in sm_rows]

su_sm_2_fine = [theoretical_speedup(K, 2, xi_sm, f_sm) for K in K_fine]
su_sm_4_fine = [theoretical_speedup(K, 4, xi_sm, f_sm) for K in K_fine]

ax.plot(K_fine, su_sm_2_fine, "b-", alpha=0.5, label=f"Eq.(6) P=2")
ax.plot(K_fine, su_sm_4_fine, "r-", alpha=0.5, label=f"Eq.(6) P=4")
ax.plot(K_emp, su2_emp, "bs", markersize=10, label="Paper P=2")
ax.plot(K_emp, su4_emp, "r^", markersize=10, label="Paper P=4")
ax.axhline(y=2, color="b", linestyle=":", alpha=0.2)
ax.axhline(y=4, color="r", linestyle=":", alpha=0.2)
ax.set_xlabel("K  (array size N = $2^{2K}$)")
ax.set_ylabel("Speedup (T$_1$ / T$_P$)")
ax.set_title("Shared Memory Speedup\n(Origin 200, $doacross directives)")
ax.set_xticks(K_values)
ax.set_xticklabels([f"{K}\n(N=$2^{{{2*K}}}$)" for K in K_values])
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)
ax.set_ylim(0, 5)

# Right: MPI
ax = axes[1]
K_emp = [row[0] for row in mpi_rows]
su2_mpi_emp = [row[5] for row in mpi_rows]
su4_mpi_emp = [row[6] for row in mpi_rows]

su_mpi_2_fine = [theoretical_speedup(K, 2, xi_mpi, f_mpi) for K in K_fine]
su_mpi_4_fine = [theoretical_speedup(K, 4, xi_mpi, f_mpi) for K in K_fine]

ax.plot(K_fine, su_mpi_2_fine, "b-", alpha=0.5, label=f"Eq.(6) P=2")
ax.plot(K_fine, su_mpi_4_fine, "r-", alpha=0.5, label=f"Eq.(6) P=4")
ax.plot(K_emp, su2_mpi_emp, "bs", markersize=10, label="Paper P=2")
ax.plot(K_emp, su4_mpi_emp, "r^", markersize=10, label="Paper P=4")
ax.axhline(y=2, color="b", linestyle=":", alpha=0.2)
ax.axhline(y=4, color="r", linestyle=":", alpha=0.2)
ax.set_xlabel("K  (array size N = $2^{2K}$)")
ax.set_ylabel("Speedup (T$_1$ / T$_P$)")
ax.set_title("Distributed Memory (MPI) Speedup\n(Origin 200, MPI 1.1)")
ax.set_xticks(K_values)
ax.set_xticklabels([f"{K}\n(N=$2^{{{2*K}}}$)" for K in K_values])
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)
ax.set_ylim(0, 5)

plt.tight_layout()
plt.savefig("../plots/fig3_speedup_analysis.png", dpi=150)
print("Saved fig3_speedup_analysis.png")
