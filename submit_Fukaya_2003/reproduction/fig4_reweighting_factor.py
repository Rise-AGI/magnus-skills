"""
Figure 4: Total reweighting factor R^N(beta, m) vs topological charge N.
Reproduces Fig. 6 from Fukaya & Onogi (2003).

R^N(beta,m) = exp(-beta*S_Gmin_pure^N) * Det^N * (integral correction)

The dominant terms are the exponential of the classical action and the
fermion determinant ratio. The integral correction from gauge fluctuations
(S_subtr) is a higher-order effect computed via MC.
"""
import sys
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(__file__))
from schwinger_luscher import (
    luscher_action_min, luscher_action, make_classical_config,
    metropolis_sweep_luscher, compute_DetN
)

os.makedirs("../data", exist_ok=True)
os.makedirs("../plots", exist_ok=True)

np.random.seed(42)

# Parameters
L = 16
L3 = 6
M = 0.9
m_f = 0.2
epsilon = 1.0
beta_target = 0.5

N_values = [0, 1, 2, 3, 4]

# --- Step 1: S_Gmin_pure^N (without beta factor) ---
print("Computing S_Gmin_pure^N...")
S_Gmin_pure = {}
for N in N_values:
    S_Gmin_pure[N] = luscher_action_min(L, N, 1.0, epsilon)
    print(f"  S_Gmin_pure^{N} = {S_Gmin_pure[N]:.8f}")

# --- Step 2: Det^N ---
print("\nComputing Det^N (m=0.2)...")
DetN = {}
for N in N_values:
    DetN[N] = compute_DetN(L, L3, M, m_f, N, n_nu=3)
    print(f"  Det^{N} = {DetN[N]:.8e}")

# --- Step 3: S_subtr from MC at a few beta values ---
print("\nComputing S_subtr^N via MC (quenched approximation)...")
beta_primes = [0.5, 1.0, 1.5, 2.0]
n_configs = 300
n_therm = 50

mean_SG_pure = {N: {} for N in N_values}
for bp in beta_primes:
    print(f"\n  beta' = {bp}:")
    for N in N_values:
        theta1, theta2 = make_classical_config(L, N)
        theta1 += 0.005 * np.random.randn(L, L)
        theta2 += 0.005 * np.random.randn(L, L)

        for _ in range(n_therm):
            theta1, theta2, _ = metropolis_sweep_luscher(theta1, theta2, L, bp, epsilon)

        S_vals = []
        for _ in range(n_configs):
            theta1, theta2, _ = metropolis_sweep_luscher(theta1, theta2, L, bp, epsilon)
            S_full = luscher_action(theta1, theta2, L, bp, epsilon)
            if S_full < np.inf:
                S_vals.append(S_full / bp)

        mean_SG_pure[N][bp] = np.mean(S_vals) if S_vals else S_Gmin_pure[N]
        print(f"    N={N}: <S_G_pure> = {mean_SG_pure[N][bp]:.4f}")

# Compute integral of S_subtr
integral_S = {}
for N in N_values:
    vals = []
    for bp in beta_primes:
        s_subtr = (mean_SG_pure[N][bp] - S_Gmin_pure[N]) - mean_SG_pure[0][bp]
        vals.append(s_subtr)
    bps_ext = beta_primes + [3.0]
    vals_ext = vals + [0.0]
    integral_S[N] = np.trapezoid(vals_ext, bps_ext)

# --- Step 4: Total R^N ---
print("\nComponents of R^N:")
R_N = {}
for N in N_values:
    exp_action = np.exp(-beta_target * S_Gmin_pure[N])
    exp_integral = np.exp(integral_S[N])
    R_N[N] = exp_action * DetN[N] * exp_integral
    print(f"  N={N}: exp(-bS_min)={exp_action:.4f}, Det={DetN[N]:.4f}, "
          f"exp(int)={exp_integral:.4f}")

# Normalize
R0 = R_N[0] if R_N[0] > 0 else 1.0
for N in N_values:
    R_N[N] /= R0

print("\nR^N(0.5, 0.2) normalized:")
for N in N_values:
    print(f"  R^{N} = {R_N[N]:.6f}")

# Paper values for comparison
paper_R = {0: 1.0, 1: 0.637, 2: 0.442, 3: 0.201, 4: 0.064}
print("\nPaper values:")
for N in N_values:
    print(f"  R^{N} = {paper_R[N]:.3f}")

# --- Save data ---
data = np.column_stack([
    N_values,
    [S_Gmin_pure[N] for N in N_values],
    [DetN[N] for N in N_values],
    [integral_S[N] for N in N_values],
    [R_N[N] for N in N_values]
])
header = ("# Total reweighting factor R^N(beta=0.5, m=0.2)\n"
          "# L=16, L3=6, M=0.9, epsilon=1.0\n"
          "# R^N = exp(-beta*S_Gmin_pure) * Det^N * exp(S_subtr_integral)\n"
          "# Columns: N, S_Gmin_pure_N, Det_N, integral_S_subtr, R_N\n"
          "N,S_Gmin_pure_N,Det_N,integral_S_subtr,R_N")
np.savetxt("../data/fig4_reweighting_factor.csv", data, delimiter=",",
           header=header, comments="", fmt="%.8e")
print("Saved data/fig4_reweighting_factor.csv")

# --- Plot ---
fig, ax = plt.subplots(figsize=(8, 6))
R_vals = [R_N[N] for N in N_values]
paper_vals = [paper_R[N] for N in N_values]
x = np.arange(len(N_values))
width = 0.35
ax.bar(x - width/2, R_vals, width, color='steelblue', alpha=0.8,
       edgecolor='black', label='This work')
ax.bar(x + width/2, paper_vals, width, color='coral', alpha=0.8,
       edgecolor='black', label='Paper (Table I)')
ax.set_xlabel(r"Topological charge $N$", fontsize=14)
ax.set_ylabel(r"$R^N(0.5, 0.2)$", fontsize=14)
ax.set_title(r"Reweighting factor $R^N(\beta=0.5, m=0.2)$")
ax.set_xticks(x)
ax.set_xticklabels(N_values)
ax.legend(fontsize=12)
ax.grid(True, alpha=0.3, axis='y')
plt.tight_layout()
plt.savefig("../plots/fig4_reweighting_factor.png", dpi=150, bbox_inches="tight")
print("Saved plots/fig4_reweighting_factor.png")
