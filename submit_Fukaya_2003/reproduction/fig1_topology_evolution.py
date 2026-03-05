"""
Figure 1: Monte Carlo evolution of topological charge.
Comparison of Wilson's gauge action vs Luscher's gauge action.
Reproduces Figs. 1 and 2 from Fukaya & Onogi (2003).

Left panel: Wilson action (beta=2.0) shows topology changes.
Right panel: Luscher action (beta=0.5, epsilon=1.0) preserves topology.
Also: evolution starting from N=2 sector with Luscher action.
"""
import sys
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(os.getcwd())))
sys.path.insert(0, os.path.dirname(__file__))
from schwinger_luscher import (
    make_classical_config, topological_charge,
    metropolis_sweep_wilson, metropolis_sweep_luscher
)

os.makedirs("../data", exist_ok=True)
os.makedirs("../plots", exist_ok=True)

np.random.seed(42)

L = 16
n_traj = 500

# --- Panel 1: Wilson action, beta=2.0, start from N=0 ---
print("Running Wilson action MC (beta=2.0)...")
beta_W = 2.0
theta1, theta2 = make_classical_config(L, 0)
# Add small random perturbation
theta1 += 0.01 * np.random.randn(L, L)
theta2 += 0.01 * np.random.randn(L, L)

Q_wilson = []
for i in range(n_traj):
    theta1, theta2, acc = metropolis_sweep_wilson(theta1, theta2, L, beta_W)
    Q = topological_charge(theta1, theta2, L)
    Q_wilson.append(Q)
    if i % 50 == 0:
        print(f"  Wilson sweep {i}: Q={Q}, acc={acc:.3f}")

# --- Panel 2: Luscher action, beta=0.5, epsilon=1.0, start from N=0 ---
print("Running Luscher action MC (beta=0.5, eps=1.0), N=0...")
beta_L = 0.5
epsilon = 1.0
theta1, theta2 = make_classical_config(L, 0)
theta1 += 0.01 * np.random.randn(L, L)
theta2 += 0.01 * np.random.randn(L, L)

Q_luscher_0 = []
for i in range(n_traj):
    theta1, theta2, acc = metropolis_sweep_luscher(theta1, theta2, L, beta_L, epsilon)
    Q = topological_charge(theta1, theta2, L)
    Q_luscher_0.append(Q)
    if i % 50 == 0:
        print(f"  Luscher(N=0) sweep {i}: Q={Q}, acc={acc:.3f}")

# --- Panel 3: Luscher action, start from N=2 ---
print("Running Luscher action MC (beta=0.5, eps=1.0), N=2...")
theta1, theta2 = make_classical_config(L, 2)
theta1 += 0.005 * np.random.randn(L, L)
theta2 += 0.005 * np.random.randn(L, L)

Q_luscher_2 = []
for i in range(n_traj):
    theta1, theta2, acc = metropolis_sweep_luscher(theta1, theta2, L, beta_L, epsilon)
    Q = topological_charge(theta1, theta2, L)
    Q_luscher_2.append(Q)
    if i % 50 == 0:
        print(f"  Luscher(N=2) sweep {i}: Q={Q}, acc={acc:.3f}")

# --- Save data ---
data = np.column_stack([
    np.arange(n_traj),
    Q_wilson,
    Q_luscher_0,
    Q_luscher_2
])
header = ("# Topological charge evolution: Wilson vs Luscher action\n"
          "# Columns: sweep, Q_wilson(beta=2.0), Q_luscher_N0(beta=0.5,eps=1.0), Q_luscher_N2(beta=0.5,eps=1.0)\n"
          "# L=16, 500 sweeps\n"
          "sweep,Q_wilson,Q_luscher_N0,Q_luscher_N2")

np.savetxt("../data/fig1_topology_evolution.csv", data, delimiter=",",
           header=header, comments="", fmt=["%d", "%d", "%d", "%d"])
print("Saved data/fig1_topology_evolution.csv")

# --- Plot ---
fig, axes = plt.subplots(1, 3, figsize=(15, 4))

axes[0].plot(Q_wilson, 'b-', linewidth=0.5)
axes[0].set_xlabel("MC sweep")
axes[0].set_ylabel("Topological charge Q")
axes[0].set_title(r"Wilson action, $\beta=2.0$")
axes[0].set_ylim(-5, 5)

axes[1].plot(Q_luscher_0, 'r-', linewidth=0.5)
axes[1].set_xlabel("MC sweep")
axes[1].set_ylabel("Topological charge Q")
axes[1].set_title(r"L\"uscher action, $\beta=0.5$, $\epsilon=1.0$, $N=0$")
axes[1].set_ylim(-5, 5)

axes[2].plot(Q_luscher_2, 'g-', linewidth=0.5)
axes[2].set_xlabel("MC sweep")
axes[2].set_ylabel("Topological charge Q")
axes[2].set_title(r"L\"uscher action, $\beta=0.5$, $\epsilon=1.0$, $N=2$")
axes[2].set_ylim(-1, 5)

plt.tight_layout()
plt.savefig("../plots/fig1_topology_evolution.png", dpi=150, bbox_inches="tight")
print("Saved plots/fig1_topology_evolution.png")
