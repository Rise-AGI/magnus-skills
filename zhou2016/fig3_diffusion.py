"""
Figure 3: Quasi-linear diffusion coefficient continuity at gamma=0.

Demonstrates the key result of the paper: the general solution
(Eq. 9, 10, 17) yields a quasi-linear diffusion coefficient (Eq. 19)
that is continuous across gamma=0, unlike the standard Vlasov result
which has a sign discontinuity (Eq. 4).

We show D(v) for several values of gamma approaching 0 from both
positive and negative sides.

Output: data/fig3_diffusion.csv, plots/fig3_diffusion.png
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(".")))
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from vlasov_landau import PlasmaParams

params = PlasmaParams(n0=1e18, T_eV=10.0)

# Parameters for a Langmuir wave with k*lambda_D = 0.3
kld = 0.3
k = kld / params.lambda_D
omega_pe = params.omega_pe
Omega = omega_pe * np.sqrt(1 + 3 * kld**2)  # real frequency
v_ph = Omega / k  # phase velocity

# Velocity grid centered on phase velocity
v_th = params.v_th
v_array = np.linspace(v_ph - 5 * v_th, v_ph + 5 * v_th, 500)

# Gamma values: positive (growing), zero, negative (damped)
gamma_values = [0.05 * omega_pe, 0.01 * omega_pe, 0.001 * omega_pe,
                -0.001 * omega_pe, -0.01 * omega_pe, -0.05 * omega_pe]
gamma_labels = [r"$\gamma = 0.05\omega_{pe}$", r"$\gamma = 0.01\omega_{pe}$",
                r"$\gamma = 0.001\omega_{pe}$",
                r"$\gamma = -0.001\omega_{pe}$", r"$\gamma = -0.01\omega_{pe}$",
                r"$\gamma = -0.05\omega_{pe}$"]

# --- Standard Vlasov diffusion (Eq. 3): D_V = gamma / ((Omega - kv)^2 + gamma^2) ---
# This has sign(gamma) discontinuity

# --- General solution diffusion (Eq. 19): continuous ---
# For gamma > 0: same as Vlasov (S=0)
# For gamma < 0: Lorentzian + 2*pi/k * delta correction
# The total is always positive and -> pi*delta(kv-Omega) as gamma->0

D_vlasov_all = []
D_general_all = []

for gamma in gamma_values:
    D_vlasov = gamma / ((Omega - k * v_array)**2 + gamma**2)

    if gamma > 0:
        D_general = D_vlasov.copy()
    else:
        # General solution correction term
        # The correction makes D positive and continuous
        # D_general = gamma/((kv-Omega)^2+gamma^2) + (2*pi/k) * delta_reg(v - Omega/k)
        # Use Gaussian regularization with width = |gamma/k|
        sigma_v = abs(gamma / k)
        delta_reg = np.exp(-0.5 * ((v_array - v_ph) / sigma_v)**2) / (
            sigma_v * np.sqrt(2 * np.pi))
        D_general = D_vlasov + (2 * np.pi / k) * delta_reg

    D_vlasov_all.append(D_vlasov)
    D_general_all.append(D_general)

# Save data: velocity + D for each gamma value
v_normalized = (v_array - v_ph) / v_th
cols = [v_normalized]
col_names = ["v_minus_vph_over_vth"]
for i, gamma in enumerate(gamma_values):
    cols.append(D_vlasov_all[i] * k)  # normalize by k for dimensionless
    cols.append(D_general_all[i] * k)
    gname = f"gamma_{gamma/omega_pe:.3f}"
    col_names.append(f"D_vlasov_{gname}")
    col_names.append(f"D_general_{gname}")

data = np.column_stack(cols)
np.savetxt("../data/fig3_diffusion.csv", data, delimiter=",",
           header=", ".join(col_names), fmt="%.8e")
print("Saved data/fig3_diffusion.csv")

# Plot
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

colors_pos = ['#1f77b4', '#2ca02c', '#ff7f0e']
colors_neg = ['#ff7f0e', '#2ca02c', '#1f77b4']

# Left: Vlasov (discontinuous)
for i in range(3):
    ax1.plot(v_normalized, D_vlasov_all[i] * k, color=colors_pos[i],
             linewidth=1.5, label=gamma_labels[i])
for i in range(3, 6):
    ax1.plot(v_normalized, D_vlasov_all[i] * k, color=colors_neg[i-3],
             linewidth=1.5, linestyle='--', label=gamma_labels[i])
ax1.axhline(y=0, color='k', linewidth=0.5)
ax1.set_xlabel(r"$(v - v_{ph}) / v_{th}$", fontsize=13)
ax1.set_ylabel(r"$k \cdot \mathcal{D}$ (normalized)", fontsize=13)
ax1.set_title("Vlasov Diffusion (Eq. 3)\nDiscontinuous at $\\gamma = 0$", fontsize=13)
ax1.legend(fontsize=9, loc='upper right')
ax1.set_xlim(-4, 4)
ax1.grid(True, alpha=0.3)

# Right: General solution (continuous)
for i in range(3):
    ax2.plot(v_normalized, D_general_all[i] * k, color=colors_pos[i],
             linewidth=1.5, label=gamma_labels[i])
for i in range(3, 6):
    ax2.plot(v_normalized, D_general_all[i] * k, color=colors_neg[i-3],
             linewidth=1.5, linestyle='--', label=gamma_labels[i])
ax2.axhline(y=0, color='k', linewidth=0.5)
ax2.set_xlabel(r"$(v - v_{ph}) / v_{th}$", fontsize=13)
ax2.set_ylabel(r"$k \cdot \mathcal{D}$ (normalized)", fontsize=13)
ax2.set_title("General Solution Diffusion (Eq. 19)\nContinuous at $\\gamma = 0$", fontsize=13)
ax2.legend(fontsize=9, loc='upper right')
ax2.set_xlim(-4, 4)
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("../plots/fig3_diffusion.png", dpi=150, bbox_inches="tight")
print("Saved plots/fig3_diffusion.png")
