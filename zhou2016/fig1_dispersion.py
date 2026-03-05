"""
Figure 1: Langmuir wave dispersion relation and Landau damping rate.

Computes omega(k) and gamma(k) for a Maxwellian plasma using:
  - Approximate analytical formulas (Bohm-Gross + weak damping)
  - Full numerical solution of the Landau dispersion relation

Output: data/fig1_dispersion.csv, plots/fig1_dispersion.png
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(".")))
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from vlasov_landau import PlasmaParams, solve_langmuir_dispersion, langmuir_dispersion_approx

params = PlasmaParams(n0=1e18, T_eV=10.0)

# Scan k*lambda_D from 0.05 to 0.55
kld_array = np.linspace(0.05, 0.55, 200)
k_array = kld_array / params.lambda_D

omega_r_approx = np.zeros(len(k_array))
gamma_approx = np.zeros(len(k_array))
omega_r_numerical = np.zeros(len(k_array))
gamma_numerical = np.zeros(len(k_array))

for i, k in enumerate(k_array):
    wr_a, ga = langmuir_dispersion_approx(k, params)
    omega_r_approx[i] = wr_a / params.omega_pe
    gamma_approx[i] = ga / params.omega_pe

    wr_n, gn = solve_langmuir_dispersion(k, params)
    omega_r_numerical[i] = wr_n / params.omega_pe
    gamma_numerical[i] = gn / params.omega_pe

# Save data
header = "k_lambda_D, omega_r_approx, gamma_approx, omega_r_numerical, gamma_numerical"
data = np.column_stack([kld_array, omega_r_approx, gamma_approx,
                        omega_r_numerical, gamma_numerical])
np.savetxt("../data/fig1_dispersion.csv", data, delimiter=",",
           header=header, fmt="%.8e")
print("Saved data/fig1_dispersion.csv")

# Plot
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8), sharex=True)

ax1.plot(kld_array, omega_r_approx, 'b-', label=r"Approx: $\omega_{pe}\sqrt{1+3k^2\lambda_D^2}$", linewidth=2)
ax1.plot(kld_array, omega_r_numerical, 'r--', label="Numerical (Landau)", linewidth=2)
ax1.set_ylabel(r"$\omega_r / \omega_{pe}$", fontsize=14)
ax1.legend(fontsize=11)
ax1.set_title("Langmuir Wave Dispersion Relation", fontsize=14)
ax1.grid(True, alpha=0.3)

ax2.plot(kld_array, np.abs(gamma_approx), 'b-', label="Approx (Eq. 27)", linewidth=2)
ax2.plot(kld_array, np.abs(gamma_numerical), 'r--', label="Numerical (Landau)", linewidth=2)
ax2.set_yscale("log")
ax2.set_xlabel(r"$k\lambda_D$", fontsize=14)
ax2.set_ylabel(r"$|\gamma| / \omega_{pe}$", fontsize=14)
ax2.legend(fontsize=11)
ax2.set_title("Landau Damping Rate", fontsize=14)
ax2.grid(True, alpha=0.3)
ax2.set_ylim(1e-15, 1)

plt.tight_layout()
plt.savefig("../plots/fig1_dispersion.png", dpi=150, bbox_inches="tight")
print("Saved plots/fig1_dispersion.png")
