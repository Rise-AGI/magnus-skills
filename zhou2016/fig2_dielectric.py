"""
Figure 2: Plasma dielectric function comparison.

Compares Vlasov dielectric function (principal value only)
vs the Landau-corrected dielectric function at real omega.

For real omega, the Vlasov treatment gives only Re(D), while
the Landau treatment gives the full D with imaginary part,
demonstrating the role of the S*Delta correction (Eq. 17).

Output: data/fig2_dielectric.csv, plots/fig2_dielectric.png
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(".")))
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from vlasov_landau import PlasmaParams, dielectric_vlasov, dielectric_landau_real_omega

params = PlasmaParams(n0=1e18, T_eV=10.0)

# Fixed k*lambda_D = 0.3
kld = 0.3
k = kld / params.lambda_D

# Scan omega/omega_pe from 0.5 to 2.0
omega_norm = np.linspace(0.5, 2.0, 400)
omega_array = omega_norm * params.omega_pe

D_vlasov = np.zeros(len(omega_array))
D_landau_re = np.zeros(len(omega_array))
D_landau_im = np.zeros(len(omega_array))

for i, w in enumerate(omega_array):
    D_vlasov[i] = dielectric_vlasov(w, k, params)
    D_L = dielectric_landau_real_omega(w, k, params)
    D_landau_re[i] = np.real(D_L)
    D_landau_im[i] = np.imag(D_L)

# Save data
header = "omega_over_wpe, D_vlasov_real, D_landau_real, D_landau_imag"
data = np.column_stack([omega_norm, D_vlasov, D_landau_re, D_landau_im])
np.savetxt("../data/fig2_dielectric.csv", data, delimiter=",",
           header=header, fmt="%.8e")
print("Saved data/fig2_dielectric.csv")

# Plot
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8), sharex=True)

ax1.plot(omega_norm, D_vlasov, 'b-', label=r"Vlasov (P.V. only)", linewidth=2)
ax1.plot(omega_norm, D_landau_re, 'r--', label=r"Landau Re$(D)$", linewidth=2)
ax1.axhline(y=0, color='k', linewidth=0.5)
ax1.set_ylabel(r"Re$(D(\omega, k))$", fontsize=14)
ax1.legend(fontsize=11)
ax1.set_title(r"Dielectric Function at $k\lambda_D = %.1f$" % kld, fontsize=14)
ax1.grid(True, alpha=0.3)
ax1.set_ylim(-1.5, 2.5)

ax2.plot(omega_norm, D_landau_im, 'r-', label=r"Landau Im$(D)$", linewidth=2)
ax2.axhline(y=0, color='k', linewidth=0.5)
ax2.set_xlabel(r"$\omega / \omega_{pe}$", fontsize=14)
ax2.set_ylabel(r"Im$(D(\omega, k))$", fontsize=14)
ax2.legend(fontsize=11)
ax2.set_title(r"Imaginary Part (Landau Correction, $S = \pi i$ at $\gamma = 0$)", fontsize=14)
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("../plots/fig2_dielectric.png", dpi=150, bbox_inches="tight")
print("Saved plots/fig2_dielectric.png")
