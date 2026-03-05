"""
Figure S_luscher: Luscher phase shift function and resonance energy levels.

Computes the kinematic function phi(q) from Luscher (Ref. S8) used in the
quantization condition for two-particle states in a finite box, and
demonstrates the pi-pi scattering phase shift in the rho-meson channel.

The basic result of Luscher (S7) is that the finite-volume energy spectrum
is given by W = 2*sqrt(Mpi^2 + k^2) with k solving:
    n*pi - delta_11(k) = phi(q),  where q = kL/(2*pi)

Output:
    ../plots/figS_luscher.png
    ../data/figS_luscher.csv
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, ".")
from lattice_qcd import (
    luscher_phi, pipi_phase_shift_rho, EXPERIMENTAL_MASSES, DECAY_WIDTHS,
)

print("=" * 60)
print("Luscher Phase Shift and Resonance Analysis")
print("=" * 60)

M_pi = EXPERIMENTAL_MASSES["pi"]   # 0.135 GeV
M_rho = EXPERIMENTAL_MASSES["rho"]  # 0.775 GeV
Gamma_rho = DECAY_WIDTHS["rho"]     # 0.149 GeV

k_rho = np.sqrt(M_rho**2 / 4.0 - M_pi**2)
print(f"M_pi = {M_pi:.3f} GeV")
print(f"M_rho = {M_rho:.3f} GeV")
print(f"Gamma_rho = {Gamma_rho:.3f} GeV")
print(f"k_rho = {k_rho:.4f} GeV")

# ---- Panel 1: phi(q) ----
q_values = np.linspace(0.01, 2.0, 500)
phi_approx = luscher_phi(q_values, use_approx=True)

# ---- Panel 2: pi-pi phase shift ----
k_values = np.linspace(0.01, 0.5, 500)
delta_11 = pipi_phase_shift_rho(k_values, M_rho, Gamma_rho, M_pi)

# ---- Panel 3: Energy levels vs L ----
L_values = np.linspace(1.0, 8.0, 200)  # fm
L_inv_GeV = L_values / 0.19733  # convert fm to 1/GeV

# Free two-pion energy levels (no interaction)
n_max = 3
free_energies = {}
for n_vec_sq in [0, 1, 2, 3]:  # |n|^2 = 0,1,2,3
    k_n = np.sqrt(n_vec_sq) * 2 * np.pi / L_inv_GeV
    W_n = 2.0 * np.sqrt(M_pi**2 + k_n**2)
    free_energies[n_vec_sq] = W_n

# ---- Plot ----
fig, axes = plt.subplots(1, 3, figsize=(16, 5))

# Panel 1: phi(q)
ax = axes[0]
ax.plot(q_values, phi_approx, "b-", linewidth=2, label=r"$\phi(q)$ (approximation)")
ax.plot(q_values, q_values**3, "r--", alpha=0.5, label=r"$q^3$ (small $q$)")
ax.plot(q_values, np.pi * q_values**2, "g--", alpha=0.5, label=r"$\pi q^2$ (large $q$)")
ax.set_xlabel(r"$q = kL/(2\pi)$", fontsize=12)
ax.set_ylabel(r"$\phi(q)$", fontsize=12)
ax.set_title(r"L\"uscher kinematic function $\phi(q)$", fontsize=13)
ax.legend(fontsize=9)
ax.set_ylim(0, 12)
ax.grid(alpha=0.3)

# Panel 2: Phase shift
ax = axes[1]
ax.plot(k_values * 1000, np.degrees(delta_11), "b-", linewidth=2)
ax.axhline(90, color="gray", linestyle="--", alpha=0.5, label=r"$\delta = 90°$")
ax.axvline(k_rho * 1000, color="red", linestyle=":", alpha=0.7,
           label=f"$k_\\rho$ = {k_rho*1000:.0f} MeV")
ax.set_xlabel(r"$k$ [MeV]", fontsize=12)
ax.set_ylabel(r"$\delta_{11}$ [degrees]", fontsize=12)
ax.set_title(r"$\pi\pi$ scattering phase shift ($I=1, J=1$)", fontsize=13)
ax.legend(fontsize=9)
ax.grid(alpha=0.3)

# Panel 3: Energy levels
ax = axes[2]
ax.axhline(M_rho, color="red", linewidth=2, linestyle="-", label=r"$M_\rho$ (resonance)")
for n_sq, W_n in free_energies.items():
    mult = {0: r"$|\mathbf{n}|^2=0$", 1: r"$|\mathbf{n}|^2=1$",
            2: r"$|\mathbf{n}|^2=2$", 3: r"$|\mathbf{n}|^2=3$"}
    ax.plot(L_values, W_n, "--", alpha=0.6, label=f"Free: {mult[n_sq]}")

ax.axhline(2 * M_pi, color="black", linestyle=":", alpha=0.3, label=r"$2M_\pi$ threshold")
ax.set_xlabel(r"$L$ [fm]", fontsize=12)
ax.set_ylabel(r"$W$ [GeV]", fontsize=12)
ax.set_title("Two-pion energy levels in finite volume", fontsize=13)
ax.set_ylim(0.2, 1.5)
ax.legend(fontsize=8, loc="upper right")
ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig("../plots/figS_luscher.png", dpi=300)
print("\nSaved: ../plots/figS_luscher.png")

# ---- Export CSV ----
with open("../data/figS_luscher.csv", "w") as f:
    f.write("# Luscher phase shift function and pi-pi scattering\n")
    f.write("# Part 1: phi(q) function\n")
    f.write("q,phi_approx\n")
    for q, p in zip(q_values, phi_approx):
        f.write(f"{q:.8f},{p:.8f}\n")

with open("../data/figS_luscher_phase.csv", "w") as f:
    f.write("# pi-pi scattering phase shift in rho channel (I=1, J=1)\n")
    f.write("# Columns: k_MeV, delta_11_deg\n")
    f.write("k_MeV,delta_11_deg\n")
    for k, d in zip(k_values, delta_11):
        f.write(f"{k*1000:.8f},{np.degrees(d):.8f}\n")

print("Saved: ../data/figS_luscher.csv")
print("Saved: ../data/figS_luscher_phase.csv")
print("Done.")
