"""
Figure S4: Finite-volume dependence of pion and nucleon masses.

Reproduces Figure S4 from Durr et al., Science 322, 1224 (2008)
Supplementary Online Material.

Shows the volume dependence at Mpi ~ 320 MeV, a ~ 0.125 fm for
the pion and nucleon channels, fitted to the Luscher type I formula:

    M_X(L) = M_X(inf) + c_X * exp(-Mpi*L) / (Mpi*L)^{3/2}

where c_X(Mpi) ~ Mpi^2 (predicted by Colangelo et al., Refs. S9, S10).

Output:
    ../plots/figS4_volume.png
    ../data/figS4_volume.csv
"""

import numpy as np
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, ".")
from lattice_qcd import (
    VOL_STUDY_A_FM, VOL_STUDY_MPI_MEV, VOL_STUDY_MPI_LAT,
    VOL_STUDY_L_LAT, VOL_STUDY_MPI_L,
    VOL_STUDY_PION_AM, VOL_STUDY_PION_AM_ERR,
    VOL_STUDY_NUCLEON_AM, VOL_STUDY_NUCLEON_AM_ERR,
    finite_volume_correction, HBARC,
)

print("=" * 60)
print("Figure S4: Finite-Volume Dependence")
print("=" * 60)
print(f"Lattice spacing: a = {VOL_STUDY_A_FM} fm")
print(f"Pion mass: Mpi = {VOL_STUDY_MPI_MEV} MeV")
print(f"Mpi (lattice units): {VOL_STUDY_MPI_LAT:.4f}")
print(f"Volumes (L): {VOL_STUDY_L_LAT}")
print(f"Mpi*L values: {VOL_STUDY_MPI_L}")

# ---- Fit function ----
def vol_model(L_arr, M_inf, c_X):
    """M_X(L) = M_inf + c_X * exp(-Mpi*L) / (Mpi*L)^{3/2}"""
    Mpi = VOL_STUDY_MPI_LAT
    x = Mpi * L_arr
    return M_inf + c_X * np.exp(-x) / x**1.5


# ---- Fit pion channel ----
print("\n--- Pion channel ---")
popt_pi, pcov_pi = curve_fit(
    vol_model, VOL_STUDY_L_LAT.astype(float), VOL_STUDY_PION_AM,
    sigma=VOL_STUDY_PION_AM_ERR, p0=[0.203, 0.5], absolute_sigma=True
)
M_pi_inf, c_pi = popt_pi
perr_pi = np.sqrt(np.diag(pcov_pi))
print(f"  aM_pi(inf) = {M_pi_inf:.6f} +/- {perr_pi[0]:.6f}")
print(f"  c_pi       = {c_pi:.4f} +/- {perr_pi[1]:.4f}")

# ---- Fit nucleon channel ----
print("\n--- Nucleon channel ---")
popt_N, pcov_N = curve_fit(
    vol_model, VOL_STUDY_L_LAT.astype(float), VOL_STUDY_NUCLEON_AM,
    sigma=VOL_STUDY_NUCLEON_AM_ERR, p0=[0.605, 2.0], absolute_sigma=True
)
M_N_inf, c_N = popt_N
perr_N = np.sqrt(np.diag(pcov_N))
print(f"  aM_N(inf)  = {M_N_inf:.6f} +/- {perr_N[0]:.6f}")
print(f"  c_N        = {c_N:.4f} +/- {perr_N[1]:.4f}")

# ---- Luscher prediction for c_X ~ Mpi^2 ----
print(f"\n  c_pi / Mpi^2 = {c_pi / VOL_STUDY_MPI_LAT**2:.2f}")
print(f"  c_N  / Mpi^2 = {c_N / VOL_STUDY_MPI_LAT**2:.2f}")

# ---- Plot ----
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

L_fine = np.linspace(12, 40, 200)

# Pion
ax = axes[0]
ax.errorbar(VOL_STUDY_L_LAT, VOL_STUDY_PION_AM, yerr=VOL_STUDY_PION_AM_ERR,
            fmt="o", color="C0", capsize=4, markersize=7, label="Lattice data")
ax.plot(L_fine, vol_model(L_fine, *popt_pi), "-", color="C0",
        label=f"Fit: $aM_\\pi(\\infty)$ = {M_pi_inf:.4f}")
ax.axhline(M_pi_inf, color="gray", linestyle="--", alpha=0.5, label=r"$L = \infty$")
ax.set_xlabel(r"$L/a$", fontsize=13)
ax.set_ylabel(r"$a M_\pi$", fontsize=13)
ax.set_title(r"$\pi$ volume dependence", fontsize=14)
ax.legend(fontsize=9)
ax.grid(alpha=0.3)

# Nucleon
ax = axes[1]
ax.errorbar(VOL_STUDY_L_LAT, VOL_STUDY_NUCLEON_AM, yerr=VOL_STUDY_NUCLEON_AM_ERR,
            fmt="o", color="C1", capsize=4, markersize=7, label="Lattice data")
ax.plot(L_fine, vol_model(L_fine, *popt_N), "-", color="C1",
        label=f"Fit: $aM_N(\\infty)$ = {M_N_inf:.4f}")
ax.axhline(M_N_inf, color="gray", linestyle="--", alpha=0.5, label=r"$L = \infty$")
ax.set_xlabel(r"$L/a$", fontsize=13)
ax.set_ylabel(r"$a M_N$", fontsize=13)
ax.set_title(r"$N$ volume dependence", fontsize=14)
ax.legend(fontsize=9)
ax.grid(alpha=0.3)

plt.suptitle(r"Volume dependence at $M_\pi \approx 320$ MeV, $a \approx 0.125$ fm",
             fontsize=14, y=1.02)
plt.tight_layout()
plt.savefig("../plots/figS4_volume.png", dpi=300, bbox_inches="tight")
print("\nSaved: ../plots/figS4_volume.png")

# ---- Export CSV ----
with open("../data/figS4_volume.csv", "w") as f:
    f.write("# Figure S4: Volume dependence at Mpi~320 MeV, a~0.125 fm\n")
    f.write("# Columns: L_over_a, MpiL, aM_pi, aM_pi_err, aM_N, aM_N_err\n")
    f.write("L_over_a,MpiL,aM_pi,aM_pi_err,aM_N,aM_N_err\n")
    for i in range(len(VOL_STUDY_L_LAT)):
        f.write(f"{VOL_STUDY_L_LAT[i]},{VOL_STUDY_MPI_L[i]:.8f},"
                f"{VOL_STUDY_PION_AM[i]:.8f},{VOL_STUDY_PION_AM_ERR[i]:.8f},"
                f"{VOL_STUDY_NUCLEON_AM[i]:.8f},{VOL_STUDY_NUCLEON_AM_ERR[i]:.8f}\n")

    # Append fit curves as comments
    f.write("# Fit curves (L/a, aM_pi_fit, aM_N_fit):\n")
    for L in L_fine:
        f.write(f"# {L:.8f},{vol_model(L, *popt_pi):.8f},{vol_model(L, *popt_N):.8f}\n")

print("Saved: ../data/figS4_volume.csv")

# ---- Export fit parameters ----
with open("../data/figS4_fit_params.csv", "w") as f:
    f.write("# Fit parameters for M_X(L) = M_inf + c_X * exp(-Mpi*L) / (Mpi*L)^{3/2}\n")
    f.write("# Mpi_lat = {:.6f}, a = {} fm\n".format(VOL_STUDY_MPI_LAT, VOL_STUDY_A_FM))
    f.write("channel,M_inf,M_inf_err,c_X,c_X_err\n")
    f.write(f"pion,{M_pi_inf:.8f},{perr_pi[0]:.8f},{c_pi:.8f},{perr_pi[1]:.8f}\n")
    f.write(f"nucleon,{M_N_inf:.8f},{perr_N[0]:.8f},{c_N:.8f},{perr_N[1]:.8f}\n")

print("Saved: ../data/figS4_fit_params.csv")
print("Done.")
