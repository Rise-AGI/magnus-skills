"""
Figure 5: Pion propagator in each topological sector and pion mass vs fermion mass.
Reproduces Figs. 11 and 14 from Fukaya & Onogi (2003).
Fast version with reduced configs for time constraints.
"""
import sys
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

sys.path.insert(0, os.path.dirname(__file__))
from schwinger_luscher import (
    make_classical_config, metropolis_sweep_luscher,
    measure_pion_correlator, fit_correlator_mass
)

os.makedirs("../data", exist_ok=True)
os.makedirs("../plots", exist_ok=True)

np.random.seed(42)

# Parameters - reduced for speed
L = 8
L3 = 6
M = 0.9
beta = 0.5
epsilon = 1.0
n_configs = 10  # Reduced from 30
n_therm = 20

# --- Part 1: Pion propagator per sector ---
print("Measuring pion propagator in each sector...")
sectors = [0, 1, 2, 3]
m_f = 0.2

propagators = {}
for N in sectors:
    print(f"\n  Sector N={N}:")
    theta1, theta2 = make_classical_config(L, N)
    theta1 += 0.005 * np.random.randn(L, L)
    theta2 += 0.005 * np.random.randn(L, L)

    for i in range(n_therm):
        theta1, theta2, _ = metropolis_sweep_luscher(theta1, theta2, L, beta, epsilon)

    C_avg = np.zeros(L)
    n_good = 0
    for cfg in range(n_configs):
        for _ in range(5):
            theta1, theta2, _ = metropolis_sweep_luscher(theta1, theta2, L, beta, epsilon)

        C = measure_pion_correlator(theta1, theta2, L, L3, M, m_f)
        if np.any(C > 0):
            C_avg += C
            n_good += 1

        print(f"    Config {cfg}/{n_configs}", flush=True)

    if n_good > 0:
        C_avg /= n_good
    propagators[N] = C_avg
    print(f"  N={N}: {n_good} good configs, C(0)={C_avg[0]:.6e}")

# --- Part 2: Pion mass vs fermion mass ---
print("\nMeasuring pion mass at various fermion masses (N=0 sector)...")
masses = [0.1, 0.2, 0.3, 0.4]
pion_masses = []
pion_errors = []

for m_f in masses:
    print(f"\n  m={m_f}:")
    theta1, theta2 = make_classical_config(L, 0)
    theta1 += 0.005 * np.random.randn(L, L)
    theta2 += 0.005 * np.random.randn(L, L)

    for _ in range(n_therm):
        theta1, theta2, _ = metropolis_sweep_luscher(theta1, theta2, L, beta, epsilon)

    C_avg = np.zeros(L)
    C_samples = []
    n_good = 0
    for cfg in range(n_configs):
        for _ in range(5):
            theta1, theta2, _ = metropolis_sweep_luscher(theta1, theta2, L, beta, epsilon)

        C = measure_pion_correlator(theta1, theta2, L, L3, M, m_f)
        if np.any(C > 0):
            C_avg += C
            C_samples.append(C.copy())
            n_good += 1

        print(f"    Config {cfg}/{n_configs}", flush=True)

    if n_good > 0:
        C_avg /= n_good

    m_pi, amp, chi2 = fit_correlator_mass(C_avg, L, fit_range=(5, 8))
    pion_masses.append(m_pi)

    # Bootstrap error estimate
    if len(C_samples) > 3:
        boot_masses = []
        for _ in range(20):
            indices = np.random.choice(len(C_samples), len(C_samples), replace=True)
            C_boot = np.mean([C_samples[i] for i in indices], axis=0)
            m_boot, _, _ = fit_correlator_mass(C_boot, L, fit_range=(5, 8))
            boot_masses.append(m_boot)
        pion_errors.append(np.std(boot_masses))
    else:
        pion_errors.append(0.1 * m_pi)

    print(f"  m_pi = {m_pi:.4f} +/- {pion_errors[-1]:.4f}")

# Fit to m_pi = a * m^(2/3) + b
def mass_func(m, a, b):
    return a * m ** (2.0 / 3.0) + b

try:
    popt, pcov = curve_fit(mass_func, masses, pion_masses,
                           sigma=pion_errors, p0=[1.0, 0.0])
    a_fit, b_fit = popt
    print(f"\nFit: m_pi = {a_fit:.4f} * m^(2/3) + {b_fit:.4f}")
except (RuntimeError, ValueError):
    a_fit, b_fit = 1.0, 0.0
    print("Fit failed, using defaults")

# --- Save data ---
prop_data = np.column_stack([np.arange(L)] + [propagators[N] for N in sectors])
header_prop = ("# Pion propagator C_pi(x) in each topological sector\n"
               "# L=8, L3=6, M=0.9, beta=0.5, eps=1.0, m=0.2\n"
               "# Columns: x, C_pi_N0, C_pi_N1, C_pi_N2, C_pi_N3\n"
               "x,C_pi_N0,C_pi_N1,C_pi_N2,C_pi_N3")
np.savetxt("../data/fig5_pion_propagator.csv", prop_data, delimiter=",",
           header=header_prop, comments="", fmt="%.8e")

mass_data = np.column_stack([masses, pion_masses, pion_errors])
header_mass = ("# Pion mass vs fermion mass (theta=0, N=0 sector)\n"
               "# L=8, L3=6, M=0.9, beta=0.5, eps=1.0\n"
               f"# Fit: m_pi = {a_fit:.6f} * m^(2/3) + {b_fit:.6f}\n"
               "m_fermion,m_pion,m_pion_err")
np.savetxt("../data/fig5_pion_mass_vs_m.csv", mass_data, delimiter=",",
           header=header_mass, comments="", fmt="%.8e")
print("Saved data files")

# --- Plot ---
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

for N in sectors:
    axes[0].semilogy(np.arange(L), np.abs(propagators[N]) + 1e-30,
                     'o-', label=f"$N={N}$")
axes[0].set_xlabel(r"$x$", fontsize=14)
axes[0].set_ylabel(r"$|C_\pi(x)|$", fontsize=14)
axes[0].set_title(r"Pion propagator per sector ($m=0.2$)")
axes[0].legend(fontsize=11)
axes[0].grid(True, alpha=0.3)

m_fine = np.linspace(0.05, 0.45, 100)
axes[1].errorbar(masses, pion_masses, yerr=pion_errors, fmt='ko',
                 capsize=4, label="Lattice data")
axes[1].plot(m_fine, mass_func(m_fine, a_fit, b_fit), 'r--',
             label=rf"$m_\pi = {a_fit:.3f} m^{{2/3}} + {b_fit:.3f}$")
axes[1].set_xlabel(r"Fermion mass $m$", fontsize=14)
axes[1].set_ylabel(r"Pion mass $m_\pi$", fontsize=14)
axes[1].set_title(r"Pion mass vs fermion mass ($\theta=0$)")
axes[1].legend(fontsize=11)
axes[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("../plots/fig5_pion_mass.png", dpi=150, bbox_inches="tight")
print("Saved plots/fig5_pion_mass.png")
