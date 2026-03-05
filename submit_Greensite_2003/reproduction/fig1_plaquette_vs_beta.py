"""
Figure 1: Average plaquette vs beta for SU(2) lattice gauge theory.

Compares Monte Carlo results with strong coupling expansion and
weak coupling perturbative prediction.

Related to Greensite (2003) Sections 3.1, 4.1 - the average plaquette
is the most basic observable in lattice gauge theory.
"""
import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) if os.path.dirname(os.path.abspath(__file__)) else '.')
from su2_lattice import SU2Lattice, strong_coupling_plaquette, weak_coupling_plaquette

# Parameters
L = 6
n_therm = 20
n_meas = 10
seed_base = 42

beta_values = [0.5, 1.0, 1.5, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.5, 4.0]

print("Figure 1: Average plaquette vs beta")
print(f"Lattice: {L}^4, thermalization: {n_therm}, measurements: {n_meas}")

mc_plaquettes = []
mc_errors = []

for i, beta in enumerate(beta_values):
    print(f"  beta = {beta:.1f} ... ", end="", flush=True)
    lattice = SU2Lattice(L, beta, seed=seed_base + i)
    lattice.initialize_cold()

    for _ in range(n_therm):
        lattice.sweep()

    measurements = []
    for _ in range(n_meas):
        lattice.sweep()
        measurements.append(lattice.average_plaquette())

    mean_p = np.mean(measurements)
    std_p = np.std(measurements) / np.sqrt(len(measurements)) if len(measurements) > 1 else 0.0
    mc_plaquettes.append(mean_p)
    mc_errors.append(std_p)
    print(f"<P> = {mean_p:.8f} +/- {std_p:.8f}")

# Analytical curves
beta_fine = np.linspace(0.1, 5.0, 200)
sc_plaq = np.array([strong_coupling_plaquette(b) for b in beta_fine])
wc_plaq = np.array([weak_coupling_plaquette(b) for b in beta_fine])

# Clip analytical curves to physical range [0, 1]
sc_plaq = np.clip(sc_plaq, 0, 1)
wc_plaq = np.clip(wc_plaq, 0, 1)

# Save data
data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data')
os.makedirs(data_dir, exist_ok=True)

with open(os.path.join(data_dir, 'fig1_plaquette_vs_beta.csv'), 'w') as f:
    f.write("# Figure 1: Average plaquette vs beta\n")
    f.write("# MC data: L=6, n_therm=20, n_meas=10\n")
    f.write("beta,mc_plaquette,mc_error\n")
    for j in range(len(beta_values)):
        f.write(f"{beta_values[j]:.1f},{mc_plaquettes[j]:.8f},{mc_errors[j]:.8f}\n")
    f.write("\n# Analytical: strong coupling expansion\n")
    f.write("# beta,sc_plaquette\n")
    for j in range(len(beta_fine)):
        f.write(f"# {beta_fine[j]:.4f},{sc_plaq[j]:.8f}\n")

# Plot
plot_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'plots')
os.makedirs(plot_dir, exist_ok=True)

fig, ax = plt.subplots(figsize=(8, 6))
ax.errorbar(beta_values, mc_plaquettes, yerr=mc_errors, fmt='ko', markersize=6,
            label='Monte Carlo', capsize=3)
ax.plot(beta_fine, sc_plaq, 'b-', label='Strong coupling', linewidth=1.5)
ax.plot(beta_fine, wc_plaq, 'r--', label='Weak coupling', linewidth=1.5)
ax.set_xlabel(r'$\beta = 4/g^2$', fontsize=14)
ax.set_ylabel('Average Plaquette', fontsize=14)
ax.set_title('Average Plaquette vs. Coupling', fontsize=14)
ax.legend(fontsize=12)
ax.set_xlim(0, 5)
ax.set_ylim(0, 1.05)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(plot_dir, 'fig1_plaquette_vs_beta.png'), dpi=150)
print("Saved fig1_plaquette_vs_beta.png")
