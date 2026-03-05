"""
Figure 4: String tension vs beta - strong coupling to asymptotic freedom.

Compares the string tension extracted from Creutz ratios with:
1. Strong coupling expansion: a^2*K = -ln(beta/4)
2. Asymptotic freedom: a^2*K ~ exp(-6*pi^2/11 * (beta-2))

This is the central result that demonstrates the confining force
follows QCD asymptotic freedom predictions.

Related to Greensite (2003) Sections 4.1-4.2 and the discussion of
Creutz (1980) results.
"""
import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) if os.path.dirname(os.path.abspath(__file__)) else '.')
from su2_lattice import SU2Lattice, strong_coupling_string_tension, asymptotic_freedom_string_tension

# Parameters
L = 8
n_therm = 20
n_meas = 8
seed = 42

# Use enough beta values across both regimes
beta_values = [1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5]

print("Figure 4: String tension vs beta")
print(f"Lattice: {L}^4, thermalization: {n_therm}, measurements: {n_meas}")

string_tensions = []

for idx, beta in enumerate(beta_values):
    print(f"  beta = {beta:.1f} ... ", flush=True)
    lattice = SU2Lattice(L, beta, seed=seed + idx * 100)
    lattice.initialize_cold()

    for _ in range(n_therm):
        lattice.sweep()

    # Measure Wilson loops
    w11_list = []
    w12_list = []
    w22_list = []

    for m in range(n_meas):
        lattice.sweep()
        w11_list.append(lattice.square_wilson_loop(1))
        w12_list.append(lattice.wilson_loop(1, 2))
        w22_list.append(lattice.square_wilson_loop(2))
        print(f"    measurement {m+1}/{n_meas}")

    w11 = np.mean(w11_list)
    w12 = np.mean(w12_list)
    w22 = np.mean(w22_list)

    # Creutz ratio chi(2,2)
    if w22 > 0 and w11 > 0 and w12 > 0:
        chi = -np.log(w22 * w11 / w12**2)
    else:
        chi = float('nan')

    string_tensions.append(chi)
    print(f"    a^2*K = {chi:.6f}")

# Save data
data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data')
os.makedirs(data_dir, exist_ok=True)

with open(os.path.join(data_dir, 'fig4_string_tension.csv'), 'w') as f:
    f.write("# Figure 4: String tension (a^2*K) vs beta\n")
    f.write("# Extracted from chi(2,2) Creutz ratio\n")
    f.write("# L=8, n_therm=20, n_meas=8\n")
    f.write("beta,string_tension\n")
    for j in range(len(beta_values)):
        f.write(f"{beta_values[j]:.2f},{string_tensions[j]:.8f}\n")

# Plot
plot_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'plots')
os.makedirs(plot_dir, exist_ok=True)

fig, ax = plt.subplots(figsize=(8, 6))

# MC data
valid = [(b, s) for b, s in zip(beta_values, string_tensions)
         if not np.isnan(s) and s > 0]
if valid:
    bv, sv = zip(*valid)
    ax.plot(bv, sv, 'ko', markersize=7, label='Monte Carlo')

# Analytical curves
beta_fine = np.linspace(0.5, 3.0, 200)
sc = np.array([strong_coupling_string_tension(b) for b in beta_fine])
af = np.array([asymptotic_freedom_string_tension(b) for b in beta_fine])

mask_sc = (sc > 0.01) & (sc < 5) & (beta_fine < 2.5)
mask_af = (af > 0.01) & (af < 5) & (beta_fine > 1.8)

ax.plot(beta_fine[mask_sc], sc[mask_sc], 'b-', label='Strong coupling', linewidth=1.5)
ax.plot(beta_fine[mask_af], af[mask_af], 'r--', label='Asymptotic freedom', linewidth=1.5)

ax.set_xlabel(r'$\beta = 4/g^2$', fontsize=14)
ax.set_ylabel(r'$a^2 K$', fontsize=14)
ax.set_title('String Tension vs. Coupling', fontsize=14)
ax.set_yscale('log')
ax.set_ylim(0.01, 5)
ax.set_xlim(0.5, 3.0)
ax.legend(fontsize=12)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(plot_dir, 'fig4_string_tension.png'), dpi=150)
print("Saved fig4_string_tension.png")
