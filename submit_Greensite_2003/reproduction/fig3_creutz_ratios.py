"""
Figure 3: Creutz ratios for string tension extraction.

The Creutz ratio chi(R,R) = -ln[W(R,R)*W(R-1,R-1) / W(R,R-1)^2]
provides an estimate of the string tension that eliminates perimeter
contributions.

Related to Greensite (2003) Section 4.1 (linearity of the confining force).
"""
import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) if os.path.dirname(os.path.abspath(__file__)) else '.')
from su2_lattice import SU2Lattice, strong_coupling_string_tension

# Parameters
L = 8
n_therm = 20
n_meas = 8
seed = 42
beta_values = [1.5, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0]

print("Figure 3: Creutz ratios")
print(f"Lattice: {L}^4, thermalization: {n_therm}, measurements: {n_meas}")

results = []

for idx, beta in enumerate(beta_values):
    print(f"  beta = {beta:.1f} ... ", flush=True)
    lattice = SU2Lattice(L, beta, seed=seed + idx * 100)
    lattice.initialize_cold()

    for _ in range(n_therm):
        lattice.sweep()

    # Measure Wilson loops for sizes 1x1, 2x2, 1x2, 2x1, 2x3, 3x3
    loop_measurements = {}
    for R in range(1, 4):
        for T in range(R, 4):
            loop_measurements[(R, T)] = []

    for m in range(n_meas):
        lattice.sweep()
        for R in range(1, 4):
            for T in range(R, 4):
                w = lattice.wilson_loop(R, T)
                loop_measurements[(R, T)].append(w)
        print(f"    measurement {m+1}/{n_meas}")

    # Compute Creutz ratios
    loop_means = {}
    for key, vals in loop_measurements.items():
        loop_means[key] = np.mean(vals)
        # Also store symmetric version
        if key[0] != key[1]:
            loop_means[(key[1], key[0])] = loop_means[key]

    creutz_ratios = {}
    for R in range(2, 4):
        try:
            w_rr = loop_means[(R, R)]
            w_r1r1 = loop_means[(R-1, R-1)]
            w_rr1 = loop_means[(R, R-1)]
            if w_rr > 0 and w_r1r1 > 0 and w_rr1 > 0:
                chi = -np.log(w_rr * w_r1r1 / w_rr1**2)
                creutz_ratios[R] = chi
                print(f"    chi({R},{R}) = {chi:.6f}")
        except (KeyError, ValueError):
            pass

    results.append({
        'beta': beta,
        'loops': loop_means,
        'creutz': creutz_ratios
    })

# Save data
data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data')
os.makedirs(data_dir, exist_ok=True)

with open(os.path.join(data_dir, 'fig3_creutz_ratios.csv'), 'w') as f:
    f.write("# Figure 3: Creutz ratios chi(R,R)\n")
    f.write("# L=8, n_therm=20, n_meas=8\n")
    f.write("beta,chi_2x2,chi_3x3\n")
    for r in results:
        chi2 = r['creutz'].get(2, float('nan'))
        chi3 = r['creutz'].get(3, float('nan'))
        f.write(f"{r['beta']:.1f},{chi2:.8f},{chi3:.8f}\n")

# Plot
plot_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'plots')
os.makedirs(plot_dir, exist_ok=True)

fig, ax = plt.subplots(figsize=(8, 6))

betas = [r['beta'] for r in results]
chi2 = [r['creutz'].get(2, float('nan')) for r in results]
chi3 = [r['creutz'].get(3, float('nan')) for r in results]

ax.plot(betas, chi2, 'bo-', label=r'$\chi(2,2)$', markersize=6)
ax.plot(betas, chi3, 'rs-', label=r'$\chi(3,3)$', markersize=6)

# Strong coupling prediction
beta_sc = np.linspace(0.5, 3.5, 100)
sc_tension = np.array([strong_coupling_string_tension(b) for b in beta_sc])
sc_tension = np.clip(sc_tension, 0, 5)
ax.plot(beta_sc[sc_tension > 0], sc_tension[sc_tension > 0], 'k--',
        label='Strong coupling', linewidth=1.5)

ax.set_xlabel(r'$\beta = 4/g^2$', fontsize=14)
ax.set_ylabel(r'$\chi(R,R) \approx a^2 K$', fontsize=14)
ax.set_title('Creutz Ratios (String Tension)', fontsize=14)
ax.set_yscale('log')
ax.set_ylim(0.01, 5)
ax.legend(fontsize=12)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(plot_dir, 'fig3_creutz_ratios.png'), dpi=150)
print("Saved fig3_creutz_ratios.png")
