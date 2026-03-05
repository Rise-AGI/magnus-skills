"""
Figure 9: Fraction of runaways from knock-on collisions (nA/nr)
as 2D parametric scan in (E/Ec, Te) space.
Plus analytic contour lines from Eq. 21.
"""
import sys
sys.path.insert(0, '.')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from physics import (growth_rate_ratio_analytic, compute_avalanche_fraction,
                     dreicer_rate, avalanche_growth_rate, thermal_collision_freq)

ln_Lambda = 15.0
n_e = 1e19
nr_over_ne = 0.01

# Parameter grid
E_over_Ec_vals = np.logspace(np.log10(2), np.log10(200), 40)
T_eV_vals = np.logspace(np.log10(50), np.log10(10000), 35)  # 0.05 to 10 keV

# Compute avalanche fraction on grid
fraction_grid = np.zeros((len(T_eV_vals), len(E_over_Ec_vals)))

print("Computing parametric scan (this may take a while)...")
for i, T_eV in enumerate(T_eV_vals):
    for j, E_over_Ec in enumerate(E_over_Ec_vals):
        try:
            frac = compute_avalanche_fraction(E_over_Ec, T_eV, n_e, ln_Lambda,
                                              target_frac=0.01)
            fraction_grid[i, j] = frac
        except Exception:
            fraction_grid[i, j] = np.nan
    if (i + 1) % 5 == 0:
        print(f"  Row {i+1}/{len(T_eV_vals)} done")

# Save data
header_line = "T_eV_rows_E_over_Ec_cols"
np.savetxt('../data/fig9_parametric_scan.csv', fraction_grid,
           delimiter=',', header=header_line, comments='', fmt='%.6e')
np.savetxt('../data/fig9_E_over_Ec.csv', E_over_Ec_vals, fmt='%.6e')
np.savetxt('../data/fig9_T_eV.csv', T_eV_vals, fmt='%.6e')

# Analytic contour lines for gamma_A/gamma_D from Eq. 21
E_cont = np.logspace(np.log10(2), np.log10(200), 100)

def find_T_for_ratio(target_ratio, E_over_Ec, ln_Lambda, nr_ne=0.01):
    """Find Te where gamma_A/gamma_D = target_ratio via bisection."""
    T_lo, T_hi = 10, 20000
    for _ in range(100):
        T_mid = (T_lo + T_hi) / 2
        ratio = growth_rate_ratio_analytic(E_over_Ec, T_mid, ln_Lambda, nr_ne)
        if ratio > target_ratio:
            T_lo = T_mid
        else:
            T_hi = T_mid
    return (T_lo + T_hi) / 2

# Compute contour lines for 5%, 50%, 90% avalanche dominance
# These correspond approximately to gamma_A/gamma_D ratios
contour_ratios = [(0.05, 'cyan', '5%'), (1.0, 'yellow', '50%'), (9.0, 'red', '90%')]
contour_data = {}
for ratio_val, color, label in contour_ratios:
    T_line = []
    E_line = []
    for Eec in E_cont:
        try:
            T = find_T_for_ratio(ratio_val, Eec, ln_Lambda, nr_over_ne)
            if 30 < T < 15000:
                T_line.append(T)
                E_line.append(Eec)
        except Exception:
            pass
    contour_data[label] = (np.array(E_line), np.array(T_line))

# Experimental data points from the paper
# (E/Ec, Te_eV, tokamak_label)
exp_data = [
    (11, 3800, 'TS #40719'),
    (85, 400, 'COMPASS #8555'),
    (94, 350, 'COMPASS #8630'),
    (12, 1000, 'Alcator C-Mod'),
    (8, 2500, 'DIII-D'),
    (6, 800, 'FTU'),
    (10, 1200, 'TEXTOR'),
    (5, 500, 'KSTAR'),
]

# Plot
fig, ax = plt.subplots(figsize=(10, 8))

# Colormesh
E_mesh, T_mesh = np.meshgrid(E_over_Ec_vals, T_eV_vals)
pcm = ax.pcolormesh(E_mesh, T_mesh / 1000, fraction_grid * 100,
                     cmap='RdYlBu_r', shading='auto', vmin=0, vmax=100)
cbar = plt.colorbar(pcm, ax=ax, label=r'$n_A / n_r$ [%]')

# Contour lines
for ratio_val, color, label in contour_ratios:
    E_line, T_line = contour_data[label]
    if len(E_line) > 0:
        ax.plot(E_line, T_line / 1000, color=color, linewidth=2.5,
                label=f'Analytic {label} avalanche')

# Experimental points
for Eec, T, tok in exp_data:
    ax.plot(Eec, T / 1000, 'ko', markersize=8)
    ax.annotate(tok, (Eec, T / 1000), textcoords="offset points",
                xytext=(5, 5), fontsize=7)

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$E / E_c$', fontsize=14)
ax.set_ylabel(r'$T_e$ [keV]', fontsize=14)
ax.set_title('Figure 9: Avalanche fraction in (E/Ec, Te) parameter space', fontsize=14)
ax.legend(fontsize=10, loc='upper right')
ax.set_xlim([2, 200])
ax.set_ylim([0.05, 10])

plt.tight_layout()
plt.savefig('../plots/fig9_parametric_scan.png', dpi=150)
print("Figure 9 saved.")
