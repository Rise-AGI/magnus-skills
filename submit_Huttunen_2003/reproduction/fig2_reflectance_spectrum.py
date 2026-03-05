"""
Figure 2: Reflectance spectrum of a finite 1D photonic crystal stack.

Shows the reflectance vs normalized frequency for finite PhC stacks with
different numbers of periods. This demonstrates how the photonic band gap
manifests as high reflectance in a realistic finite structure.

Parameters: eps_rod=13 (GaAs), eps_air=1, l/P=0.2
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(sys.argv[0])) if sys.argv[0] else '.')
from photonic_crystal import finite_stack_reflectance, find_bandgaps

script_dir = os.path.dirname(os.path.abspath(sys.argv[0])) if sys.argv[0] else '.'
data_dir = os.path.join(script_dir, '..', 'data')
plots_dir = os.path.join(script_dir, '..', 'plots')
os.makedirs(data_dir, exist_ok=True)
os.makedirs(plots_dir, exist_ok=True)

# Parameters
eps_rod = 13.0
eps_air = 1.0
l_frac = 0.2
d_air = 1.0 - l_frac

# Frequency range
n_omega = 500
omega_max = 0.6
omegas = np.linspace(0.01, omega_max, n_omega)

# Number of periods to compare
n_periods_list = [3, 5, 10, 20]
colors = ['lightblue', 'blue', 'darkblue', 'black']

# Find band gaps
gaps = find_bandgaps(eps_rod, eps_air, l_frac, d_air, n_omega=3000, omega_max=omega_max)

# Compute reflectance for each n_periods
all_reflectances = {}
for n_per in n_periods_list:
    R_values = []
    for omega in omegas:
        R = finite_stack_reflectance(omega, eps_rod, eps_air, l_frac, d_air,
                                      n_per, eps_in=1.0, eps_out=1.0)
        R_values.append(R)
    all_reflectances[n_per] = np.array(R_values)
    print(f"N={n_per}: R(0.3)={all_reflectances[n_per][np.argmin(np.abs(omegas-0.3))]:.4f}")

# Save data
save_arr = [omegas]
header_parts = ['omega_normalized']
for n_per in n_periods_list:
    save_arr.append(all_reflectances[n_per])
    header_parts.append(f'R_N{n_per}')

save_data = np.column_stack(save_arr)
np.savetxt(os.path.join(data_dir, 'fig2_reflectance_spectrum.csv'),
           save_data, delimiter=',',
           header=','.join(header_parts),
           comments='# Reflectance spectrum of finite 1D PhC stack\n'
                    '# eps_rod=13, eps_air=1, l/P=0.2\n# ')

# Plot
fig, ax = plt.subplots(1, 1, figsize=(10, 6))

for i, n_per in enumerate(n_periods_list):
    ax.plot(omegas, all_reflectances[n_per], color=colors[i],
            linewidth=1.0 + 0.3 * i, label=f'N = {n_per} periods')

# Shade band gaps
for gl, gh in gaps:
    ax.axvspan(gl, gh, alpha=0.1, color='red', label='Band gap (infinite)' if gl == gaps[0][0] else '')

ax.axvline(x=0.3, color='green', linestyle='--', alpha=0.5, linewidth=1.5,
           label=r'$\omega P/(2\pi c) = 0.3$')

ax.set_xlabel(r'Normalized frequency $\omega P / (2\pi c)$')
ax.set_ylabel('Reflectance')
ax.set_title('Reflectance spectrum of finite 1D photonic crystal\n'
             r'$\varepsilon_{rod} = 13$ (GaAs), $\varepsilon_{air} = 1$, $l/P = 0.2$')
ax.legend(loc='upper right')
ax.grid(True, alpha=0.3)
ax.set_xlim(0, omega_max)
ax.set_ylim(-0.02, 1.05)

fig.savefig(os.path.join(plots_dir, 'fig2_reflectance_spectrum.png'), dpi=150, bbox_inches='tight')
plt.close(fig)

print("Done: fig2_reflectance_spectrum.png")
