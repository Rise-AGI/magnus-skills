"""
Figure 10: Time for 1% of electrons to run away vs E/Ec
For different temperatures.
"""
import sys
sys.path.insert(0, '.')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from physics import time_to_fraction

n_e = 1e19
ln_Lambda = 15.0

E_over_Ec_vals = np.logspace(np.log10(2), np.log10(200), 60)

temperatures = [
    (500, 'b--', r'$T_e = 0.5$ keV'),
    (2000, 'g-', r'$T_e = 2$ keV'),
    (5000, 'r-', r'$T_e = 5$ keV'),
]

fig, ax = plt.subplots(figsize=(8, 6))

all_data = [E_over_Ec_vals]
header_parts = ["E_over_Ec"]

for T_eV, style, label in temperatures:
    times = []
    for Eec in E_over_Ec_vals:
        t = time_to_fraction(Eec, T_eV, n_e, ln_Lambda,
                            target_frac=0.01, include_avalanche=True,
                            max_tau_th=1e15)
        times.append(t)
    times = np.array(times)
    all_data.append(times)
    header_parts.append(f"time_tau_th_Te_{T_eV}eV")

    # Filter finite values
    mask = np.isfinite(times) & (times > 0)
    ax.semilogy(E_over_Ec_vals[mask], times[mask], style, linewidth=2, label=label)

# Save data
np.savetxt('../data/fig10_timescale.csv',
           np.column_stack(all_data),
           delimiter=',', header=','.join(header_parts), comments='', fmt='%.8e')

ax.set_xlabel(r'$E / E_c$', fontsize=14)
ax.set_ylabel(r'Time to 1% runaway $[t / \tau_{th}]$', fontsize=14)
ax.set_title('Figure 10: Time for 1% runaway generation', fontsize=14)
ax.legend(fontsize=12)
ax.grid(True, alpha=0.3)
ax.set_xlim([2, 200])
ax.set_xscale('log')

plt.tight_layout()
plt.savefig('../plots/fig10_timescale.png', dpi=150)
print("Figure 10 saved.")
