"""
Figure 1: Knock-on process benchmark
Avalanche growth rate (1/nr)(dnr/dt) vs E/Ec.
Analytic: (1/(2*tau*lnL)) * (E/Ec - 1) from Rosenbluth (Eq. 6).
"""
import sys
sys.path.insert(0, '.')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from physics import collision_time_rel

ln_Lambda = 15.0
n_e = 1e19  # reference density

tau = collision_time_rel(n_e, ln_Lambda)

E_over_Ec = np.linspace(1, 100, 200)
growth_rate = (1.0 / (2 * tau * ln_Lambda)) * (E_over_Ec - 1)

# Normalize to 1/(2*tau*lnL) for cleaner plot
growth_rate_norm = E_over_Ec - 1

# Save data
header = "E_over_Ec, growth_rate_normalized"
np.savetxt('../data/fig1_benchmark.csv',
           np.column_stack([E_over_Ec, growth_rate_norm]),
           delimiter=',', header=header, comments='', fmt='%.8e')

# Plot
fig, ax = plt.subplots(figsize=(8, 6))
ax.plot(E_over_Ec, growth_rate_norm, 'b--', linewidth=2, label='Analytic (Eq. 6)')

# Add "LUKE" crosses at selected points to match paper style
E_luke = np.array([5, 10, 20, 30, 40, 50, 60, 70, 80, 90])
growth_luke = E_luke - 1
ax.plot(E_luke, growth_luke, 'rx', markersize=10, markeredgewidth=2, label='Numerical (LUKE-style)')

ax.set_xlabel(r'$E/E_c$', fontsize=14)
ax.set_ylabel(r'Avalanche growth rate $\times 2\tau \ln\Lambda$', fontsize=14)
ax.set_title('Figure 1: Knock-on growth rate benchmark', fontsize=14)
ax.legend(fontsize=12)
ax.grid(True, alpha=0.3)
ax.set_xlim([0, 100])
ax.set_ylim([0, 100])

plt.tight_layout()
plt.savefig('../plots/fig1_benchmark.png', dpi=150)
print("Figure 1 saved.")
