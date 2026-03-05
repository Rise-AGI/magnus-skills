"""
Figure 8: Avalanche multiplication factor vs inverse aspect ratio
With analytic estimate from Eq. A.4
"""
import sys
sys.path.insert(0, '.')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from physics import avalanche_toroidicity_factor

E_over_Ec = 5.0
epsilon = np.linspace(0, 0.5, 200)
factor = avalanche_toroidicity_factor(epsilon, E_over_Ec)

# Save data
header = "epsilon, gamma_A_over_gamma_A_cyl"
np.savetxt('../data/fig8_avalanche_toroidicity.csv',
           np.column_stack([epsilon, factor]),
           delimiter=',', header=header, comments='', fmt='%.8e')

fig, ax = plt.subplots(figsize=(8, 6))
ax.plot(epsilon, factor, 'b-', linewidth=2,
        label=f'Analytic (Eq. A.4), $E/E_c = {E_over_Ec}$')

# Also plot for other E/Ec values
for Eec, col in [(10, 'green'), (20, 'orange'), (40, 'red')]:
    f = avalanche_toroidicity_factor(epsilon, Eec)
    ax.plot(epsilon, f, color=col, linewidth=1.5, linestyle='--',
            label=f'$E/E_c = {Eec}$')

ax.set_xlabel(r'Inverse aspect ratio $\epsilon = r/R$', fontsize=14)
ax.set_ylabel(r'$\bar{\gamma}_A / \bar{\gamma}_{A,cyl}$', fontsize=14)
ax.set_title('Figure 8: Avalanche multiplication factor vs inverse aspect ratio', fontsize=14)
ax.legend(fontsize=11)
ax.grid(True, alpha=0.3)
ax.set_xlim([0, 0.5])
ax.set_ylim([0, 1.1])

plt.tight_layout()
plt.savefig('../plots/fig8_avalanche_toroidicity.png', dpi=150)
print("Figure 8 saved.")
