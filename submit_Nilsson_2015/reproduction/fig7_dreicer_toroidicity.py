"""
Figure 7: Dreicer growth rate vs inverse aspect ratio
With fit: gamma_D/gamma_D,cyl = 1 - 1.2*sqrt(2*epsilon/(1+epsilon))
"""
import sys
sys.path.insert(0, '.')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from physics import dreicer_toroidicity_factor

epsilon = np.linspace(0, 0.6, 200)
factor = np.array([dreicer_toroidicity_factor(e) for e in epsilon])

# Save data
header = "epsilon, gamma_D_over_gamma_D_cyl"
np.savetxt('../data/fig7_dreicer_toroidicity.csv',
           np.column_stack([epsilon, factor]),
           delimiter=',', header=header, comments='', fmt='%.8e')

fig, ax = plt.subplots(figsize=(8, 6))
ax.plot(epsilon, factor, 'b-', linewidth=2,
        label=r'Fit: $1 - 1.2\sqrt{2\epsilon/(1+\epsilon)}$')

# Mark Tore Supra value
eps_ts = 0.3
ax.axvline(x=eps_ts, color='gray', linestyle=':', alpha=0.5)
ax.plot(eps_ts, dreicer_toroidicity_factor(eps_ts), 'ro', markersize=8,
        label=f'Tore Supra ($\\epsilon \\approx {eps_ts}$)')

ax.set_xlabel(r'Inverse aspect ratio $\epsilon = r/R$', fontsize=14)
ax.set_ylabel(r'$\gamma_D / \gamma_{D,cyl}$', fontsize=14)
ax.set_title('Figure 7: Dreicer growth rate vs inverse aspect ratio', fontsize=14)
ax.legend(fontsize=12)
ax.grid(True, alpha=0.3)
ax.set_xlim([0, 0.6])
ax.set_ylim([0, 1.1])

plt.tight_layout()
plt.savefig('../plots/fig7_dreicer_toroidicity.png', dpi=150)
print("Figure 7 saved.")
