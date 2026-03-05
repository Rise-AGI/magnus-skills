"""
Fig 4: Percentage error of Im R_{3,0}(zeta) and Im R_{3,1}(zeta)
Hunana et al. (2019), Figure 4
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, '.')
from plasma_dispersion import R_exact, R_3_0, R_3_1

zeta = np.linspace(0.01, 4, 500)
R_ex = R_exact(zeta + 0j)

errors = {}
for name, func, color in [('R_{3,0}', R_3_0, 'r'), ('R_{3,1}', R_3_1, 'g')]:
    v = func(zeta + 0j)
    err_im = np.abs(np.imag(v) - np.imag(R_ex)) / np.abs(R_ex) * 100
    errors[name] = err_im

fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(zeta, errors['R_{3,0}'], 'r-', linewidth=2, label=r'$R_{3,0}(\zeta)$')
ax.plot(zeta, errors['R_{3,1}'], 'g-', linewidth=2, label=r'$R_{3,1}(\zeta)$')
ax.set_xlabel(r'$\zeta$', fontsize=13)
ax.set_ylabel('% error of Im part', fontsize=13)
ax.legend(fontsize=12)
ax.set_title(r'Fig. 4: % error of Im $R_{n,n\prime}(\zeta)$', fontsize=13)
ax.grid(True, alpha=0.3)
ax.set_xlim(0, 4)
ax.set_ylim(0, max(np.max(errors['R_{3,0}']), np.max(errors['R_{3,1}'])) * 1.1)

plt.tight_layout()
plt.savefig('../plots/fig4_pade_error.png', dpi=200, bbox_inches='tight')
plt.close()
print("Saved ../plots/fig4_pade_error.png")

# Export CSV
with open('../data/fig4_pade_error.csv', 'w') as f:
    f.write("# Figure 4: Percentage error of Im R_{3,0} and R_{3,1}\n")
    f.write("zeta,error_Im_R_3_0,error_Im_R_3_1\n")
    for i in range(len(zeta)):
        f.write(f"{zeta[i]:.8f},{errors['R_{3,0}'][i]:.8f},{errors['R_{3,1}'][i]:.8f}\n")
print("Exported ../data/fig4_pade_error.csv")
