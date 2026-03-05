"""
Fig 3: 4-pole and 5-pole Pade approximants of R(zeta)
Hunana et al. (2019), Figure 3

Top: 4-pole (R_{4,0} to R_{4,4})
Bottom: 5-pole (R_{5,0} to R_{5,6})
Left: Im R(zeta), Right: Re R(zeta) for real zeta
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, '.')
from plasma_dispersion import (R_exact, R_4_0, R_4_1, R_4_2, R_4_3, R_4_4,
                                R_5_0, R_5_1, R_5_2, R_5_3, R_5_4, R_5_5, R_5_6)

zeta = np.linspace(0, 4, 500)
R_ex = R_exact(zeta + 0j)

approx_4 = [
    (R_4_0, r'$R_{4,0}$', 'r--'),
    (R_4_1, r'$R_{4,1}$', 'g:'),
    (R_4_2, r'$R_{4,2}$', 'b-.'),
    (R_4_3, r'$R_{4,3}$', 'm--'),
    (R_4_4, r'$R_{4,4}$', 'c:'),
]
approx_5 = [
    (R_5_0, r'$R_{5,0}$', 'r--'),
    (R_5_2, r'$R_{5,2}$', 'b-.'),
    (R_5_3, r'$R_{5,3}$', 'g:'),
    (R_5_5, r'$R_{5,5}$', 'm--'),
    (R_5_6, r'$R_{5,6}$', 'c:'),
]

fig, axes = plt.subplots(2, 2, figsize=(12, 9))

for row, approxes, title_prefix in [(0, approx_4, '4-pole'), (1, approx_5, '5-pole')]:
    for col, part, label in [(0, 'imag', r'Im $R(\zeta)$'), (1, 'real', r'Re $R(\zeta)$')]:
        ax = axes[row, col]
        vals = np.imag(R_ex) if part == 'imag' else np.real(R_ex)
        ax.plot(zeta, vals, 'k-', linewidth=2, label=r'Exact $R(\zeta)$')
        for func, name, style in approxes:
            v = func(zeta + 0j)
            vv = np.imag(v) if part == 'imag' else np.real(v)
            ax.plot(zeta, vv, style, linewidth=1.5, label=name)
        ax.set_xlabel(r'$\zeta$', fontsize=12)
        ax.set_ylabel(label, fontsize=12)
        ax.legend(fontsize=9, ncol=2)
        ax.set_title(f'{title_prefix}: {label}', fontsize=11)
        ax.grid(True, alpha=0.3)

fig.suptitle('Fig. 3: Pad\u00e9 approximants of R(\u03b6) - 4 and 5 poles', fontsize=13)
plt.tight_layout()
plt.savefig('../plots/fig3_pade_high.png', dpi=200, bbox_inches='tight')
plt.close()
print("Saved ../plots/fig3_pade_high.png")

# Export CSV
with open('../data/fig3_pade_high.csv', 'w') as f:
    f.write("# Figure 3: Pade approximants of R(zeta), 4-5 poles\n")
    f.write("zeta,ImR_exact,ReR_exact")
    for n, np_ in [(4,0),(4,1),(4,2),(4,3),(4,4),(5,0),(5,2),(5,3),(5,5),(5,6)]:
        f.write(f",ImR_{n}_{np_},ReR_{n}_{np_}")
    f.write("\n")

    funcs_list = [R_4_0, R_4_1, R_4_2, R_4_3, R_4_4, R_5_0, R_5_2, R_5_3, R_5_5, R_5_6]
    for i in range(len(zeta)):
        vals = [zeta[i], np.imag(R_ex[i]), np.real(R_ex[i])]
        for func in funcs_list:
            v = func(zeta[i] + 0j)
            vals.extend([np.imag(v), np.real(v)])
        f.write(",".join(f"{v:.8f}" for v in vals) + "\n")
print("Exported ../data/fig3_pade_high.csv")
