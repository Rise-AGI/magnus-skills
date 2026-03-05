"""
Fig 2: 1-pole, 2-pole, and 3-pole Pade approximants of R(zeta)
Hunana et al. (2019), Figure 2

Top: 1-pole (R_1) and 2-pole (R_{2,0})
Bottom: 3-pole (R_{3,0}, R_{3,1}, R_{3,2})
Left: Im R(zeta), Right: Re R(zeta) for real zeta
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, '.')
from plasma_dispersion import R_exact, R_1, R_2_0, R_3_0, R_3_1, R_3_2

zeta = np.linspace(0, 4, 500)

R_ex = R_exact(zeta + 0j)

approx_top = [
    (R_1, r'$R_1(\zeta)$', 'r--'),
    (R_2_0, r'$R_{2,0}(\zeta)$', 'b-.'),
]
approx_bot = [
    (R_3_0, r'$R_{3,0}(\zeta)$', 'r--'),
    (R_3_1, r'$R_{3,1}(\zeta)$', 'g:'),
    (R_3_2, r'$R_{3,2}(\zeta)$', 'b-.'),
]

fig, axes = plt.subplots(2, 2, figsize=(12, 9))

# Top row: 1-pole and 2-pole
for col, part, label in [(0, 'imag', r'Im $R(\zeta)$'), (1, 'real', r'Re $R(\zeta)$')]:
    ax = axes[0, col]
    vals = np.imag(R_ex) if part == 'imag' else np.real(R_ex)
    ax.plot(zeta, vals, 'k-', linewidth=2, label=r'Exact $R(\zeta)$')
    for func, name, style in approx_top:
        v = func(zeta + 0j)
        vv = np.imag(v) if part == 'imag' else np.real(v)
        ax.plot(zeta, vv, style, linewidth=1.5, label=name)
    ax.set_xlabel(r'$\zeta$', fontsize=12)
    ax.set_ylabel(label, fontsize=12)
    ax.legend(fontsize=10)
    ax.set_title(f'1-pole and 2-pole: {label}', fontsize=11)
    ax.grid(True, alpha=0.3)

# Bottom row: 3-pole
for col, part, label in [(0, 'imag', r'Im $R(\zeta)$'), (1, 'real', r'Re $R(\zeta)$')]:
    ax = axes[1, col]
    vals = np.imag(R_ex) if part == 'imag' else np.real(R_ex)
    ax.plot(zeta, vals, 'k-', linewidth=2, label=r'Exact $R(\zeta)$')
    for func, name, style in approx_bot:
        v = func(zeta + 0j)
        vv = np.imag(v) if part == 'imag' else np.real(v)
        ax.plot(zeta, vv, style, linewidth=1.5, label=name)
    ax.set_xlabel(r'$\zeta$', fontsize=12)
    ax.set_ylabel(label, fontsize=12)
    ax.legend(fontsize=10)
    ax.set_title(f'3-pole: {label}', fontsize=11)
    ax.grid(True, alpha=0.3)

fig.suptitle('Fig. 2: Pad\u00e9 approximants of R(\u03b6) - 1, 2, 3 poles', fontsize=13)
plt.tight_layout()
plt.savefig('../plots/fig2_pade_low.png', dpi=200, bbox_inches='tight')
plt.close()
print("Saved ../plots/fig2_pade_low.png")

# Export CSV
with open('../data/fig2_pade_low.csv', 'w') as f:
    f.write("# Figure 2: Pade approximants of R(zeta), 1-3 poles\n")
    f.write("# zeta is real\n")
    f.write("zeta,ImR_exact,ReR_exact,ImR_1,ReR_1,ImR_2_0,ReR_2_0,ImR_3_0,ReR_3_0,ImR_3_1,ReR_3_1,ImR_3_2,ReR_3_2\n")
    for i in range(len(zeta)):
        vals = [zeta[i], np.imag(R_ex[i]), np.real(R_ex[i])]
        for func in [R_1, R_2_0, R_3_0, R_3_1, R_3_2]:
            v = func(zeta[i] + 0j)
            vals.extend([np.imag(v), np.real(v)])
        f.write(",".join(f"{v:.8f}" for v in vals) + "\n")
print("Exported ../data/fig2_pade_low.csv")
