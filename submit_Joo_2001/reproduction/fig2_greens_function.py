"""
Fig 2: Green's function G(tau) for coexisting metallic and insulating solutions.
Analogous to Fig. 3 of Joo & Oudovenko (2001).
"""
import sys, os
sys.path.insert(0, ".")
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from dmft_ipt import dmft_loop

beta = 50.0
U = 2.8
n_freq = 256
n_tau = 1024

print(f"Computing coexistence at beta={beta}, U={U}...")
r_met = dmft_loop(U, beta, n_freq=n_freq, n_tau=n_tau, max_iter=80, tol=1e-7, seed='metallic')
r_ins = dmft_loop(U, beta, n_freq=n_freq, n_tau=n_tau, max_iter=80, tol=1e-7, seed='insulating')

tau = r_met['tau']

fig, ax = plt.subplots(figsize=(8, 6))
ax.plot(tau, -r_met['G_tau'], 'b-', linewidth=2, label='Metallic')
ax.plot(tau, -r_ins['G_tau'], 'r--', linewidth=2, label='Insulating')
ax.set_xlabel(r'$\tau$', fontsize=14)
ax.set_ylabel(r'$-G(\tau)$', fontsize=14)
ax.set_title(f'Fig. 2: Coexisting $G(\\tau)$ ($\\beta={int(beta)}$, $U={U}$)', fontsize=13)
ax.legend(fontsize=12)
ax.grid(True, alpha=0.3)
ax.set_xlim(0, beta)

plt.tight_layout()
plt.savefig('../plots/fig2_greens_function.png', dpi=150, bbox_inches='tight')
print("Saved ../plots/fig2_greens_function.png")
plt.close()

n_out = 200
indices = np.linspace(0, len(tau) - 1, n_out, dtype=int)
with open('../data/fig2_greens_function.csv', 'w') as f:
    f.write(f"# Figure 2: G(tau) coexistence, beta={beta}, U={U}\n")
    f.write("tau,neg_G_tau_metallic,neg_G_tau_insulating\n")
    for i in indices:
        f.write(f"{tau[i]:.8f},{-r_met['G_tau'][i]:.8f},{-r_ins['G_tau'][i]:.8f}\n")
print("Exported ../data/fig2_greens_function.csv")
