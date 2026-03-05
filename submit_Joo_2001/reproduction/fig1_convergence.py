"""
Fig 1: DMFT iteration convergence at a point in the coexistence region.
Shows how metallic and insulating seeds converge to distinct solutions.
Analogous to Fig. 1 of Joo & Oudovenko (2001).
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

print(f"Running DMFT-IPT: beta={beta}, U={U}")
r_met = dmft_loop(U, beta, n_freq=n_freq, n_tau=n_tau,
                   max_iter=80, tol=1e-8, mixing=0.5,
                   seed='metallic', verbose=True)
r_ins = dmft_loop(U, beta, n_freq=n_freq, n_tau=n_tau,
                   max_iter=80, tol=1e-8, mixing=0.5,
                   seed='insulating', verbose=True)

print(f"Metallic: n_iter={r_met['n_iter']}, ImSig={r_met['im_sigma_history'][-1]:.6f}")
print(f"Insulating: n_iter={r_ins['n_iter']}, ImSig={r_ins['im_sigma_history'][-1]:.6f}")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

iters_met = np.arange(1, len(r_met['im_sigma_history']) + 1)
iters_ins = np.arange(1, len(r_ins['im_sigma_history']) + 1)

ax1.plot(iters_met, r_met['im_sigma_history'], 'b-o', markersize=3, label='Metallic seed')
ax1.set_xlabel('Iteration', fontsize=12)
ax1.set_ylabel(r'Im$\Sigma(i\pi T)$', fontsize=12)
ax1.set_title(f'(a) Metallic solution ($\\beta={int(beta)}$, $U={U}$)', fontsize=12)
ax1.legend(fontsize=10)
ax1.grid(True, alpha=0.3)

ax2.plot(iters_ins, r_ins['im_sigma_history'], 'r-s', markersize=3, label='Insulating seed')
ax2.set_xlabel('Iteration', fontsize=12)
ax2.set_ylabel(r'Im$\Sigma(i\pi T)$', fontsize=12)
ax2.set_title(f'(b) Insulating solution ($\\beta={int(beta)}$, $U={U}$)', fontsize=12)
ax2.legend(fontsize=10)
ax2.grid(True, alpha=0.3)

plt.suptitle('Fig. 1: Convergence of DMFT-IPT iterations', fontsize=14, y=1.02)
plt.tight_layout()
plt.savefig('../plots/fig1_convergence.png', dpi=150, bbox_inches='tight')
print("Saved ../plots/fig1_convergence.png")
plt.close()

with open('../data/fig1_convergence.csv', 'w') as f:
    f.write(f"# Figure 1: DMFT-IPT convergence, beta={beta}, U={U}\n")
    f.write("iteration,ImSigma_metallic,ImSigma_insulating\n")
    n_max = max(len(r_met['im_sigma_history']), len(r_ins['im_sigma_history']))
    for i in range(n_max):
        met_val = f"{r_met['im_sigma_history'][i]:.8f}" if i < len(r_met['im_sigma_history']) else ''
        ins_val = f"{r_ins['im_sigma_history'][i]:.8f}" if i < len(r_ins['im_sigma_history']) else ''
        f.write(f"{i+1},{met_val},{ins_val}\n")
print("Exported ../data/fig1_convergence.csv")
