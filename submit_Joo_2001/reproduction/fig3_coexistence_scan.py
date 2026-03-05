"""
Fig 3: Spectral weight (ImG at lowest Matsubara freq) scan across coexistence.
Shows distinct metallic and insulating branches.
"""
import sys, os
sys.path.insert(0, ".")
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from dmft_ipt import dmft_loop

beta = 50.0
n_freq = 256
n_tau = 1024

U_values = np.arange(2.0, 3.41, 0.05)
print(f"Scanning U from {U_values[0]:.2f} to {U_values[-1]:.2f} at beta={beta}")

im_g_met = []
im_g_ins = []
z_met = []
z_ins = []

for U in U_values:
    r_m = dmft_loop(U, beta, n_freq=n_freq, n_tau=n_tau, max_iter=80, tol=1e-7, seed='metallic')
    r_i = dmft_loop(U, beta, n_freq=n_freq, n_tau=n_tau, max_iter=80, tol=1e-7, seed='insulating')
    im_g_met.append(r_m['im_g_lowest'])
    im_g_ins.append(r_i['im_g_lowest'])
    z_met.append(r_m['z_qp'])
    z_ins.append(r_i['z_qp'])
    print(f"  U={U:.2f}: met ImG={r_m['im_g_lowest']:.4f} Z={r_m['z_qp']:.4f} | ins ImG={r_i['im_g_lowest']:.4f} Z={r_i['z_qp']:.4f}")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

ax1.plot(U_values, np.abs(im_g_met), 'bo-', markersize=4, label='Metallic')
ax1.plot(U_values, np.abs(im_g_ins), 'rs-', markersize=4, label='Insulating')
ax1.set_xlabel('$U / D$', fontsize=13)
ax1.set_ylabel(r'$|$Im$G(i\pi T)|$', fontsize=13)
ax1.set_title(f'(a) Spectral weight ($\\beta={int(beta)}$)', fontsize=12)
ax1.legend(fontsize=11)
ax1.grid(True, alpha=0.3)
ax1.axvspan(2.55, 3.05, alpha=0.1, color='green', label='Coexistence')

ax2.plot(U_values, z_met, 'bo-', markersize=4, label='Metallic')
ax2.plot(U_values, z_ins, 'rs-', markersize=4, label='Insulating')
ax2.set_xlabel('$U / D$', fontsize=13)
ax2.set_ylabel('$Z$ (quasiparticle weight)', fontsize=13)
ax2.set_title(f'(b) Quasiparticle weight ($\\beta={int(beta)}$)', fontsize=12)
ax2.legend(fontsize=11)
ax2.grid(True, alpha=0.3)
ax2.set_ylim(-0.02, 0.5)

plt.suptitle('Fig. 3: Coexistence region scan', fontsize=14, y=1.02)
plt.tight_layout()
plt.savefig('../plots/fig3_coexistence_scan.png', dpi=150, bbox_inches='tight')
print("Saved ../plots/fig3_coexistence_scan.png")
plt.close()

with open('../data/fig3_coexistence_scan.csv', 'w') as f:
    f.write(f"# Figure 3: Coexistence scan, beta={beta}\n")
    f.write("U,ImG_metallic,ImG_insulating,Z_metallic,Z_insulating\n")
    for i in range(len(U_values)):
        f.write(f"{U_values[i]:.8f},{im_g_met[i]:.8f},{im_g_ins[i]:.8f},{z_met[i]:.8f},{z_ins[i]:.8f}\n")
print("Exported ../data/fig3_coexistence_scan.csv")
