"""
Fig 3: Double occupancy scan across coexistence region.

Scans U through the coexistence region at fixed temperature, showing
distinct metallic and insulating branches of the double occupancy
<n_up n_down>. This is analogous to the double occupancy data
discussed in Fig. 4 of Joo & Oudovenko (2001).

Parameters: beta=64.
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(".")))
sys.path.insert(0, ".")

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from dmft_ipt import dmft_loop

beta = 64.0
n_freq = 512
n_tau = 2048

# Scan U values across the coexistence region
U_values = np.arange(2.0, 3.01, 0.05)

print(f"Scanning U from {U_values[0]:.2f} to {U_values[-1]:.2f} at beta={beta}")

d_met = []
d_ins = []
im_sig_met = []
im_sig_ins = []

for U in U_values:
    print(f"  U = {U:.2f} ...", end='', flush=True)
    r_m = dmft_loop(U, beta, n_freq=n_freq, n_tau=n_tau,
                     max_iter=100, tol=1e-7, seed='metallic')
    r_i = dmft_loop(U, beta, n_freq=n_freq, n_tau=n_tau,
                     max_iter=100, tol=1e-7, seed='insulating')
    d_met.append(r_m['double_occ'])
    d_ins.append(r_i['double_occ'])
    im_sig_met.append(r_m['im_sigma_history'][-1])
    im_sig_ins.append(r_i['im_sigma_history'][-1])
    print(f" d_met={r_m['double_occ']:.5f}, d_ins={r_i['double_occ']:.5f}")

d_met = np.array(d_met)
d_ins = np.array(d_ins)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

# Double occupancy
ax1.plot(U_values, d_met, 'bo-', markersize=5, label='Metallic')
ax1.plot(U_values, d_ins, 'rs-', markersize=5, label='Insulating')
ax1.set_xlabel('$U$', fontsize=13)
ax1.set_ylabel(r'$\langle n_\uparrow n_\downarrow \rangle$', fontsize=13)
ax1.set_title(f'(a) Double occupancy ($\\beta={int(beta)}$)', fontsize=12)
ax1.legend(fontsize=11)
ax1.grid(True, alpha=0.3)

# ImSigma(i*pi*T) — order parameter
ax2.plot(U_values, np.abs(im_sig_met), 'bo-', markersize=5, label='Metallic')
ax2.plot(U_values, np.abs(im_sig_ins), 'rs-', markersize=5, label='Insulating')
ax2.set_xlabel('$U$', fontsize=13)
ax2.set_ylabel(r'$|$Im$\Sigma(i\pi T)|$', fontsize=13)
ax2.set_title(f'(b) Self-energy at first Matsubara freq.', fontsize=12)
ax2.legend(fontsize=11)
ax2.grid(True, alpha=0.3)

plt.suptitle('Fig. 3: Coexistence region scan', fontsize=14, y=1.02)
plt.tight_layout()
plt.savefig('../plots/fig3_double_occupancy.png', dpi=300, bbox_inches='tight')
print("Saved ../plots/fig3_double_occupancy.png")
plt.close()

# Export CSV
with open('../data/fig3_double_occupancy.csv', 'w') as f:
    f.write(f"# Figure 3: Double occupancy scan, beta={beta}\n")
    f.write("# Columns: U, d_metallic, d_insulating, ImSig_metallic, ImSig_insulating\n")
    f.write("U,d_metallic,d_insulating,ImSig_metallic,ImSig_insulating\n")
    for i in range(len(U_values)):
        f.write(f"{U_values[i]:.8f},{d_met[i]:.8f},{d_ins[i]:.8f},"
                f"{im_sig_met[i]:.8f},{im_sig_ins[i]:.8f}\n")
print("Exported ../data/fig3_double_occupancy.csv")
