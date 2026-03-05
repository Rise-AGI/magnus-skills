"""
Fig 4: Phase diagram of the paramagnetic metal-insulator transition.
Computes IPT coexistence boundaries at multiple temperatures and overlays
QMC data from Joo & Oudovenko (2001).
"""
import sys, os
sys.path.insert(0, ".")
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from dmft_ipt import dmft_loop


def find_boundaries(beta, U_range, n_freq=256, n_tau=1024):
    """Find Uc1 and Uc2 by detecting where metallic and insulating solutions diverge."""
    im_g_met = []
    im_g_ins = []
    for U in U_range:
        r_m = dmft_loop(U, beta, n_freq=n_freq, n_tau=n_tau, max_iter=60, tol=1e-6, seed='metallic')
        r_i = dmft_loop(U, beta, n_freq=n_freq, n_tau=n_tau, max_iter=60, tol=1e-6, seed='insulating')
        im_g_met.append(abs(r_m['im_g_lowest']))
        im_g_ins.append(abs(r_i['im_g_lowest']))

    im_g_met = np.array(im_g_met)
    im_g_ins = np.array(im_g_ins)
    diff = np.abs(im_g_met - im_g_ins)
    threshold = 0.05 * np.max(im_g_met)
    coex = diff > threshold
    if np.any(coex):
        idx = np.where(coex)[0]
        return U_range[idx[0]], U_range[idx[-1]]
    return None, None


beta_values = [30, 35, 40, 45, 50, 60, 80, 100, 150]
U_scan = np.arange(2.0, 3.6, 0.05)

ipt_Uc1, ipt_Uc2, ipt_T = [], [], []

print("Computing IPT coexistence boundaries...")
for beta in beta_values:
    T = 1.0 / beta
    print(f"  beta={beta} (T={T:.4f})...", end='', flush=True)
    uc1, uc2 = find_boundaries(beta, U_scan)
    if uc1 is not None:
        ipt_Uc1.append(uc1)
        ipt_Uc2.append(uc2)
        ipt_T.append(T)
        print(f" Uc1={uc1:.2f}, Uc2={uc2:.2f}")
    else:
        print(" no coexistence")

ipt_Uc1 = np.array(ipt_Uc1)
ipt_Uc2 = np.array(ipt_Uc2)
ipt_T = np.array(ipt_T)

# QMC data from the paper
qmc_data = {
    57:  (2.36, 2.40, 2.42, 2.53),
    64:  (2.36, 2.41, 2.48, 2.52),
    100: (2.35, 2.45, 2.58, 2.62),
}

fig, ax = plt.subplots(figsize=(9, 7))

# IPT coexistence region
if len(ipt_T) > 0:
    sort_idx = np.argsort(ipt_T)
    ax.fill_betweenx(ipt_T[sort_idx], ipt_Uc1[sort_idx], ipt_Uc2[sort_idx],
                      alpha=0.15, color='blue', label='IPT coexistence')
    ax.plot(ipt_Uc1[sort_idx], ipt_T[sort_idx], 'b-', linewidth=2, label=r'IPT $U_{c1}$')
    ax.plot(ipt_Uc2[sort_idx], ipt_T[sort_idx], 'b--', linewidth=2, label=r'IPT $U_{c2}$')

# QMC data
for beta_val, (uc1_lo, uc1_hi, uc2_lo, uc2_hi) in qmc_data.items():
    T_val = 1.0 / beta_val
    uc1_mid = (uc1_lo + uc1_hi) / 2
    uc2_mid = (uc2_lo + uc2_hi) / 2
    ax.errorbar(uc1_mid, T_val, xerr=(uc1_hi-uc1_lo)/2, fmt='ro', markersize=8,
                capsize=4, elinewidth=2,
                label=r'QMC $U_{c1}$' if beta_val == 57 else '')
    ax.errorbar(uc2_mid, T_val, xerr=(uc2_hi-uc2_lo)/2, fmt='r^', markersize=8,
                capsize=4, elinewidth=2,
                label=r'QMC $U_{c2}$' if beta_val == 57 else '')

# T=0 reference points
ax.plot(2.96, 0, 'g*', markersize=15, label=r'Projective $U_{c2}(T{=}0)$',
        zorder=5, markeredgecolor='black')
ax.plot(2.94, 0, 'mp', markersize=10, label=r'NRG $U_{c2}(T{=}0)$',
        zorder=5, markeredgecolor='black')
ax.plot(2.38, 1.0/47, 'ks', markersize=10, label=r'QMC Mott endpoint $T_c$', zorder=5)

ax.set_xlabel('$U / D$', fontsize=14)
ax.set_ylabel('$T / D$', fontsize=14)
ax.set_title('Fig. 4: Phase diagram — Mott-Hubbard transition', fontsize=13)
ax.set_xlim(2.0, 3.4)
ax.set_ylim(0, 0.04)
ax.legend(fontsize=9, loc='upper right')
ax.grid(True, alpha=0.3)
ax.annotate('Metal', xy=(2.15, 0.025), fontsize=14, fontweight='bold', color='blue')
ax.annotate('Insulator', xy=(3.05, 0.025), fontsize=14, fontweight='bold', color='red')

plt.tight_layout()
plt.savefig('../plots/fig4_phase_diagram.png', dpi=150, bbox_inches='tight')
print("Saved ../plots/fig4_phase_diagram.png")
plt.close()

with open('../data/fig4_phase_diagram.csv', 'w') as f:
    f.write("# Figure 4: Phase diagram\n")
    f.write("# IPT boundaries\n")
    f.write("T,Uc1_IPT,Uc2_IPT\n")
    for i in range(len(ipt_T)):
        f.write(f"{ipt_T[i]:.8f},{ipt_Uc1[i]:.8f},{ipt_Uc2[i]:.8f}\n")
    f.write("\n# QMC data from paper\n")
    f.write("beta,T,Uc1_low,Uc1_high,Uc2_low,Uc2_high\n")
    for b, (u1l, u1h, u2l, u2h) in qmc_data.items():
        f.write(f"{b},{1.0/b:.8f},{u1l},{u1h},{u2l},{u2h}\n")
print("Exported ../data/fig4_phase_diagram.csv")
