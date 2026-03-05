"""
Figure 6: Conservation quantity N comparison between Split-Step and Finite Difference.

Tests conservation of N = integral(|psi|^2 dx) for a single soliton
advanced 8 time steps, for both methods and multiple S values.
Reproduces the comparison described in paper Section 2.6.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, '.')
from nls_solver import (init_one_soliton, splitstep_advance,
                         fd_forward, fd_central, compute_norm)

S_values = [-0.1, 0.0, 0.4, -10.0, 2.0]

# Split-Step parameters
ss_L, ss_N, ss_tau = 64.0, 512, 0.01
# Finite Difference parameters
fd_L, fd_N, fd_tau = 30.0, 512, 0.001

n_advance = 8
results = []

for S in S_values:
    # Split-Step
    psi_ss, x_ss = init_one_soliton(8.0, 10.0, ss_L, ss_N, S)
    N0_ss = compute_norm(psi_ss, ss_L / ss_N)
    psi_run = psi_ss.copy()
    for _ in range(n_advance):
        psi_run = splitstep_advance(psi_run, S, ss_tau, ss_L, ss_N)
    N_final_ss = compute_norm(psi_run, ss_L / ss_N)
    delta_ss = abs(N0_ss - N_final_ss)

    # Finite Difference
    psi_fd, x_fd = init_one_soliton(8.0, 10.0, fd_L, fd_N, S)
    N0_fd = compute_norm(psi_fd, fd_L / fd_N)
    h_fd = fd_L / fd_N
    psi_prev = psi_fd.copy()
    psi_curr = fd_forward(psi_prev, S, fd_tau, h_fd)
    for _ in range(n_advance - 1):
        psi_next = fd_central(psi_curr, psi_prev, S, fd_tau, h_fd)
        psi_prev = psi_curr
        psi_curr = psi_next
    N_final_fd = compute_norm(psi_curr, fd_L / fd_N)
    delta_fd = abs(N0_fd - N_final_fd)

    results.append((S, delta_ss, delta_fd))
    print(f"S={S:6.1f}: SS dN={delta_ss:.6f}, FD dN={delta_fd:.6f}")

# Plot bar chart comparison
plt.rcParams.update({'font.size': 12})
fig, ax = plt.subplots(figsize=(10, 6))
x_pos = np.arange(len(S_values))
width = 0.35
ss_vals = [r[1] for r in results]
fd_vals = [r[2] for r in results]
ax.bar(x_pos - width/2, ss_vals, width, label='Split-Step', color='steelblue')
ax.bar(x_pos + width/2, fd_vals, width, label='Finite Difference', color='coral')
ax.set_xlabel('S parameter')
ax.set_ylabel('|Delta N| (conservation error)')
ax.set_title('Conservation of N = integral(|psi|^2 dx) after 8 time steps')
ax.set_xticks(x_pos)
ax.set_xticklabels([str(s) for s in S_values])
ax.legend()
ax.set_yscale('log')
plt.tight_layout()
plt.savefig('../plots/fig6_conservation.png', dpi=150)
plt.close()

# Export CSV
with open('../data/fig6_conservation.csv', 'w') as f:
    f.write('# Conservation quantity comparison: |N_initial - N_final| after 8 time steps\n')
    f.write('# SS: L=64, N=512, tau=0.01; FD: L=30, N=512, tau=0.001\n')
    f.write('S,delta_N_splitstep,delta_N_finitediff\n')
    for S, dss, dfd in results:
        f.write(f'{S:.8f},{dss:.8f},{dfd:.8f}\n')

print("Figure 6 done: Conservation comparison")
