"""
Figure 4: Split-Step soliton collision, S=0.4, v1=20, v2=-20
Reproduces paper Figure 6.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import sys
sys.path.insert(0, '.')
from nls_solver import init_two_soliton, run_splitstep

S = 0.4
v1, v2 = 20.0, -20.0
L = 64.0
N = 512
T = 1.0
tau = 0.01
x1_off, x2_off = 8.0, 18.0

psi0, x = init_two_soliton(x1_off, v1, x2_off, v2, L, N, S)
psi_ev, times = run_splitstep(psi0, S, T, tau, L, N, store_every=1)

plt.rcParams.update({'font.size': 12})
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')
X, Y = np.meshgrid(x, times)
surf = ax.plot_surface(X, Y, psi_ev, cmap=cm.jet, linewidth=0, antialiased=True, rstride=1, cstride=2)
ax.set_xlabel('x')
ax.set_ylabel('t')
ax.set_zlabel('|psi|')
ax.set_zlim(0, max(2.0, psi_ev.max() * 1.1))
ax.set_title('Split-Step: S=0.4, v1=20, v2=-20')
fig.colorbar(surf, shrink=0.5, aspect=10)
plt.tight_layout()
plt.savefig('../plots/fig4_ss_small_pos_v20.png', dpi=150)
plt.close()

with open('../data/fig4_ss_small_pos_v20.csv', 'w') as f:
    f.write('# Split-Step soliton collision: S=0.4, v1=20, v2=-20\n')
    f.write('# L=64, N=512, T=1.0, tau=0.01\n')
    t_indices = np.linspace(0, len(times)-1, 51, dtype=int)
    x_indices = np.arange(0, N, 4)
    header = 'time,' + ','.join([f'x_{x[xi]:.4f}' for xi in x_indices])
    f.write(header + '\n')
    for ti in t_indices:
        row = f'{times[ti]:.8f},' + ','.join([f'{psi_ev[ti, xi]:.8f}' for xi in x_indices])
        f.write(row + '\n')

print("Figure 4 done: Split-Step S=0.4, v=+/-20")
