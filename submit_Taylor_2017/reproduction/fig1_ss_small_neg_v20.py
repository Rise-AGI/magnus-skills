"""
Figure 1: Split-Step soliton collision, S=-0.1, v1=20, v2=-20
Reproduces paper Figure 2.
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

# Parameters (from paper Section 2.6)
S = -0.1
v1, v2 = 20.0, -20.0
L = 64.0
N = 512
T = 1.0
tau = 0.01
x1_off, x2_off = 8.0, 18.0

# Initialize and run
psi0, x = init_two_soliton(x1_off, v1, x2_off, v2, L, N, S)
psi_ev, times = run_splitstep(psi0, S, T, tau, L, N, store_every=1)

# Plot 3D surface
plt.rcParams.update({'font.size': 12})
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')
X, Y = np.meshgrid(x, times)
surf = ax.plot_surface(X, Y, psi_ev, cmap=cm.jet, linewidth=0, antialiased=True, rstride=1, cstride=2)
ax.set_xlabel('x')
ax.set_ylabel('t')
ax.set_zlabel('|psi|')
ax.set_zlim(0, max(2.0, psi_ev.max() * 1.1))
ax.set_title('Split-Step: S=-0.1, v1=20, v2=-20')
fig.colorbar(surf, shrink=0.5, aspect=10)
plt.tight_layout()
plt.savefig('../plots/fig1_ss_small_neg_v20.png', dpi=150)
plt.close()

# Export CSV
with open('../data/fig1_ss_small_neg_v20.csv', 'w') as f:
    f.write('# Split-Step soliton collision: S=-0.1, v1=20, v2=-20\n')
    f.write('# L=64, N=512, T=1.0, tau=0.01\n')
    # Store a subset of time steps for manageable file size
    t_indices = np.linspace(0, len(times)-1, 51, dtype=int)
    x_indices = np.arange(0, N, 4)  # every 4th point
    header = 'time,' + ','.join([f'x_{x[xi]:.4f}' for xi in x_indices])
    f.write(header + '\n')
    for ti in t_indices:
        row = f'{times[ti]:.8f},' + ','.join([f'{psi_ev[ti, xi]:.8f}' for xi in x_indices])
        f.write(row + '\n')

print("Figure 1 done: Split-Step S=-0.1, v=+/-20")
