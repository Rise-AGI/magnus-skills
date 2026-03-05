"""
Figure 5: Random center vortex model.

Demonstrates that the Z_2 random center vortex model produces
Wilson loops with area-law behavior, providing a microscopic
mechanism for confinement.

W(C) = (1 - 2f)^A, where f is the vortex piercing probability
and A is the minimal area of the loop.

Related to Greensite (2003) Section 5 - the center vortex confinement mechanism.
"""
import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) if os.path.dirname(os.path.abspath(__file__)) else '.')
from center_vortex import vortex_wilson_loop_analytic, vortex_string_tension, simulate_vortex_model_2d

# Parameters
L = 16
n_configs = 500
rng = np.random.default_rng(42)
f_values = [0.05, 0.10, 0.15, 0.20, 0.25]

print("Figure 5: Random center vortex model")
print(f"2D lattice: {L}x{L}, configurations: {n_configs}")

all_results = {}

for f in f_values:
    print(f"  f = {f:.2f} ... ", end="", flush=True)
    mc_results = simulate_vortex_model_2d(L, L, f, n_configs, rng)
    all_results[f] = mc_results
    sigma = vortex_string_tension(f)
    print(f"sigma = {sigma:.6f}")
    for S, w in sorted(mc_results.items()):
        w_anal = vortex_wilson_loop_analytic(S*S, f)
        print(f"    W({S}x{S}): MC = {w:.6f}, Analytic = {w_anal:.6f}")

# Save data
data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data')
os.makedirs(data_dir, exist_ok=True)

with open(os.path.join(data_dir, 'fig5_vortex_model.csv'), 'w') as f_out:
    f_out.write("# Figure 5: Random center vortex model Wilson loops\n")
    f_out.write("# 2D lattice L=16, n_configs=500\n")
    f_out.write("f,loop_size,area,wilson_mc,wilson_analytic,string_tension\n")
    for f in f_values:
        sigma = vortex_string_tension(f)
        for S in sorted(all_results[f].keys()):
            w_mc = all_results[f][S]
            w_anal = vortex_wilson_loop_analytic(S*S, f)
            f_out.write(f"{f:.2f},{S},{S*S},{w_mc:.8f},{w_anal:.8f},{sigma:.8f}\n")

# Plot
plot_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'plots')
os.makedirs(plot_dir, exist_ok=True)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# Left panel: Wilson loops vs area for different f
markers = ['o', 's', '^', 'D', 'v']
colors = ['blue', 'red', 'green', 'purple', 'orange']

for i, f in enumerate(f_values):
    areas = []
    mc_vals = []
    anal_vals = []
    for S in sorted(all_results[f].keys()):
        areas.append(S * S)
        mc_vals.append(all_results[f][S])
        anal_vals.append(vortex_wilson_loop_analytic(S*S, f))
    ax1.plot(areas, mc_vals, markers[i], color=colors[i], markersize=6,
             label=f'MC f={f:.2f}')
    area_fine = np.linspace(0, max(areas), 100)
    ax1.plot(area_fine, vortex_wilson_loop_analytic(area_fine, f),
             '-', color=colors[i], linewidth=1, alpha=0.7)

ax1.set_yscale('log')
ax1.set_xlabel(r'Loop area $A = S^2$', fontsize=13)
ax1.set_ylabel(r'$W(S \times S)$', fontsize=13)
ax1.set_title('Vortex Model: Wilson Loops', fontsize=13)
ax1.legend(fontsize=10)
ax1.grid(True, alpha=0.3)

# Right panel: String tension vs f
f_fine = np.linspace(0.01, 0.45, 100)
sigma_anal = np.array([vortex_string_tension(ff) for ff in f_fine])

# Extract MC string tensions from -ln(W(1x1))
sigma_mc = []
for f in f_values:
    if 1 in all_results[f] and all_results[f][1] > 0:
        sigma_mc.append(-np.log(all_results[f][1]))
    else:
        sigma_mc.append(float('nan'))

ax2.plot(f_fine, sigma_anal, 'k-', linewidth=2, label=r'$-\ln(1-2f)$')
ax2.plot(f_values, sigma_mc, 'ro', markersize=8, label=r'MC: $-\ln W(1\times1)$')
ax2.set_xlabel('Vortex density $f$', fontsize=13)
ax2.set_ylabel(r'String tension $\sigma$', fontsize=13)
ax2.set_title('Vortex Model: String Tension', fontsize=13)
ax2.legend(fontsize=12)
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(plot_dir, 'fig5_vortex_model.png'), dpi=150)
print("Saved fig5_vortex_model.png")
