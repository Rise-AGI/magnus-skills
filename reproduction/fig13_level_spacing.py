"""
Fig 13: Level spacing distribution p(s) for the 3D Anderson model.
Left: metallic regime (Gaussian W=2), compared with Wigner surmise.
Right: localized regime (Box W=32), compared with Poisson.
Parameters: L=14 for metallic, L=10 for localized. Band center |E|<0.5.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.linalg import eigvalsh
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from anderson_model import anderson_hamiltonian_3d, level_spacing_ratios, wigner_surmise

# Metallic: Gaussian W=2, L=10, Nstat=500
L_met = 10
W_met = 2.0
n_stat_met = 500

all_spacings_met = []
for r in range(n_stat_met):
    if r % 50 == 0:
        print(f"Metallic realization {r}/{n_stat_met}")
    H = anderson_hamiltonian_3d(L_met, W_met, disorder="gaussian", seed=42 + r)
    eigs = eigvalsh(H.toarray())
    s = level_spacing_ratios(eigs, energy_window=(-0.5, 0.5))
    all_spacings_met.extend(s)

all_spacings_met = np.array(all_spacings_met)
print(f"Metallic: {len(all_spacings_met)} spacings")

# Localized: Box W=32, L=10
L_loc = 10
W_loc = 32.0
n_stat_loc = 500

all_spacings_loc = []
for r in range(n_stat_loc):
    if r % 100 == 0:
        print(f"Localized realization {r}/{n_stat_loc}")
    H = anderson_hamiltonian_3d(L_loc, W_loc, disorder="box", seed=1000 + r)
    eigs = eigvalsh(H.toarray())
    s = level_spacing_ratios(eigs, energy_window=(-0.5, 0.5))
    all_spacings_loc.extend(s)

all_spacings_loc = np.array(all_spacings_loc)
print(f"Localized: {len(all_spacings_loc)} spacings")

# Histograms
bins_s = np.linspace(0, 5, 100)
s_centers = 0.5 * (bins_s[:-1] + bins_s[1:])

hist_met, _ = np.histogram(all_spacings_met, bins=bins_s, density=True)
hist_loc, _ = np.histogram(all_spacings_loc, bins=bins_s, density=True)

# Theoretical curves
s_th = np.linspace(0, 5, 300)
wigner = wigner_surmise(s_th, beta=1)
poisson = np.exp(-s_th)

# Save data
data_met = np.column_stack([s_centers, hist_met])
np.savetxt("../data/fig13_level_spacing_metallic.csv", data_met,
           delimiter=",", header="s,p_s_metallic", fmt="%.8e", comments="")

data_loc = np.column_stack([s_centers, hist_loc])
np.savetxt("../data/fig13_level_spacing_localized.csv", data_loc,
           delimiter=",", header="s,p_s_localized", fmt="%.8e", comments="")

data_th = np.column_stack([s_th, wigner, poisson])
np.savetxt("../data/fig13_level_spacing_theory.csv", data_th,
           delimiter=",", header="s,wigner_GOE,poisson", fmt="%.8e", comments="")

# Plot
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

ax1.bar(s_centers, hist_met, width=s_centers[1]-s_centers[0], alpha=0.6, label=f"Gaussian W={W_met}, L={L_met}")
ax1.plot(s_th, wigner, 'k-', lw=2, label="Wigner (GOE)")
ax1.set_xlabel("s")
ax1.set_ylabel("p(s)")
ax1.set_title("Metallic regime")
ax1.legend()
ax1.set_xlim(0, 4)

ax2.bar(s_centers, hist_loc, width=s_centers[1]-s_centers[0], alpha=0.6, color='orange',
        label=f"Box W={W_loc}, L={L_loc}")
ax2.plot(s_th, poisson, 'k--', lw=2, label="Poisson")
ax2.set_xlabel("s")
ax2.set_ylabel("p(s)")
ax2.set_title("Localized regime")
ax2.legend()
ax2.set_xlim(0, 4)

plt.tight_layout()
plt.savefig("../plots/fig13_level_spacing.png", dpi=150)
print("Fig 13 saved.")
