"""
Fig 14: Critical level spacing distribution p_c(s) for the 3D Anderson model.
Parameters: Box disorder W_c = 16.5, Gaussian disorder W_c = 6.1, L=10.
Band center |E| < 0.5.
Compare with Wigner surmise and Poisson distribution.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.linalg import eigvalsh
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from anderson_model import anderson_hamiltonian_3d, level_spacing_ratios, wigner_surmise

# Critical: Box disorder W_c ~ 16.5, L=10
L = 10
W_box_crit = 16.5
W_gauss_crit = 6.1
n_stat = 1000

# Box disorder
all_spacings_box = []
for r in range(n_stat):
    if r % 200 == 0:
        print(f"Box critical: realization {r}/{n_stat}")
    H = anderson_hamiltonian_3d(L, W_box_crit, disorder="box", seed=2000 + r)
    eigs = eigvalsh(H.toarray())
    s = level_spacing_ratios(eigs, energy_window=(-0.5, 0.5))
    all_spacings_box.extend(s)

all_spacings_box = np.array(all_spacings_box)
print(f"Box critical: {len(all_spacings_box)} spacings")

# Gaussian disorder
all_spacings_gauss = []
for r in range(n_stat):
    if r % 200 == 0:
        print(f"Gaussian critical: realization {r}/{n_stat}")
    H = anderson_hamiltonian_3d(L, W_gauss_crit, disorder="gaussian", seed=3000 + r)
    eigs = eigvalsh(H.toarray())
    s = level_spacing_ratios(eigs, energy_window=(-0.5, 0.5))
    all_spacings_gauss.extend(s)

all_spacings_gauss = np.array(all_spacings_gauss)
print(f"Gaussian critical: {len(all_spacings_gauss)} spacings")

# Also do L=12, Gaussian disorder for universality check
L14 = 12
all_spacings_gauss_L14 = []
n_stat_14 = 300
for r in range(n_stat_14):
    if r % 100 == 0:
        print(f"Gaussian critical L=12: realization {r}/{n_stat_14}")
    H = anderson_hamiltonian_3d(L14, W_gauss_crit, disorder="gaussian", seed=4000 + r)
    eigs = eigvalsh(H.toarray())
    s = level_spacing_ratios(eigs, energy_window=(-0.5, 0.5))
    all_spacings_gauss_L14.extend(s)

all_spacings_gauss_L14 = np.array(all_spacings_gauss_L14)

# Histograms
bins_s = np.linspace(0, 5, 100)
s_centers = 0.5 * (bins_s[:-1] + bins_s[1:])

hist_box, _ = np.histogram(all_spacings_box, bins=bins_s, density=True)
hist_gauss, _ = np.histogram(all_spacings_gauss, bins=bins_s, density=True)
hist_gauss_L14, _ = np.histogram(all_spacings_gauss_L14, bins=bins_s, density=True)

# Theoretical curves
s_th = np.linspace(0, 5, 300)
wigner = wigner_surmise(s_th, beta=1)
poisson = np.exp(-s_th)

# Save data
data = np.column_stack([s_centers, hist_box, hist_gauss, hist_gauss_L14])
np.savetxt("../data/fig14_critical_spacing.csv", data,
           delimiter=",",
           header="s,p_c_box_L10,p_c_gauss_L10,p_c_gauss_L14",
           fmt="%.8e", comments="")

# Plot
fig, ax = plt.subplots(figsize=(8, 6))
ax.plot(s_centers, hist_box, 'b-', lw=1.5, label=f"Box W={W_box_crit}, L=10")
ax.plot(s_centers, hist_gauss, 'r-', lw=1.5, label=f"Gauss W={W_gauss_crit}, L=10")
ax.plot(s_centers, hist_gauss_L14, 'go', markersize=3, label=f"Gauss W={W_gauss_crit}, L=12")
ax.plot(s_th, wigner, 'k-', lw=2, label="Wigner (GOE)")
ax.plot(s_th, poisson, 'k--', lw=2, label="Poisson")
ax.set_xlabel("s")
ax.set_ylabel(r"$p_c(s)$")
ax.set_title("Critical level spacing distribution")
ax.legend()
ax.set_xlim(0, 4)
plt.tight_layout()
plt.savefig("../plots/fig14_critical_spacing.png", dpi=150)
print("Fig 14 saved.")
