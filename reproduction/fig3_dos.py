"""
Fig 3: Density of states for the 3D Anderson model with box and Gaussian disorder.
Parameters: L=10, W = 0, 5, 10, 16.5 (box) and W = 0, 2, 4, 6.1 (Gaussian), V=1.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.sparse.linalg import eigsh
from scipy.linalg import eigvalsh
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from anderson_model import anderson_hamiltonian_3d, density_of_states

L = 10
n_realizations = 5
seed_base = 42

# Box disorder
W_box = [0, 5, 10, 16.5]
# Gaussian disorder
W_gauss = [0, 2, 4, 6.1]

energy_range = (-15, 15)
n_bins = 200

results_box = {}
results_gauss = {}

for W in W_box:
    all_eigs = []
    for r in range(n_realizations):
        H = anderson_hamiltonian_3d(L, W, disorder="box", seed=seed_base + r)
        eigs = eigvalsh(H.toarray())
        all_eigs.extend(eigs)
    all_eigs = np.array(all_eigs)
    centers, dos = density_of_states(all_eigs, bins=n_bins, energy_range=energy_range)
    results_box[W] = (centers, dos)
    print(f"Box W={W}: done, {len(all_eigs)} eigenvalues")

for W in W_gauss:
    all_eigs = []
    for r in range(n_realizations):
        H = anderson_hamiltonian_3d(L, W, disorder="gaussian", seed=seed_base + 100 + r)
        eigs = eigvalsh(H.toarray())
        all_eigs.extend(eigs)
    all_eigs = np.array(all_eigs)
    centers, dos = density_of_states(all_eigs, bins=n_bins, energy_range=energy_range)
    results_gauss[W] = (centers, dos)
    print(f"Gaussian W={W}: done, {len(all_eigs)} eigenvalues")

# Save data
header_box = "energy"
data_box = [results_box[W_box[0]][0]]
for W in W_box:
    header_box += f",dos_box_W{W}"
    data_box.append(results_box[W][1])
data_box = np.column_stack(data_box)
np.savetxt("../data/fig3_dos_box.csv", data_box, delimiter=",",
           header=header_box, fmt="%.8e", comments="")

header_gauss = "energy"
data_gauss = [results_gauss[W_gauss[0]][0]]
for W in W_gauss:
    header_gauss += f",dos_gauss_W{W}"
    data_gauss.append(results_gauss[W][1])
data_gauss = np.column_stack(data_gauss)
np.savetxt("../data/fig3_dos_gauss.csv", data_gauss, delimiter=",",
           header=header_gauss, fmt="%.8e", comments="")

# Plot
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

for W in W_box:
    c, d = results_box[W]
    ax1.plot(c, d, label=f"W={W}")
ax1.set_xlabel("E")
ax1.set_ylabel(r"$\rho(E)$")
ax1.set_title("Box disorder")
ax1.legend()
ax1.set_xlim(-15, 15)

for W in W_gauss:
    c, d = results_gauss[W]
    ax2.plot(c, d, label=f"W={W}")
ax2.set_xlabel("E")
ax2.set_ylabel(r"$\rho(E)$")
ax2.set_title("Gaussian disorder")
ax2.legend()
ax2.set_xlim(-15, 15)

plt.tight_layout()
plt.savefig("../plots/fig3_dos.png", dpi=150)
print("Fig 3 saved.")
