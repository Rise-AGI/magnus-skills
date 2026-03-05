"""
Figure 6: Density of states with second nearest neighbor interactions.
Shows multiple curves: 1st NN only (black), 1st+2nd NN with Table I (red),
1st+2nd NN with Table II (blue), and Ref [1] (green).
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(".")))
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from graphene_tb import compute_dos, band_nearest, band_2nn, band_3nn

N_K = 800
N_BINS = 600
E_RANGE = (-10, 12)

# 1st NN only (black)
e1, d1 = compute_dos(band_nearest, dict(E2p=0.0, gamma0=-2.74, s0=0.065),
                     n_kpoints=N_K, n_bins=N_BINS, energy_range=E_RANGE)

# 1st+2nd NN Table I (red)
e2, d2 = compute_dos(band_2nn,
                     dict(E2p=-0.21, gamma0=-2.74, gamma1=-0.07, s0=0.065, s1=0.002),
                     n_kpoints=N_K, n_bins=N_BINS, energy_range=E_RANGE)

# 1st+2nd NN Table II (blue)
e3, d3 = compute_dos(band_2nn,
                     dict(E2p=-0.30, gamma0=-2.77, gamma1=-0.10, s0=0.095, s1=0.003),
                     n_kpoints=N_K, n_bins=N_BINS, energy_range=E_RANGE)

# Reference [1] with 3rd NN (green)
e4, d4 = compute_dos(band_3nn,
                     dict(E2p=-0.36, gamma0=-2.78, gamma1=-0.12, gamma2=-0.068,
                          s0=0.106, s1=0.001, s2=0.003),
                     n_kpoints=N_K, n_bins=N_BINS, energy_range=E_RANGE)

header = "energy,dos_1nn,dos_2nn_tableI,dos_2nn_tableII,dos_ref"
data = np.column_stack([e1, d1, d2, d3, d4])
np.savetxt("../data/fig6_dos_2nn.csv", data, delimiter=",",
           header=header, fmt="%.8f", comments="")

fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(e1, d1, "k-", linewidth=1.0, label="1st NN")
ax.plot(e2, d2, "r-", linewidth=1.0, label="1st+2nd NN (Table I)")
ax.plot(e3, d3, "b-", linewidth=1.0, label="1st+2nd NN (Table II)")
ax.plot(e4, d4, "g-", linewidth=1.0, label="Ref [1]")
ax.set_xlabel("Energy in eV")
ax.set_ylabel("Density of states (arb. units)")
ax.set_xlim(E_RANGE)
ax.set_ylim(bottom=0)
ax.legend(fontsize=9)
ax.set_title("Fig 6: DOS with 2nd NN interactions")
fig.tight_layout()
fig.savefig("../plots/fig6_dos_2nn.png", dpi=150)
plt.close(fig)
print("Fig 6 done.")
