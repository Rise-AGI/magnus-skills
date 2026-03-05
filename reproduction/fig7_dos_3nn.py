"""
Figure 7: Density of states with third nearest neighbor interactions.
Shows: 1st NN only (black), 3rd NN Table I (red), 3rd NN Table II (blue), Ref [1] (green).
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(".")))
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from graphene_tb import compute_dos, band_nearest, band_3nn

N_K = 800
N_BINS = 600
E_RANGE = (-10, 12)

# 1st NN only (black)
e1, d1 = compute_dos(band_nearest, dict(E2p=0.0, gamma0=-2.74, s0=0.065),
                     n_kpoints=N_K, n_bins=N_BINS, energy_range=E_RANGE)

# 3rd NN Table I (red)
e2, d2 = compute_dos(band_3nn,
                     dict(E2p=-0.21, gamma0=-2.74, gamma1=-0.07, gamma2=-0.015,
                          s0=0.065, s1=0.002, s2=0.001),
                     n_kpoints=N_K, n_bins=N_BINS, energy_range=E_RANGE)

# 3rd NN Table II (blue)
e3, d3 = compute_dos(band_3nn,
                     dict(E2p=-0.45, gamma0=-2.78, gamma1=-0.15, gamma2=-0.095,
                          s0=0.117, s1=0.004, s2=0.002),
                     n_kpoints=N_K, n_bins=N_BINS, energy_range=E_RANGE)

# Reference [1] (green)
e4, d4 = compute_dos(band_3nn,
                     dict(E2p=-0.36, gamma0=-2.78, gamma1=-0.12, gamma2=-0.068,
                          s0=0.106, s1=0.001, s2=0.003),
                     n_kpoints=N_K, n_bins=N_BINS, energy_range=E_RANGE)

header = "energy,dos_1nn,dos_3nn_tableI,dos_3nn_tableII,dos_ref"
data = np.column_stack([e1, d1, d2, d3, d4])
np.savetxt("../data/fig7_dos_3nn.csv", data, delimiter=",",
           header=header, fmt="%.8f", comments="")

fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(e1, d1, "k-", linewidth=1.0, label="1st NN")
ax.plot(e2, d2, "r-", linewidth=1.0, label="3rd NN (Table I)")
ax.plot(e3, d3, "b-", linewidth=1.0, label="3rd NN (Table II)")
ax.plot(e4, d4, "g-", linewidth=1.0, label="Ref [1]")
ax.set_xlabel("Energy in eV")
ax.set_ylabel("Density of states (arb. units)")
ax.set_xlim(E_RANGE)
ax.set_ylim(bottom=0)
ax.legend(fontsize=9)
ax.set_title("Fig 7: DOS with 3rd NN interactions")
fig.tight_layout()
fig.savefig("../plots/fig7_dos_3nn.png", dpi=150)
plt.close(fig)
print("Fig 7 done.")
