"""
Figure 5: Density of states for nearest neighbor only.
Black curve: without overlap (s0=0)
Red curve: with overlap integral (s0=0.065)
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(".")))
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from graphene_tb import compute_dos, band_nearest, band_nearest_no_overlap

N_K = 800
N_BINS = 600
E_RANGE = (-10, 12)

# Without overlap
energy_no_s, dos_no_s = compute_dos(
    band_nearest_no_overlap,
    dict(E2p=0.0, gamma0=-2.74),
    n_kpoints=N_K, n_bins=N_BINS, energy_range=E_RANGE
)

# With overlap
energy_s, dos_s = compute_dos(
    band_nearest,
    dict(E2p=0.0, gamma0=-2.74, s0=0.065),
    n_kpoints=N_K, n_bins=N_BINS, energy_range=E_RANGE
)

header = "energy,dos_no_overlap,dos_with_overlap"
data = np.column_stack([energy_no_s, dos_no_s, dos_s])
np.savetxt("../data/fig5_dos_nn.csv", data, delimiter=",",
           header=header, fmt="%.8f", comments="")

fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(energy_no_s, dos_no_s, "k-", linewidth=1.0, label="No overlap")
ax.plot(energy_s, dos_s, "r-", linewidth=1.0, label=r"With $s_0$")
ax.set_xlabel("Energy in eV")
ax.set_ylabel("Density of states (arb. units)")
ax.set_xlim(E_RANGE)
ax.set_ylim(bottom=0)
ax.legend()
ax.set_title("Fig 5: DOS (nearest neighbor)")
fig.tight_layout()
fig.savefig("../plots/fig5_dos_nn.png", dpi=150)
plt.close(fig)
print("Fig 5 done.")
