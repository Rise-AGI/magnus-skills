"""
Figure 4: Band structure of graphene along Gamma-K-M-Gamma
with Table II parameters (freely chosen parameters).
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(".")))
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from graphene_tb import (kpath_GKMGamma, get_tick_positions,
                         band_nearest, band_2nn, band_3nn)

# Table II parameters
PARAMS_1NN = dict(E2p=0.0, gamma0=-2.74, s0=0.065)
PARAMS_2NN = dict(E2p=-0.30, gamma0=-2.77, gamma1=-0.10, s0=0.095, s1=0.003)
PARAMS_3NN = dict(E2p=-0.45, gamma0=-2.78, gamma1=-0.15, gamma2=-0.095,
                  s0=0.117, s1=0.004, s2=0.002)

# Reference [1] parameters
PARAMS_REF = dict(E2p=-0.36, gamma0=-2.78, gamma1=-0.12, gamma2=-0.068,
                  s0=0.106, s1=0.001, s2=0.003)

npts = 300
kx, ky, dist = kpath_GKMGamma(npts)
ticks = get_tick_positions(npts)

Ec_1nn, Ev_1nn = band_nearest(kx, ky, **PARAMS_1NN)
Ec_2nn, Ev_2nn = band_2nn(kx, ky, **PARAMS_2NN)
Ec_3nn, Ev_3nn = band_3nn(kx, ky, **PARAMS_3NN)
Ec_ref, Ev_ref = band_3nn(kx, ky, **PARAMS_REF)

header = "dist,Ec_1nn,Ev_1nn,Ec_2nn,Ev_2nn,Ec_3nn,Ev_3nn,Ec_ref,Ev_ref"
data = np.column_stack([dist, Ec_1nn, Ev_1nn, Ec_2nn, Ev_2nn,
                        Ec_3nn, Ev_3nn, Ec_ref, Ev_ref])
np.savetxt("../data/fig4_bands_tableII.csv", data, delimiter=",",
           header=header, fmt="%.8f", comments="")

fig, ax = plt.subplots(figsize=(8, 6))
ax.plot(dist, Ec_1nn, "k-", linewidth=1.2, label="1st NN")
ax.plot(dist, Ev_1nn, "k-", linewidth=1.2)
ax.plot(dist, Ec_2nn, "r-", linewidth=1.2, label="1st+2nd NN")
ax.plot(dist, Ev_2nn, "r-", linewidth=1.2)
ax.plot(dist, Ec_3nn, "g-", linewidth=1.2, label="1st+2nd+3rd NN")
ax.plot(dist, Ev_3nn, "g-", linewidth=1.2)
ax.plot(dist, Ec_ref, "m-", linewidth=1.0, label="Ref [1] params")
ax.plot(dist, Ev_ref, "m-", linewidth=1.0)

ax.set_ylabel("Energy in eV")
ax.set_xlabel(r"Wave vector in $\AA^{-1}$")
ax.set_xlim(dist[0], dist[-1])
ax.set_ylim(-8, 12)
ax.set_xticks(ticks)
ax.set_xticklabels([r"$\Gamma$", "K", "M", r"$\Gamma$"])
for t in ticks[1:-1]:
    ax.axvline(t, color="k", linestyle="--", linewidth=0.5)
ax.legend(loc="upper right", fontsize=9)
ax.set_title("Fig 4: Band structure (Table II)")
fig.tight_layout()
fig.savefig("../plots/fig4_bands_tableII.png", dpi=150)
plt.close(fig)
print("Fig 4 done.")
