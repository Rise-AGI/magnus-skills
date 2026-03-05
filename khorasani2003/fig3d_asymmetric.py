"""
Figure 3d: Band structure for asymmetric ramp profile with jump.
"""
import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.realpath(sys.argv[0])))
from dtmm import n_ramp, compute_bandstructure_full

OMEGA_MAX = 1.6
N_POINTS = 800
N_LAYERS = 600

cases = [
    ("TE", 0.0, "-", "TE@theta=0"),
    ("TE", np.pi / 4, "--", "TE@theta=pi/4"),
    ("TM", 0.0, ":", "TM@theta=0"),
    ("TM", np.pi / 4, "-.", "TM@theta=pi/4"),
]

fig, ax = plt.subplots(figsize=(6, 7))
all_data = {}

for pol, theta, ls, label in cases:
    omega, kappa, cos_kL = compute_bandstructure_full(
        OMEGA_MAX, n_ramp, theta, pol, N=N_LAYERS, n_points=N_POINTS
    )
    mask = ~np.isnan(kappa)
    ax.plot(kappa[mask], omega[mask], ls, label=label, linewidth=1.2)
    key = label.replace("@", "_").replace("=", "").replace("/", "")
    all_data[f"omega_{key}"] = omega
    all_data[f"kappa_{key}"] = kappa

ax.set_xlabel("Normalized Phase kappa*L/pi")
ax.set_ylabel("Normalized Frequency omega*L/c")
ax.set_xlim(0, 1)
ax.set_ylim(0, OMEGA_MAX)
ax.legend(fontsize=8)
ax.set_title("(d) Asymmetric ramp profile")
ax.grid(True, alpha=0.3)

outdir_plots = os.path.join(os.path.dirname(sys.argv[0]), "..", "plots")
outdir_data = os.path.join(os.path.dirname(sys.argv[0]), "..", "data")
os.makedirs(outdir_plots, exist_ok=True)
os.makedirs(outdir_data, exist_ok=True)

fig.savefig(os.path.join(outdir_plots, "fig3d_asymmetric.png"), dpi=150, bbox_inches="tight")
plt.close(fig)

omega_col = all_data["omega_TE_theta0"]
cols = [omega_col]
hdr = ["omega_norm"]
for pol, theta, ls, label in cases:
    key = label.replace("@", "_").replace("=", "").replace("/", "")
    cols.append(all_data[f"kappa_{key}"])
    hdr.append(f"kappa_{key}")

arr = np.column_stack(cols)
np.savetxt(os.path.join(outdir_data, "fig3d_asymmetric.csv"), arr, delimiter=",",
           header=",".join(hdr), comments="", fmt="%.8f")
print("Fig 3d done.")
