"""
Figure 2: Refractive index profiles.
Four periodic profiles varying between n=1 and n=3:
  (1) Sinusoidal  (2) Triangular  (3) Square  (4) Asymmetric ramp with jump
"""
import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(os.path.join(os.getcwd(), __file__))) if False else os.path.dirname(os.path.realpath(sys.argv[0])))
from dtmm import n_sinusoidal, n_triangular, n_square, n_ramp

L = 1.0
x = np.linspace(-L, L, 2000)  # two periods

profiles = [
    ("Sinusoidal", n_sinusoidal, "-"),
    ("Triangular", n_triangular, "--"),
    ("Square", n_square, "-."),
    ("Ramp (asymmetric)", n_ramp, ":"),
]

# Compute profiles (wrap x into [-L/2, L/2])
x_wrapped = ((x + L / 2) % L) - L / 2

fig, ax = plt.subplots(figsize=(8, 5))
data_dict = {"x_over_L": x / L}
for name, func, ls in profiles:
    n_vals = func(x_wrapped, L)
    ax.plot(x / L, n_vals, ls, label=name, linewidth=1.5)
    key = name.split()[0].lower()
    data_dict[f"n_{key}"] = n_vals

ax.set_xlabel("x/L")
ax.set_ylabel("n(x)")
ax.set_xlim(-1, 1)
ax.set_ylim(0.8, 3.2)
ax.legend()
ax.set_title("Periodic refractive index profiles")
ax.grid(True, alpha=0.3)

outdir_plots = os.path.join(os.path.dirname(sys.argv[0]), "..", "plots")
outdir_data = os.path.join(os.path.dirname(sys.argv[0]), "..", "data")
os.makedirs(outdir_plots, exist_ok=True)
os.makedirs(outdir_data, exist_ok=True)

fig.savefig(os.path.join(outdir_plots, "fig2_profiles.png"), dpi=150, bbox_inches="tight")
plt.close(fig)

# Save data
header = "x_over_L,n_sinusoidal,n_triangular,n_square,n_ramp"
arr = np.column_stack([data_dict["x_over_L"], data_dict["n_sinusoidal"],
                       data_dict["n_triangular"], data_dict["n_square"],
                       data_dict["n_ramp"]])
np.savetxt(os.path.join(outdir_data, "fig2_profiles.csv"), arr, delimiter=",",
           header=header, comments="", fmt="%.8f")
print("Fig 2 done.")
