"""
Figure 5: Frontiers between instability domains in (gamma_b, alpha) space.
Computed for various OmegaB values and R=1/1836, R=1/100.

Uses analytical growth rate expressions from Tables 1 & 2 to determine
which mode dominates at each point, then plots domain boundaries.

Reproduces Fig. 5 of Bret (2009).
"""
import sys
import os
import numpy as np

sys.path.insert(0, os.path.dirname(__file__) if "__file__" in dir() else ".")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def analytical_growth_rates_dilute(alpha, gamma_b, OmegaB, R):
    """Compute analytical growth rates for dilute beam regime (Tables 1 & 2).

    Returns dict of mode_name -> growth_rate.
    """
    rates = {}

    # Two-Stream (Eq. 10): delta = sqrt(3)/2^(4/3) * alpha^(1/3) / gamma_b
    rates["TwoStream"] = np.sqrt(3) / 2**(4./3) * alpha**(1./3) / gamma_b

    # Buneman (Eq. 11): delta = sqrt(3)/2^(4/3) * R^(1/3)
    rates["Buneman"] = np.sqrt(3) / 2**(4./3) * R**(1./3)

    # Filamentation (Zx -> inf, Zz=0, OmegaB=0):
    # delta ~ sqrt(alpha) * beta / gamma_b  (from Bret et al. 2004)
    beta = np.sqrt(1 - 1.0/gamma_b**2) if gamma_b > 1 else 0.0
    if OmegaB == 0 and beta > 0:
        rates["Filamentation"] = np.sqrt(alpha) * beta / gamma_b
    elif OmegaB > 0 and OmegaB < np.sqrt(gamma_b):
        # Filamentation stabilized for OmegaB > sqrt(gamma_b)
        rates["Filamentation"] = max(0, np.sqrt(alpha) * beta / gamma_b *
                                     (1 - OmegaB**2 / gamma_b))
    else:
        rates["Filamentation"] = 0.0

    # Oblique modes (OmegaB=0, Table 1):
    # delta = sqrt(3)/2^(4/3) * (alpha/gamma_b)^(1/3)
    if OmegaB < 1:
        rates["Oblique"] = np.sqrt(3) / 2**(4./3) * (alpha / gamma_b)**(1./3)
    else:
        # For OmegaB > 1 (Table 1): delta ~ alpha^(1/3) / gamma_b
        rates["Oblique"] = 0.25 * alpha**(1./3) / gamma_b

    # Upper-Hybrid-Like (OmegaB > 1, Table 1):
    if OmegaB > 1:
        rates["UHL"] = 0.5 * np.sqrt(alpha) / OmegaB
    else:
        rates["UHL"] = 0.0

    return rates


def dominant_mode(alpha, gamma_b, OmegaB, R):
    """Return the name of the dominant mode."""
    rates = analytical_growth_rates_dilute(alpha, gamma_b, OmegaB, R)
    if not rates:
        return "None"
    return max(rates, key=rates.get)


# Compute hierarchy maps
os.makedirs("../data", exist_ok=True)
os.makedirs("../plots", exist_ok=True)

alpha_arr = np.logspace(-1, 0, 80)  # 0.1 to 1
gamma_b_arr = np.logspace(0, 3, 100)  # 1 to 1000
OmegaB_values = [0, 1, 3, 10]
R_values = [1.0/1836, 1.0/100]

mode_colors = {
    "TwoStream": "blue",
    "Buneman": "red",
    "Filamentation": "green",
    "Oblique": "orange",
    "UHL": "purple",
    "None": "white"
}
mode_ids = {"TwoStream": 0, "Buneman": 1, "Filamentation": 2, "Oblique": 3, "UHL": 4, "None": 5}

fig, axes = plt.subplots(len(R_values), len(OmegaB_values),
                         figsize=(16, 8), squeeze=False)

for ri, R in enumerate(R_values):
    for oi, OmegaB in enumerate(OmegaB_values):
        ax = axes[ri, oi]
        dom_map = np.zeros((len(gamma_b_arr), len(alpha_arr)))

        for i, gb in enumerate(gamma_b_arr):
            for j, al in enumerate(alpha_arr):
                mode = dominant_mode(al, gb, OmegaB, R)
                dom_map[i, j] = mode_ids.get(mode, 5)

        im = ax.pcolormesh(alpha_arr, gamma_b_arr, dom_map, cmap="Set2",
                           vmin=0, vmax=5, shading="auto")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(r"$\alpha = n_b/n_p$")
        ax.set_ylabel(r"$\gamma_b$")
        R_label = "1/1836" if abs(R - 1/1836) < 1e-6 else "1/100"
        ax.set_title(rf"$\Omega_B={OmegaB},\ R={R_label}$", fontsize=9)

fig.suptitle("Fig 5: Instability hierarchy in parameter space", fontsize=14)
fig.tight_layout()
fig.savefig("../plots/fig5_hierarchy.png", dpi=150)
plt.close()

# Save hierarchy data
for R in R_values:
    R_label = "R1836" if abs(R - 1/1836) < 1e-6 else "R100"
    for OmegaB in OmegaB_values:
        dom_arr = []
        for gb in gamma_b_arr:
            for al in alpha_arr:
                mode = dominant_mode(al, gb, OmegaB, R)
                rates = analytical_growth_rates_dilute(al, gb, OmegaB, R)
                max_rate = max(rates.values())
                dom_arr.append([al, gb, mode_ids.get(mode, 5), max_rate])
        data = np.array(dom_arr)
        np.savetxt(f"../data/fig5_hierarchy_{R_label}_OB{OmegaB}.csv",
                   data, delimiter=",",
                   header="alpha,gamma_b,dominant_mode_id,max_growth_rate",
                   fmt="%.6e")

print("Done. Saved fig5 data and plots.")
print("Mode IDs: TwoStream=0, Buneman=1, Filamentation=2, Oblique=3, UHL=4")
