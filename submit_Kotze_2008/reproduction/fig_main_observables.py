"""
Figures 7-9, 14: Main MC observables vs temperature for different lattice sizes.

Produces:
  - ../data/fig7_energy.csv
  - ../data/fig8_heat_capacity.csv
  - ../data/fig9_magnetization.csv
  - ../data/fig14_susceptibility_prime.csv
  - ../plots/fig7_energy.png
  - ../plots/fig8_heat_capacity.png
  - ../plots/fig9_magnetization.png
  - ../plots/fig14_susceptibility_prime.png
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from ising_model import run_simulation

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(SCRIPT_DIR)
DATA = os.path.join(BASE, "data")
PLOTS = os.path.join(BASE, "plots")
os.makedirs(DATA, exist_ok=True)
os.makedirs(PLOTS, exist_ok=True)

temperatures = np.arange(0.5, 5.05, 0.1)
lattice_sizes = [2, 4, 8, 16]
n_eq = 2000
n_mc_map = {2: 50000, 4: 20000, 8: 10000, 16: 5000}

all_results = {}
for L in lattice_sizes:
    print(f"Running L={L} ...")
    all_results[L] = run_simulation(L, temperatures, n_eq=n_eq,
                                     n_mc=n_mc_map[L], seed=42)

# ---- Save CSVs ----
header_base = "T," + ",".join([f"L{L}" for L in lattice_sizes])

for key, fname, ylabel in [
    ('E_avg', 'fig7_energy', 'E/N'),
    ('C', 'fig8_heat_capacity', 'C/N'),
    ('Mabs_avg', 'fig9_magnetization', '|M|/N'),
    ('chi_prime', 'fig14_susceptibility_prime', "chi'/N"),
]:
    data_arr = np.column_stack([temperatures] +
                               [all_results[L][key] for L in lattice_sizes])
    np.savetxt(os.path.join(DATA, f"{fname}.csv"), data_arr,
               delimiter=',', header=header_base, fmt='%.8f', comments='')

    fig, ax = plt.subplots(figsize=(8, 6))
    for L in lattice_sizes:
        ax.plot(temperatures, all_results[L][key], label=f'L={L}')
    ax.set_xlabel('Temperature (T)')
    ax.set_ylabel(ylabel)
    ax.legend()
    ax.set_title(fname.replace('_', ' ').title())
    fig.tight_layout()
    fig.savefig(os.path.join(PLOTS, f"{fname}.png"), dpi=150)
    plt.close(fig)
    print(f"  Saved {fname}")

print("Main observables done.")
