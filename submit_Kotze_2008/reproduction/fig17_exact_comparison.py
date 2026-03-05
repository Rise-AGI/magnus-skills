"""
Figure 17: Exact vs MC comparison for 2x2 lattice.

Produces:
  - ../data/fig17_exact_comparison.csv
  - ../plots/fig17_exact_comparison.png
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from ising_model import run_simulation, exact_2x2

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(SCRIPT_DIR)
DATA = os.path.join(BASE, "data")
PLOTS = os.path.join(BASE, "plots")
os.makedirs(DATA, exist_ok=True)
os.makedirs(PLOTS, exist_ok=True)

temperatures = np.arange(0.5, 5.05, 0.1)

print("Running L=2 MC simulation...")
mc = run_simulation(2, temperatures, n_eq=5000, n_mc=100000, seed=42)
ex = exact_2x2(temperatures)

data_arr = np.column_stack([
    temperatures,
    mc['C'], ex['C'],
    mc['chi_prime'], ex['chi_prime'],
])
header = "T,C_mc,C_exact,chi_prime_mc,chi_prime_exact"
np.savetxt(os.path.join(DATA, "fig17_exact_comparison.csv"), data_arr,
           delimiter=',', header=header, fmt='%.8f', comments='')

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

ax1.plot(temperatures, mc['C'], 'o', ms=4, label='MC')
ax1.plot(temperatures, ex['C'], '-', label='Exact')
ax1.set_xlabel('Temperature (T)')
ax1.set_ylabel('C/N')
ax1.legend()
ax1.set_title('Heat Capacity (L=2)')

ax2.plot(temperatures, mc['chi_prime'], 'o', ms=4, label='MC')
ax2.plot(temperatures, ex['chi_prime'], '-', label='Exact')
ax2.set_xlabel('Temperature (T)')
ax2.set_ylabel("chi'/N")
ax2.legend()
ax2.set_title('Susceptibility (L=2)')

fig.tight_layout()
fig.savefig(os.path.join(PLOTS, "fig17_exact_comparison.png"), dpi=150)
plt.close(fig)
print("Fig 17 done.")
