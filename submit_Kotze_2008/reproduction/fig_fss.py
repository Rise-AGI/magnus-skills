"""
Figures 19, 20, 21: Finite size scaling, reduced-unit plot, and cumulant intersection.

Produces:
  - ../data/fig19_critical_exponents.csv
  - ../data/fig20_reduced_unit.csv
  - ../data/fig21_cumulants.csv
  - ../plots/fig19_critical_exponents.png
  - ../plots/fig20_reduced_unit.png
  - ../plots/fig21_cumulants.png
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

Tc_exact = 2.0 / np.log(1 + np.sqrt(2))  # ~2.269

temps_fine = np.arange(1.5, 3.01, 0.02)
temps_coarse = np.arange(0.5, 5.05, 0.1)

lattice_sizes = [4, 8, 16, 32]
n_eq = 3000
n_mc_map = {4: 30000, 8: 15000, 16: 8000, 32: 4000}

# ---------- Coarse grid simulations ----------
print("=== Coarse grid simulations ===")
coarse_results = {}
for L in lattice_sizes:
    print(f"Running L={L} (coarse)...")
    coarse_results[L] = run_simulation(L, temps_coarse, n_eq=n_eq,
                                        n_mc=n_mc_map[L], seed=42)

# ---------- Fig 19: Critical exponents ----------
print("\n=== Fig 19: Critical exponents ===")
Ls = np.array(lattice_sizes, dtype=float)
chi_peaks = []
mabs_at_Tc = []
C_peaks = []

for L in lattice_sizes:
    r = coarse_results[L]
    idx_chi = np.argmax(r['chi_prime'])
    chi_peaks.append(r['chi_prime'][idx_chi])
    mabs_at_Tc.append(r['Mabs_avg'][idx_chi])
    C_peaks.append(r['C'][idx_chi])

chi_peaks = np.array(chi_peaks)
mabs_at_Tc = np.array(mabs_at_Tc)
C_peaks = np.array(C_peaks)

ln_L = np.log(Ls)
fit_gamma = np.polyfit(ln_L, np.log(chi_peaks), 1)
gamma_nu = fit_gamma[0]
fit_beta = np.polyfit(ln_L, np.log(mabs_at_Tc), 1)
beta_nu = -fit_beta[0]

data19 = np.column_stack([Ls, chi_peaks, mabs_at_Tc, C_peaks])
np.savetxt(os.path.join(DATA, "fig19_critical_exponents.csv"), data19,
           delimiter=',', header="L,chi_peak,M_at_Tc,C_peak", fmt='%.8f',
           comments='')

fig, axes = plt.subplots(1, 3, figsize=(15, 5))
axes[0].plot(ln_L, np.log(mabs_at_Tc), 'o-')
axes[0].set_xlabel('ln(L)')
axes[0].set_ylabel('ln(|M| at Tc)')
axes[0].set_title(f'beta/nu = {beta_nu:.3f} (exact: 0.125)')

axes[1].plot(ln_L, np.log(chi_peaks), 'o-')
axes[1].set_xlabel('ln(L)')
axes[1].set_ylabel("ln(chi' peak)")
axes[1].set_title(f'gamma/nu = {gamma_nu:.3f} (exact: 1.75)')

axes[2].plot(ln_L, np.log(C_peaks), 'o-')
axes[2].set_xlabel('ln(L)')
axes[2].set_ylabel('ln(C peak)')
axes[2].set_title('C ~ ln(L) (alpha=0)')

fig.tight_layout()
fig.savefig(os.path.join(PLOTS, "fig19_critical_exponents.png"), dpi=150)
plt.close(fig)
print(f"  gamma/nu = {gamma_nu:.3f}, beta/nu = {beta_nu:.3f}")

# ---------- Fig 20: Reduced unit plot ----------
print("\n=== Fig 20: Reduced unit plot ===")
fig, ax = plt.subplots(figsize=(8, 6))
rows_20 = []
for L in lattice_sizes:
    r = coarse_results[L]
    t_reduced = (r['T'] - Tc_exact) / Tc_exact * L
    chi_scaled = r['chi_prime'] / L**(1.75)
    ax.plot(t_reduced, chi_scaled, 'o-', ms=3, label=f'L={L}')
    for i in range(len(r['T'])):
        rows_20.append([L, r['T'][i], t_reduced[i], chi_scaled[i]])

ax.set_xlabel('(T - Tc) / Tc * L')
ax.set_ylabel("chi' / L^(7/4)")
ax.set_xlim(-10, 10)
ax.legend()
ax.set_title('Reduced unit plot for susceptibility')
fig.tight_layout()
fig.savefig(os.path.join(PLOTS, "fig20_reduced_unit.png"), dpi=150)
plt.close(fig)

data20 = np.array(rows_20)
np.savetxt(os.path.join(DATA, "fig20_reduced_unit.csv"), data20,
           delimiter=',', header="L,T,t_reduced,chi_scaled", fmt='%.8f',
           comments='')

# ---------- Fig 21: Cumulant intersection ----------
print("\n=== Fig 21: Cumulant intersection ===")
print("Running fine-grid simulations around Tc...")
fine_results = {}
for L in lattice_sizes:
    print(f"  L={L} (fine grid)...")
    fine_results[L] = run_simulation(L, temps_fine, n_eq=n_eq,
                                      n_mc=n_mc_map[L], seed=123)

fig, ax = plt.subplots(figsize=(8, 6))
rows_21 = []
for L in lattice_sizes:
    r = fine_results[L]
    ax.plot(r['T'], r['U_L'], 'o-', ms=3, label=f'L={L}')
    for i in range(len(r['T'])):
        rows_21.append([L, r['T'][i], r['U_L'][i]])

ax.axvline(Tc_exact, color='k', ls='--', alpha=0.5, label=f'Tc={Tc_exact:.3f}')
ax.set_xlabel('Temperature (T)')
ax.set_ylabel('U_L (Binder cumulant)')
ax.legend()
ax.set_title('Cumulant intersection')
fig.tight_layout()
fig.savefig(os.path.join(PLOTS, "fig21_cumulants.png"), dpi=150)
plt.close(fig)

data21 = np.array(rows_21)
np.savetxt(os.path.join(DATA, "fig21_cumulants.csv"), data21,
           delimiter=',', header="L,T,U_L", fmt='%.8f', comments='')

print("\nFinite size scaling done.")
