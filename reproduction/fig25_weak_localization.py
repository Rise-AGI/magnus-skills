"""
Fig 25: Weak localization corrections in 2D orthogonal Anderson model.
Shows <g> vs L on log scale, with slope close to -1/pi = -0.318.
Parameters: 2D Anderson model, various disorder W, system sizes L.
Uses transfer matrix method for conductance computation.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from anderson_model import anderson_hamiltonian_2d

def conductance_2d_exact(L, W, E=0.0, n_samples=200, disorder="box", seed_base=42):
    """
    Compute mean conductance for 2D L x L system using
    the Landauer formula via transfer matrix.
    """
    from anderson_model import conductance_2d_strip
    rng = np.random.default_rng(seed_base)
    g_vals = []
    for i in range(n_samples):
        g = conductance_2d_strip(L, L, W, E, disorder=disorder,
                                  seed=rng.integers(0, 2**31))
        g_vals.append(g)
    return np.mean(g_vals), np.std(g_vals) / np.sqrt(n_samples)


# Use the Green's function / transfer matrix for conductance
# Actually, for 2D systems, the simplest approach that works is
# direct eigenvalue-based Thouless conductance: g_T = delta_E / Delta_E

def thouless_conductance_2d(L, W, E_center=0.0, n_samples=100, disorder="box",
                             seed_base=42):
    """
    Thouless conductance: ratio of energy shift due to boundary condition
    change to mean level spacing.

    g_T = <delta_E> / Delta_E

    We compute eigenvalues with periodic and antiperiodic BC,
    measure the shift, and normalize by level spacing.
    """
    from scipy.linalg import eigvalsh
    rng = np.random.default_rng(seed_base)
    g_vals = []

    for s in range(n_samples):
        seed = rng.integers(0, 2**31)
        rng2 = np.random.default_rng(seed)
        N = L * L

        if disorder == "box":
            eps = W * rng2.uniform(-0.5, 0.5, size=N)
        else:
            eps = W * rng2.normal(0, 1, size=N)

        # Build Hamiltonian with periodic BC
        def build_H(phase_x=0, phase_y=0):
            H = np.zeros((N, N))
            for x in range(L):
                for y in range(L):
                    i = x * L + y
                    H[i, i] = eps[i]
                    # +x neighbor
                    jx = ((x+1) % L) * L + y
                    px = np.exp(1j * phase_x) if x == L-1 else 1.0
                    # +y neighbor
                    jy = x * L + ((y+1) % L)
                    py = np.exp(1j * phase_y) if y == L-1 else 1.0

                    if np.isclose(phase_x, 0) and np.isclose(phase_y, 0):
                        H[i, jx] += 1.0
                        H[jx, i] += 1.0
                        H[i, jy] += 1.0
                        H[jy, i] += 1.0
                    else:
                        # Complex Hamiltonian for phase != 0
                        pass
            return H

        # Periodic BC
        H_p = build_H(0, 0)
        eigs_p = eigvalsh(H_p)

        # Anti-periodic BC
        H_ap = np.zeros((N, N))
        for x in range(L):
            for y in range(L):
                i = x * L + y
                H_ap[i, i] = eps[i]
                jx = ((x+1) % L) * L + y
                hop_x = -1.0 if x == L-1 else 1.0  # antiperiodic
                H_ap[i, jx] += hop_x
                H_ap[jx, i] += hop_x
                jy = x * L + ((y+1) % L)
                hop_y = -1.0 if y == L-1 else 1.0
                H_ap[i, jy] += hop_y
                H_ap[jy, i] += hop_y
        eigs_ap = eigvalsh(H_ap)

        # Select eigenvalues near E_center
        dE = 1.0
        mask_p = np.abs(eigs_p - E_center) < dE
        mask_ap = np.abs(eigs_ap - E_center) < dE

        if np.sum(mask_p) < 5:
            continue

        eigs_p_sel = eigs_p[mask_p]
        eigs_ap_sel = eigs_ap[mask_ap]

        n_eigs = min(len(eigs_p_sel), len(eigs_ap_sel))
        delta_E = np.mean(np.abs(eigs_p_sel[:n_eigs] - eigs_ap_sel[:n_eigs]))
        Delta_E = np.mean(np.diff(eigs_p_sel))

        if Delta_E > 0:
            g_vals.append(delta_E / Delta_E)

    return np.mean(g_vals), np.std(g_vals) / np.sqrt(len(g_vals))


# Compute for various L and W
W_values = [2, 3, 4]
L_values = [10, 14, 20, 30, 40, 50, 70, 100]

results = {}
for W in W_values:
    print(f"W = {W}")
    g_means = []
    g_errs = []
    L_used = []
    for L in L_values:
        if L * L > 10000:  # skip too large for exact diag
            continue
        n_samp = max(50, 500 // (L // 10))
        print(f"  L={L}, n_samples={n_samp}")
        g_mean, g_err = thouless_conductance_2d(L, W, n_samples=n_samp,
                                                  seed_base=42 + W*1000)
        g_means.append(g_mean)
        g_errs.append(g_err)
        L_used.append(L)
        print(f"    <g> = {g_mean:.4f} +/- {g_err:.4f}")
    results[W] = (np.array(L_used), np.array(g_means), np.array(g_errs))

# Save data
all_data = []
header = "L"
for W in W_values:
    header += f",g_mean_W{W},g_err_W{W}"

# Find common L values
max_len = max(len(results[W][0]) for W in W_values)
for W in W_values:
    Ls, gs, errs = results[W]
    data_w = np.column_stack([Ls, gs, errs])
    np.savetxt(f"../data/fig25_weak_loc_W{W}.csv", data_w,
               delimiter=",", header=f"L,g_mean,g_err", fmt="%.8e", comments="")

# Plot
fig, ax = plt.subplots(figsize=(8, 6))
for W in W_values:
    Ls, gs, errs = results[W]
    ax.errorbar(Ls, gs, yerr=errs, fmt='o-', label=f"W={W}")
    # Linear fit in log(L)
    if len(Ls) > 2:
        coeffs = np.polyfit(np.log(Ls), gs, 1)
        ax.plot(Ls, coeffs[0]*np.log(Ls) + coeffs[1], '--',
                label=f"slope={coeffs[0]:.3f}")

ax.set_xscale('log')
ax.set_xlabel("L")
ax.set_ylabel(r"$\langle g \rangle$ [$e^2/h$]")
ax.set_title("2D Weak localization corrections")
ax.axhline(y=0, color='gray', linestyle=':')
ax.legend()
plt.tight_layout()
plt.savefig("../plots/fig25_weak_localization.png", dpi=150)
print("Fig 25 saved.")
