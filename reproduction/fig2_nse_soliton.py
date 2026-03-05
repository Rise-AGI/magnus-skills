"""
Figure 2: NSE soliton error convergence.

Reproduces the RMS error vs stepsize plot for different propagation schemes
(SiSSM, SySSM, IFM, LEM) applied to the single-soliton solution of the
nonlinear Schroedinger equation (NSE).

Reference: Melchert & Demircan, CPC (2022), Section 9.1, Figure 2.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.fft import fft, ifft, fftfreq
import csv
import os

# --- NSE soliton parameters ---
# The NSE: i dE/dz + 0.5 * d2E/dt2 + |E|^2 E = 0
# Exact soliton: E(z,t) = sech(t) * exp(i*z/2)

# Domain parameters (matching typical py-fmas setup)
t_max = 25.0       # half-period
t_num = 2**12       # number of temporal grid points
z_max = 1.0         # propagation distance (in soliton periods ~ pi/2)

# Temporal grid
dt = 2 * t_max / t_num
t = np.linspace(-t_max, t_max - dt, t_num)

# Frequency grid
w = 2 * np.pi * fftfreq(t_num, d=dt)

# Linear operator for NSE: i*u_z + (1/2)*u_tt + |u|^2*u = 0
# => u_z = -i*w^2/2 * u_w + i*FT(|u|^2*u)
# So L_hat = -i*w^2/2
L_hat = -1j * w**2 / 2.0

# Exact soliton solution
def exact_soliton(z, t):
    return 1.0 / np.cosh(t) * np.exp(1j * z / 2.0)

# Nonlinear operator for NSE: N(E) = i |E|^2 E
def N_func(E):
    return 1j * np.abs(E)**2 * E

def N_hat_func(E):
    return fft(N_func(E))

# --- Propagation methods ---

def propagate_SiSSM(E0_w, z_max, z_steps, L_hat, rk_order=2):
    """Simple split-step Fourier method."""
    dz = z_max / z_steps
    E_w = E0_w.copy()
    for _ in range(z_steps):
        # Nonlinear step (RK2 or RK4)
        E_t = ifft(E_w)
        if rk_order == 2:
            k1 = dz * N_hat_func(E_t)
            E_mid = ifft(E_w + 0.5 * k1)
            k2 = dz * N_hat_func(E_mid)
            E_w = E_w + k2
        else:
            k1 = dz * N_hat_func(E_t)
            E2 = ifft(E_w + 0.5 * k1)
            k2 = dz * N_hat_func(E2)
            E3 = ifft(E_w + 0.5 * k2)
            k3 = dz * N_hat_func(E3)
            E4 = ifft(E_w + k3)
            k4 = dz * N_hat_func(E4)
            E_w = E_w + (k1 + 2*k2 + 2*k3 + k4) / 6.0
        # Linear step
        E_w = E_w * np.exp(L_hat * dz)
    return E_w

def propagate_SySSM(E0_w, z_max, z_steps, L_hat, rk_order=2):
    """Symmetric split-step Fourier method."""
    dz = z_max / z_steps
    E_w = E0_w.copy()
    for _ in range(z_steps):
        # Half linear step
        E_w = E_w * np.exp(L_hat * dz / 2.0)
        # Full nonlinear step
        E_t = ifft(E_w)
        if rk_order == 2:
            k1 = dz * N_hat_func(E_t)
            E_mid = ifft(E_w + 0.5 * k1)
            k2 = dz * N_hat_func(E_mid)
            E_w = E_w + k2
        else:
            k1 = dz * N_hat_func(E_t)
            E2 = ifft(E_w + 0.5 * k1)
            k2 = dz * N_hat_func(E2)
            E3 = ifft(E_w + 0.5 * k2)
            k3 = dz * N_hat_func(E3)
            E4 = ifft(E_w + k3)
            k4 = dz * N_hat_func(E4)
            E_w = E_w + (k1 + 2*k2 + 2*k3 + k4) / 6.0
        # Half linear step
        E_w = E_w * np.exp(L_hat * dz / 2.0)
    return E_w

def propagate_IFM_RK4IP(E0_w, z_max, z_steps, L_hat):
    """Integrating factor method with RK4 in interaction picture (Eq. 18-19).

    Reference point z0 = z + dz/2 (midpoint). The IP field is:
    psi_w(z) = exp(-L*(z - z0)) * E_w(z)
    At z_n:   psi_n = exp(L*dz/2) * E_w(z_n)
    At z_n+h: A_{n+1} = exp(L*dz/2) * psi_{n+1}
    """
    dz = z_max / z_steps
    E_w = E0_w.copy()
    half = np.exp(L_hat * dz / 2.0)
    inv_half = np.exp(-L_hat * dz / 2.0)

    for _ in range(z_steps):
        # Transform to IP
        psi_w = half * E_w

        # k1 at z_n (f = -dz/2): G = exp(L*dz/2) * N_hat(exp(-L*dz/2)*psi)
        # exp(-L*dz/2)*psi_n = exp(-L*dz/2)*exp(L*dz/2)*E_w = E_w
        k1 = dz * half * N_hat_func(ifft(E_w))

        # k2 at z_n+dz/2 (f = 0): G = N_hat(psi + k1/2)
        k2 = dz * N_hat_func(ifft(psi_w + 0.5 * k1))

        # k3 at z_n+dz/2 (f = 0): G = N_hat(psi + k2/2)
        k3 = dz * N_hat_func(ifft(psi_w + 0.5 * k2))

        # k4 at z_n+dz (f = dz/2): G = exp(-L*dz/2) * N_hat(exp(L*dz/2)*(psi+k3))
        k4 = dz * inv_half * N_hat_func(ifft(half * (psi_w + k3)))

        psi_w = psi_w + (k1 + 2*k2 + 2*k3 + k4) / 6.0

        # Transform back: A_{n+1} = exp(L*dz/2) * psi_{n+1}
        E_w = half * psi_w

    return E_w

def propagate_LEM(E0_w, z_max, z_steps, L_hat, goal_error=1e-7):
    """Local error method with adaptive stepsize."""
    dz = z_max / z_steps
    E_w = E0_w.copy()
    z_current = 0.0

    for _ in range(z_steps):
        z_target = z_current + dz
        h = dz  # initial substep size
        E_w_local = E_w.copy()
        z_local = z_current

        while z_local < z_target - 1e-15:
            h = min(h, z_target - z_local)

            # Take one full step of size h
            E_full = _ifm_rk4_step(E_w_local, h, L_hat)
            # Take two half steps
            E_half1 = _ifm_rk4_step(E_w_local, h/2, L_hat)
            E_half2 = _ifm_rk4_step(E_half1, h/2, L_hat)

            # Estimate local error
            err = np.sqrt(np.mean(np.abs(E_half2 - E_full)**2)) / np.sqrt(np.mean(np.abs(E_half2)**2))

            if err > 2 * goal_error:
                h = h / 2
                continue
            elif err > goal_error:
                E_w_local = E_half2
                z_local += h
                h = h * 2**(-1/3)
            elif err < goal_error / 2:
                E_w_local = E_half2
                z_local += h
                h = h * 2**(1/3)
            else:
                E_w_local = E_half2
                z_local += h

        E_w = E_w_local
        z_current = z_target

    return E_w

def _ifm_rk4_step(E_w, h, L_hat):
    """Single IFM-RK4IP step of size h."""
    half = np.exp(L_hat * h / 2.0)
    inv_half = np.exp(-L_hat * h / 2.0)

    psi_w = half * E_w

    k1 = h * half * N_hat_func(ifft(E_w))
    k2 = h * N_hat_func(ifft(psi_w + 0.5 * k1))
    k3 = h * N_hat_func(ifft(psi_w + 0.5 * k2))
    k4 = h * inv_half * N_hat_func(ifft(half * (psi_w + k3)))

    psi_w = psi_w + (k1 + 2*k2 + 2*k3 + k4) / 6.0

    return half * psi_w


def compute_rms_error(E_numerical_w, E_exact_w):
    """Compute RMS error between numerical and exact solution."""
    E_num = ifft(E_numerical_w)
    E_ex = ifft(E_exact_w)
    return np.sqrt(np.mean(np.abs(E_num - E_ex)**2))


# --- Main computation ---
print("Computing NSE soliton error convergence (Figure 2)...")

# Initial condition
E0 = exact_soliton(0, t)
E0_w = fft(E0)

# Exact solution at z = z_max
E_exact = exact_soliton(z_max, t)
E_exact_w = fft(E_exact)

# Stepsize values to test
z_steps_list = [2**n for n in range(3, 14)]
dz_list = [z_max / ns for ns in z_steps_list]

methods = {
    "SiSSM": lambda E0w, zmax, zsteps: propagate_SiSSM(E0w, zmax, zsteps, L_hat, rk_order=2),
    "SySSM": lambda E0w, zmax, zsteps: propagate_SySSM(E0w, zmax, zsteps, L_hat, rk_order=2),
    "IFM-RK4IP": lambda E0w, zmax, zsteps: propagate_IFM_RK4IP(E0w, zmax, zsteps, L_hat),
    "LEM": lambda E0w, zmax, zsteps: propagate_LEM(E0w, zmax, zsteps, L_hat, goal_error=1e-7),
}

results = {name: [] for name in methods}

for name, method in methods.items():
    print(f"  Method: {name}")
    for i, (ns, dz) in enumerate(zip(z_steps_list, dz_list)):
        try:
            E_num_w = method(E0_w.copy(), z_max, ns)
            err = compute_rms_error(E_num_w, E_exact_w)
            results[name].append((dz, err))
            print(f"    dz={dz:.6f}, error={err:.2e}")
        except Exception as e:
            print(f"    dz={dz:.6f}, FAILED: {e}")
            results[name].append((dz, np.nan))

# Save data
os.makedirs("../data", exist_ok=True)
with open("../data/fig2_nse_soliton_errors.csv", "w", newline="") as f:
    writer = csv.writer(f)
    header = ["dz"]
    for name in methods:
        header.append(f"rms_error_{name}")
    writer.writerow(header)

    max_len = max(len(results[n]) for n in methods)
    for i in range(max_len):
        row = [results[list(methods.keys())[0]][i][0] if i < len(results[list(methods.keys())[0]]) else ""]
        for name in methods:
            if i < len(results[name]):
                row.append(results[name][i][1])
            else:
                row.append("")
        writer.writerow(row)

# Plot
fig, ax = plt.subplots(1, 1, figsize=(6, 5))

markers = {"SiSSM": "o", "SySSM": "s", "IFM-RK4IP": "^", "LEM": "d"}
colors = {"SiSSM": "C0", "SySSM": "C1", "IFM-RK4IP": "C2", "LEM": "C3"}

for name in methods:
    dzs = [r[0] for r in results[name] if not np.isnan(r[1])]
    errs = [r[1] for r in results[name] if not np.isnan(r[1])]
    if dzs:
        ax.loglog(dzs, errs, marker=markers[name], color=colors[name], label=name, markersize=5)

# Reference slopes
dz_ref = np.array([1e-3, 1e-1])
ax.loglog(dz_ref, 1e-3 * (dz_ref/dz_ref[0])**1, "k--", alpha=0.3, label=r"$O(\Delta z)$")
ax.loglog(dz_ref, 1e-5 * (dz_ref/dz_ref[0])**2, "k-.", alpha=0.3, label=r"$O(\Delta z^2)$")
ax.loglog(dz_ref, 1e-8 * (dz_ref/dz_ref[0])**4, "k:", alpha=0.3, label=r"$O(\Delta z^4)$")

ax.set_xlabel(r"Step size $\Delta z$")
ax.set_ylabel("Average RMS error")
ax.set_title("NSE Soliton: Error Convergence (Fig. 2)")
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig("../plots/fig2_nse_soliton_errors.pdf", dpi=150)
fig.savefig("../plots/fig2_nse_soliton_errors.png", dpi=150)
plt.close(fig)

print("Done! Saved fig2_nse_soliton_errors.csv and fig2_nse_soliton_errors.pdf")
