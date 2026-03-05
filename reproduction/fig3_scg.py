"""
Figure 3: Supercontinuum generation in a photonic crystal fiber (PCF).

Reproduces the SCG simulation using the simplified forward model for the
analytic signal with Raman effect (FMAS-S-R), with polynomial dispersion.

Reference: Melchert & Demircan, CPC (2022), Section 9.2, Figure 3.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.fft import fft, ifft, fftfreq, fftshift, ifftshift
import math
import csv
import os

# --- Physical parameters ---
# Propagation constant: polynomial approximation Eq. (29)
# beta(w) = sum_{n=2}^{10} beta_n / n! * (w - w0)^n
# with beta_n from Table 2

beta_n = {
    2: -0.011830,
    3: 0.081038,
    4: -0.095205,
    5: 0.20737,
    6: -0.53943,
    7: 1.34860,
    8: -2.5495,
    9: 3.0524,
    10: -1.7140,
}  # units: fs^n/micron

# Reference frequency
w0 = np.pi  # rad/fs (approximately 375 THz, ~800 nm)

# Nonlinear parameters
gamma = 0.11  # W^-1 m^-1 => 0.11e-6 W^-1 um^-1
# n2*w/c = gamma, but we work in consistent units
# Paper: n2*c = w0 with chi = 0.11e-6 W^-1 m^-1

# Raman parameters
f_R = 0.18
tau1 = 12.2  # fs
tau2 = 32.0  # fs

# --- Domain parameters ---
t_max = 2000.0   # fs, half-period
t_num = 2**14     # 16384 grid points
z_max = 0.14e6    # 14 cm in microns (0.14 m = 140000 um)
z_steps = 3500    # number of z-steps (dz ~ 40 um)
z_skip = 10       # store every 10th step

# --- Grid setup ---
dt = 2 * t_max / t_num
t = np.linspace(-t_max, t_max - dt, t_num)
dw = np.pi / t_max
w = 2 * np.pi * fftfreq(t_num, d=dt)

# --- Dispersion relation ---
def beta_poly(w_arr, w0, beta_n_dict):
    """Compute beta(w) using polynomial expansion."""
    dw = w_arr - w0
    result = np.zeros_like(w_arr, dtype=float)
    for n, bn in beta_n_dict.items():
        result += bn / math.factorial(n) * dw**n
    return result

beta_w = beta_poly(w, w0, beta_n)

# Linear operator: L_hat = i * k(w) where k(w) = -beta(w) for the model
# For FMAS-S-R: i*dz E_w + beta(w)*E_w + NL = 0
# => dz E_w = i*beta(w)*E_w + NL_contribution
# Actually: dz E_w = -i*beta(w)*E_w - i*NL
# The linear operator is: L_hat(w) = -i * beta(w)
# But since beta already has the right sign convention from the expansion:
L_hat = 1j * beta_w  # This is -i*k(w) = i*beta(w) with the right sign

# For the FMAS_S model, the nonlinear part is:
# N(E) = -i * n2 * w / c * |E|^2 * E  (for positive freq only)
# We simplify: gamma = n2 * w0 / c
# N_hat = -i * gamma * FFT(|E|^2 * E) * (w > 0 factor)

# --- Raman response in frequency domain ---
# h(w) = (tau1^-2 + tau2^-2) * tau1^2 / (tau1^(-2) - (w + i*tau2^(-1))^2)
def raman_response_freq(w_arr, tau1, tau2):
    """Raman response in frequency domain, Eq. (10)."""
    h_w = (tau1**(-2) + tau2**(-2)) * tau1**2 / (tau1**(-2) - (w_arr + 1j / tau2)**2)
    return h_w

h_w = raman_response_freq(w, tau1, tau2)

# --- Initial condition ---
# Soliton-like pulse for SCG:
# E(0,t) = sqrt(P0) * sech(t/t0) * exp(-i*w0*t)
# From the paper context, typical SCG uses:
P0 = 10000.0  # W (10 kW peak power)
t0 = 28.4     # fs (pulse duration parameter for N~8 soliton)

E0_real = np.sqrt(P0) / np.cosh(t / t0)
# The field is the real-valued optical field
E0 = E0_real.copy()

# Compute analytic signal
E0_w = fft(E0)
# Set negative frequencies to zero (analytic signal)
E0_w[w < 0] = 0
E0_w[w == 0] *= 0.5  # handle DC
E0_w *= 2  # normalize

# --- Nonlinear operator ---
def compute_nonlinear(E_w, w_arr, gamma_val, f_R, h_w, w0_ref):
    """Compute nonlinear term for FMAS-S-R model."""
    E_t = ifft(E_w)
    intensity = np.abs(E_t)**2

    # Instantaneous Kerr + Raman
    # NL = gamma * w/w0 * [(1-fR)*|E|^2*E + fR*E*IFT(h_w * FT(|E|^2))]
    raman_conv = ifft(h_w * fft(intensity))
    NL_t = (1 - f_R) * intensity * E_t + f_R * E_t * raman_conv

    NL_w = fft(NL_t)

    # Frequency weighting: w/w0 for each frequency component
    w_factor = np.abs(w_arr) / w0_ref
    w_factor[w_arr == 0] = 0

    # Apply to positive frequencies only
    mask = w_arr > 0
    result = np.zeros_like(E_w)
    result[mask] = 1j * gamma_val * w_factor[mask] * NL_w[mask]

    return result


# --- Propagation using IFM-RK4IP ---
print("Computing supercontinuum generation (Figure 3)...")
print(f"  Domain: t in [{-t_max}, {t_max}] fs, {t_num} points")
print(f"  Propagation: z_max = {z_max/1e6:.4f} m, {z_steps} steps")

# Convert gamma to consistent units (um)
# gamma = 0.11 W^-1 m^-1 = 0.11e-6 W^-1 um^-1
gamma_um = 0.11e-6  # W^-1 um^-1

dz = z_max / z_steps
half_lin = np.exp(L_hat * dz / 2.0)
inv_half_lin = np.exp(-L_hat * dz / 2.0)

E_w = E0_w.copy()

# Storage
n_stored = z_steps // z_skip + 1
z_stored = np.zeros(n_stored)
intensity_stored = np.zeros((n_stored, t_num))
spectrum_stored = np.zeros((n_stored, t_num))

# Store initial
E_t = ifft(E_w)
z_stored[0] = 0
intensity_stored[0] = np.abs(E_t)**2
spectrum_stored[0] = np.abs(E_w)**2

store_idx = 1

for step in range(1, z_steps + 1):
    # IFM-RK4IP step (correct IP with z0 = z + dz/2)
    half = half_lin
    inv_half = inv_half_lin
    psi_w = half * E_w

    # k1 at z_n: G = half * N(E_w)
    k1 = dz * half * compute_nonlinear(E_w, w, gamma_um, f_R, h_w, w0)
    # k2 at midpoint: G = N(psi + k1/2)
    k2 = dz * compute_nonlinear(psi_w + 0.5 * k1, w, gamma_um, f_R, h_w, w0)
    # k3 at midpoint
    k3 = dz * compute_nonlinear(psi_w + 0.5 * k2, w, gamma_um, f_R, h_w, w0)
    # k4 at z_n+dz: G = inv_half * N(half * (psi + k3))
    k4 = dz * inv_half * compute_nonlinear(half * (psi_w + k3), w, gamma_um, f_R, h_w, w0)

    psi_w = psi_w + (k1 + 2*k2 + 2*k3 + k4) / 6.0
    E_w = half * psi_w

    if step % z_skip == 0 and store_idx < n_stored:
        E_t = ifft(E_w)
        z_stored[store_idx] = step * dz
        intensity_stored[store_idx] = np.abs(E_t)**2
        spectrum_stored[store_idx] = np.abs(E_w)**2
        store_idx += 1

    if step % 500 == 0:
        print(f"  Step {step}/{z_steps}, z = {step*dz/1e6:.4f} m")

print(f"  Propagation complete. Stored {store_idx} snapshots.")

# --- Save data ---
os.makedirs("../data", exist_ok=True)

# Save final field and spectrum
E_final_t = ifft(E_w)
with open("../data/fig3_scg_final_field.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["t_fs", "Re_E", "Im_E", "intensity"])
    for i in range(0, t_num, 4):  # downsample for file size
        writer.writerow([t[i], np.real(E_final_t[i]), np.imag(E_final_t[i]), np.abs(E_final_t[i])**2])

with open("../data/fig3_scg_final_spectrum.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["w_rad_fs", "spectrum_dB"])
    spec = np.abs(E_w)**2
    spec_dB = 10 * np.log10(spec / np.max(spec) + 1e-30)
    w_sorted = fftshift(w)
    spec_sorted = fftshift(spec_dB)
    mask = w_sorted > 0
    for i in range(0, np.sum(mask), 4):
        idx = np.where(mask)[0][i]
        writer.writerow([w_sorted[idx], spec_sorted[idx]])

# --- Plot ---
fig, axes = plt.subplots(2, 3, figsize=(14, 8))

# (a) Dispersion curves
ax = axes[0, 0]
w_plot = np.linspace(0.5, 5.0, 1000)
beta1 = np.zeros_like(w_plot)
beta2 = np.zeros_like(w_plot)
for i, wp in enumerate(w_plot):
    dw_val = wp - w0
    for n, bn in beta_n.items():
        if n >= 2:
            beta1[i] += bn / math.factorial(n-1) * dw_val**(n-1)
        if n >= 2:
            beta2[i] += bn / math.factorial(n-2) * dw_val**(n-2) if n >= 2 else 0
    # Recompute beta2 properly
    beta2_val = 0
    for n, bn in beta_n.items():
        if n >= 2:
            beta2_val += bn / math.factorial(n-2) * dw_val**(n-2)
    beta2[i] = beta2_val

ax.plot(w_plot, beta1, 'b-', label=r'$\beta_1(\omega)$')
ax2 = ax.twinx()
ax2.plot(w_plot, beta2, 'r-', label=r'$\beta_2(\omega)$')
ax2.axhline(0, color='gray', ls='--', lw=0.5)
ax.set_xlabel(r'$\omega$ (rad/fs)')
ax.set_ylabel(r'$\beta_1$ (fs/$\mu$m)', color='b')
ax2.set_ylabel(r'$\beta_2$ (fs$^2$/$\mu$m)', color='r')
ax.set_title('(a) Dispersion')

# (b) Final optical field
ax = axes[0, 1]
E_final_real = np.real(E_final_t)
ax.plot(t, E_final_real, 'b-', lw=0.3)
ax.set_xlim(-500, 2000)
ax.set_xlabel('t (fs)')
ax.set_ylabel(r'$E$ (a.u.)')
ax.set_title(f'(b) Field at z={z_max/1e6:.2f} m')

# (c) Final spectrum
ax = axes[0, 2]
spec = np.abs(E_w)**2
spec_dB = 10 * np.log10(spec / np.max(spec) + 1e-30)
w_s = fftshift(w)
spec_s = fftshift(spec_dB)
mask = (w_s > 0.5) & (w_s < 6)
ax.plot(w_s[mask], spec_s[mask], 'b-', lw=0.5)
ax.set_xlabel(r'$\omega$ (rad/fs)')
ax.set_ylabel('Spectrum (dB)')
ax.set_ylim(-80, 5)
ax.set_title(f'(c) Spectrum at z={z_max/1e6:.2f} m')

# (d) Spectrogram placeholder - compute STFT
ax = axes[1, 0]
# Simple spectrogram
sigma = 20.0  # fs, window width
n_t_spec = 200
n_w_spec = 200
t_spec = np.linspace(-500, 2000, n_t_spec)
w_spec = np.linspace(0.5, 5.0, n_w_spec)
P_spec = np.zeros((n_w_spec, n_t_spec))

for it, t0_val in enumerate(t_spec):
    window = np.exp(-(t - t0_val)**2 / (2 * sigma**2))
    windowed = E_final_t * window
    windowed_w = fft(windowed)
    for iw, w0_val in enumerate(w_spec):
        idx = np.argmin(np.abs(w - w0_val))
        P_spec[iw, it] = np.abs(windowed_w[idx])**2

P_spec_dB = 10 * np.log10(P_spec / np.max(P_spec) + 1e-30)
ax.pcolormesh(t_spec, w_spec, P_spec_dB, vmin=-40, vmax=0, cmap='jet', shading='auto')
ax.set_xlabel('t (fs)')
ax.set_ylabel(r'$\omega$ (rad/fs)')
ax.set_title(f'(d) Spectrogram at z={z_max/1e6:.2f} m')

# (e) Intensity evolution
ax = axes[1, 1]
z_cm = z_stored[:store_idx] / 1e4  # um to cm
t_lim = (-200, 800)
t_mask = (t > t_lim[0]) & (t < t_lim[1])
ax.pcolormesh(t[t_mask], z_cm, intensity_stored[:store_idx][:, t_mask],
              cmap='hot', shading='auto')
ax.set_xlabel('t (fs)')
ax.set_ylabel('z (cm)')
ax.set_title('(e) Intensity evolution')

# (f) Spectrum evolution
ax = axes[1, 2]
w_s = fftshift(w)
spec_evol = fftshift(spectrum_stored[:store_idx], axes=1)
spec_evol_dB = 10 * np.log10(spec_evol / np.max(spec_evol) + 1e-30)
w_mask = (w_s > 0.5) & (w_s < 5.5)
ax.pcolormesh(w_s[w_mask], z_cm, spec_evol_dB[:, w_mask],
              vmin=-60, vmax=0, cmap='jet', shading='auto')
ax.set_xlabel(r'$\omega$ (rad/fs)')
ax.set_ylabel('z (cm)')
ax.set_title('(f) Spectrum evolution')

fig.suptitle("Supercontinuum Generation in PCF (Fig. 3)", fontsize=14)
fig.tight_layout()
fig.savefig("../plots/fig3_scg.pdf", dpi=150)
fig.savefig("../plots/fig3_scg.png", dpi=150)
plt.close(fig)

print("Done! Saved fig3 data and plots.")
