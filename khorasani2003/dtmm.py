"""
Core module for Differential Transfer Matrix Method (DTMM) computation
of photonic crystal band structures.

Implements the formalism from:
  Khorasani & Adibi, "Analytical determination of bandgaps in photonic crystals"

Methods:
  - Transfer matrix computation for continuous index profiles via thin-layer slicing
  - Transfer matrix for stratified (step) profiles with proper interface conditions
  - Band structure computation for TE and TM polarizations
"""

import numpy as np
from scipy.linalg import expm


# --- Refractive index profiles ---

def n_sinusoidal(x, L=1.0):
    """Sinusoidal profile: n(x) = 2 + cos(2*pi*x/L), varies 1 to 3, even symmetric."""
    return 2.0 + np.cos(2.0 * np.pi * x / L)


def n_triangular(x, L=1.0):
    """Triangular profile: peaks at center (n=3), min at edges (n=1), even symmetric."""
    return 3.0 - 4.0 * np.abs(x) / L


def n_square(x, L=1.0):
    """Square (step) profile: n=3 for |x|<L/4, n=1 otherwise, even symmetric."""
    return np.where(np.abs(x) < L / 4.0, 3.0, 1.0)


def n_ramp(x, L=1.0):
    """Asymmetric ramp with jump: n(x) = 2 + 2x/L for -L/2 <= x < L/2.
    Ramps from 1 to 3, then jumps back to 1. NOT even symmetric."""
    return 2.0 + 2.0 * x / L


# --- Transfer matrix computation ---

def compute_transfer_matrix_te(omega_norm, n_func, theta, L=1.0, N=500):
    """
    Compute the one-period transfer matrix for TE polarization using thin-layer slicing.

    Parameters:
        omega_norm: normalized frequency omega*L/c
        n_func: refractive index function n(x)
        theta: angle of incidence (radians)
        L: period length
        N: number of thin layers

    Returns:
        T: 2x2 transfer matrix for one period
    """
    k0 = omega_norm / L  # k0 = omega/c, with omega*L/c = omega_norm
    beta = k0 * np.sin(theta)  # propagation constant along z
    dx = L / N
    x_vals = np.linspace(-L / 2 + dx / 2, L / 2 - dx / 2, N)

    T = np.eye(2, dtype=complex)
    for x in x_vals:
        n = n_func(x, L)
        kx2 = (k0 * n) ** 2 - beta ** 2
        if kx2 > 0:
            kx = np.sqrt(kx2)
            c = np.cos(kx * dx)
            s = np.sin(kx * dx)
            Ti = np.array([[c, s / kx], [-kx * s, c]], dtype=complex)
        else:
            kappa = np.sqrt(-kx2)
            ch = np.cosh(kappa * dx)
            sh = np.sinh(kappa * dx)
            Ti = np.array([[ch, sh / kappa], [kappa * sh, ch]], dtype=complex)
        T = Ti @ T
    return T


def compute_transfer_matrix_tm(omega_norm, n_func, theta, L=1.0, N=500):
    """
    Compute the one-period transfer matrix for TM polarization.

    For TM, the state vector is [H, (1/n^2)*dH/dx], and interface conditions differ.
    Within each thin layer, the transfer matrix uses kx but the off-diagonal terms
    are scaled by 1/n^2.
    """
    k0 = omega_norm / L
    beta = k0 * np.sin(theta)
    dx = L / N
    x_vals = np.linspace(-L / 2 + dx / 2, L / 2 - dx / 2, N)

    T = np.eye(2, dtype=complex)
    for x in x_vals:
        n = n_func(x, L)
        kx2 = (k0 * n) ** 2 - beta ** 2
        if kx2 > 0:
            kx = np.sqrt(kx2)
            c = np.cos(kx * dx)
            s = np.sin(kx * dx)
            # TM: state vector [H, (1/n^2)*dH/dx]
            # The transfer matrix in this basis:
            Ti = np.array([
                [c, n ** 2 * s / kx],
                [-kx * s / n ** 2, c]
            ], dtype=complex)
        else:
            kappa = np.sqrt(-kx2)
            ch = np.cosh(kappa * dx)
            sh = np.sinh(kappa * dx)
            Ti = np.array([
                [ch, n ** 2 * sh / kappa],
                [kappa * sh / n ** 2, ch]
            ], dtype=complex)
        T = Ti @ T
    return T


def compute_transfer_matrix_step_te(omega_norm, n_vals, d_vals, theta):
    """
    Transfer matrix for a stratified medium (piecewise constant n) in one period, TE.

    Parameters:
        omega_norm: normalized frequency omega*L/c
        n_vals: list of refractive indices [n1, n2, ...]
        d_vals: list of layer thicknesses [d1, d2, ...] (must sum to L)
        theta: angle of incidence

    Returns:
        T: 2x2 transfer matrix
    """
    L = sum(d_vals)
    k0 = omega_norm / L
    beta = k0 * np.sin(theta)

    T = np.eye(2, dtype=complex)
    for n, d in zip(n_vals, d_vals):
        kx2 = (k0 * n) ** 2 - beta ** 2
        if kx2 > 0:
            kx = np.sqrt(kx2)
            c = np.cos(kx * d)
            s = np.sin(kx * d)
            Ti = np.array([[c, s / kx], [-kx * s, c]], dtype=complex)
        else:
            kappa = np.sqrt(-kx2)
            ch = np.cosh(kappa * d)
            sh = np.sinh(kappa * d)
            Ti = np.array([[ch, sh / kappa], [kappa * sh, ch]], dtype=complex)
        T = Ti @ T
    return T


def compute_transfer_matrix_step_tm(omega_norm, n_vals, d_vals, theta):
    """Transfer matrix for stratified medium, TM polarization."""
    L = sum(d_vals)
    k0 = omega_norm / L
    beta = k0 * np.sin(theta)

    T = np.eye(2, dtype=complex)
    for n, d in zip(n_vals, d_vals):
        kx2 = (k0 * n) ** 2 - beta ** 2
        if kx2 > 0:
            kx = np.sqrt(kx2)
            c = np.cos(kx * d)
            s = np.sin(kx * d)
            Ti = np.array([
                [c, n ** 2 * s / kx],
                [-kx * s / n ** 2, c]
            ], dtype=complex)
        else:
            kappa = np.sqrt(-kx2)
            ch = np.cosh(kappa * d)
            sh = np.sinh(kappa * d)
            Ti = np.array([
                [ch, n ** 2 * sh / kappa],
                [kappa * sh / n ** 2, ch]
            ], dtype=complex)
        T = Ti @ T
    return T


# --- Band structure computation ---

def compute_band_structure(omega_range, n_func, theta, polarization='TE',
                           L=1.0, N=500, is_step=False, n_vals=None, d_vals=None):
    """
    Compute cos(kappa*L) for a range of normalized frequencies.

    Parameters:
        omega_range: array of normalized frequencies omega*L/c
        n_func: refractive index function (ignored if is_step=True)
        theta: angle of incidence
        polarization: 'TE' or 'TM'
        L: period
        N: thin-layer count for continuous profiles
        is_step: if True, use stratified medium formulas
        n_vals, d_vals: layer parameters for step profile

    Returns:
        kappa_norm: array of kappa*L/pi values (NaN where in bandgap)
        cos_kL: array of cos(kappa*L) values
    """
    cos_kL = np.zeros(len(omega_range))

    for i, omega_n in enumerate(omega_range):
        if omega_n <= 0:
            cos_kL[i] = 1.0
            continue

        if is_step:
            if polarization == 'TE':
                T = compute_transfer_matrix_step_te(omega_n, n_vals, d_vals, theta)
            else:
                T = compute_transfer_matrix_step_tm(omega_n, n_vals, d_vals, theta)
        else:
            if polarization == 'TE':
                T = compute_transfer_matrix_te(omega_n, n_func, theta, L, N)
            else:
                T = compute_transfer_matrix_tm(omega_n, n_func, theta, L, N)

        cos_kL[i] = np.real(0.5 * np.trace(T))

    # Compute kappa*L/pi
    kappa_norm = np.full_like(omega_range, np.nan)
    mask = np.abs(cos_kL) <= 1.0
    kappa_norm[mask] = np.arccos(cos_kL[mask]) / np.pi

    return kappa_norm, cos_kL


def compute_bandstructure_full(omega_max, n_func, theta, polarization='TE',
                                L=1.0, N=500, n_points=1000,
                                is_step=False, n_vals=None, d_vals=None):
    """
    Compute the full band structure up to omega_max.
    Returns arrays suitable for plotting frequency vs. wavenumber.
    """
    omega_range = np.linspace(0.001, omega_max, n_points)
    kappa_norm, cos_kL = compute_band_structure(
        omega_range, n_func, theta, polarization, L, N, is_step, n_vals, d_vals
    )
    return omega_range, kappa_norm, cos_kL
