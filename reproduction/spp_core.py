"""
Core physics module for SPP propagation around bends.

Implements the model from:
  Hasegawa, Nockel, Deutsch, Appl. Phys. Lett. 84, 1835 (2004)

Uses mpmath-based solver for the boundary matching equation (Eq. 2)
with continuation for efficient scanning over parameter ranges.
"""

import numpy as np
import mpmath

# Physical constants
c = 2.998e8        # speed of light [m/s]

# Silver Drude parameters (Johnson & Christy fit)
EPS_INF = 3.7      # high-frequency dielectric constant
OMEGA_P = 9.01     # plasma frequency [eV]
GAMMA_D = 0.048    # damping rate [eV]


def silver_epsilon(wavelength_nm):
    """Silver dielectric function using Drude model."""
    E_eV = 1240.0 / np.asarray(wavelength_nm, dtype=float)
    eps = EPS_INF - OMEGA_P**2 / (E_eV * (E_eV + 1j * GAMMA_D))
    return eps


def spp_wavevector(omega, eps_i, eps_o):
    """SPP wavevector k on a flat interface."""
    return (omega / c) * np.sqrt(eps_i * eps_o / (eps_i + eps_o) + 0j)


def decay_constants(omega, eps_i, eps_o):
    """Decay constants gamma_i (metal) and gamma_o (dielectric)."""
    gamma_i = -omega * eps_i * np.sqrt(-1.0 / (eps_i + eps_o) + 0j) / c
    gamma_o = omega * eps_o * np.sqrt(-1.0 / (eps_i + eps_o) + 0j) / c
    return gamma_i, gamma_o


def solve_mode_index(omega, R, eps_i, eps_o, n_guess=None, dps=15):
    """Solve boundary matching equation (Eq. 2) for mode index m.

    Uses mpmath Bessel functions with Newton iteration.
    Continuation-friendly: pass n_guess from previous solution.
    """
    mpmath.mp.dps = dps

    ki = omega * np.sqrt(eps_i + 0j) / c
    ko = omega * np.sqrt(eps_o + 0j) / c
    k = spp_wavevector(omega, eps_i, eps_o)

    ki_R = mpmath.mpc(float((ki * R).real), float((ki * R).imag))
    ko_R = mpmath.mpc(float((ko * R).real), float((ko * R).imag))
    ki_mp = mpmath.mpc(float(ki.real), float(ki.imag))
    ko_mp = mpmath.mpc(float(ko.real), float(ko.imag))

    if n_guess is None:
        n_guess = k * R
    m = mpmath.mpc(float(n_guess.real), float(n_guess.imag))

    def boundary_residual(n):
        jn = mpmath.besselj(n, ki_R)
        jnm1 = mpmath.besselj(n - 1, ki_R)
        lhs = (jnm1 / jn - n / ki_R) / ki_mp

        jn_ko = mpmath.besselj(n, ko_R)
        yn_ko = mpmath.bessely(n, ko_R)
        hn = jn_ko + 1j * yn_ko
        jnm1_ko = mpmath.besselj(n - 1, ko_R)
        ynm1_ko = mpmath.bessely(n - 1, ko_R)
        hnm1 = jnm1_ko + 1j * ynm1_ko

        rhs = (hnm1 / hn - n / ko_R) / ko_mp
        return lhs - rhs

    tol = mpmath.mpf(10)**(-dps + 3)
    for iteration in range(25):
        f = boundary_residual(m)
        if abs(f) < tol:
            break
        h = m * mpmath.mpf('1e-8')
        if abs(h) < 1e-12:
            h = mpmath.mpc('1e-8', '1e-8')
        fp = (boundary_residual(m + h) - f) / h
        if abs(fp) < 1e-300:
            break
        step = f / fp
        m = m - step

    return complex(m)


def transmittance(m, k, R, theta):
    """Energy transmittance T (Eq. 3)."""
    kR = k * R
    num = 4.0 * m * kR
    denom = (-np.exp(1j * m * theta) * (m - kR)**2
             + np.exp(-1j * m * theta) * (m + kR)**2)
    return abs(num / denom)**2


def reflectance(m, k, R, theta):
    """Energy reflectance R from 1D potential well analogy."""
    kR = k * R
    denom = (-np.exp(1j * m * theta) * (m - kR)**2
             + np.exp(-1j * m * theta) * (m + kR)**2)
    r = -(m**2 - kR**2) * 2j * np.sin(m * theta) / denom
    return abs(r)**2


def upper_bound_transmittance(m, theta):
    """Upper bound transmittance Tu (Eq. 6)."""
    return np.exp(-2.0 * m.imag * theta)


def flat_surface_loss(omega, eps_i, eps_o, arc_length):
    """Energy remaining after SPP propagation on flat surface."""
    k = spp_wavevector(omega, eps_i, eps_o)
    return np.exp(-2.0 * k.imag * arc_length)


def compute_Tu_vs_R(wavelength_nm, R_array, theta, eps_o=1.0):
    """Compute upper bound transmittance vs bend radius with continuation."""
    eps_i = silver_epsilon(wavelength_nm)
    omega = 2.0 * np.pi * c / (wavelength_nm * 1e-9)
    k = spp_wavevector(omega, eps_i, eps_o)

    Tu_array = np.zeros(len(R_array))
    m_array = np.zeros(len(R_array), dtype=complex)

    m_prev = None
    for i, R in enumerate(R_array):
        if m_prev is not None:
            # Use continuation: scale previous m by ratio of R values
            n_guess = m_prev * (R / R_array[max(0, i-1)])
        else:
            n_guess = k * R

        try:
            m = solve_mode_index(omega, R, eps_i, eps_o, n_guess=n_guess)
            Tu_array[i] = upper_bound_transmittance(m, theta)
            m_array[i] = m
            m_prev = m
        except Exception as e:
            Tu_array[i] = np.nan
            m_array[i] = np.nan

    return Tu_array, m_array
