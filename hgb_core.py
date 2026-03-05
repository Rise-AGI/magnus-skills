"""
Core computation module for hollow Gaussian beam propagation on curved surfaces.

Paper: "The hollow Gaussian beam propagation on curved surface based on matrix optics method"
Authors: Weifeng Ding, Zhaoying Wang (2021)

Key equations implemented:
- ABCD transfer matrix for CGCS (Eq. 5)
- Collins integral analytic solution (Eq. 12)
- Gaussian beam waist width (Eq. 14)
- Second-order moment beam width (Eq. 16)
- On-axis intensity for integer n (Eq. 18) and general n (Eq. 21)
- Dark spot size type A (Eq. 19) and type B (Eq. 20)
"""

import math
import numpy as np
from scipy.special import gamma as gamma_func, hyp1f1


def abcd_matrix(r0, r, theta):
    """
    ABCD transfer matrix for constant Gaussian curvature surface (Eq. 5).

    Parameters:
        r0: transverse parameter (radius at equator)
        r: curvature parameter
        theta: rotation angle (propagation coordinate)

    Returns:
        A, B, C, D matrix elements
    """
    phi = r0 * theta / r
    A = np.cos(phi)
    B = r * np.sin(phi)
    C = -np.sin(phi) / r
    D = np.cos(phi)
    return A, B, C, D


def abcd_free_space(z):
    """ABCD matrix for free-space propagation over distance z."""
    return 1.0, z, 0.0, 1.0


def electric_field_hgb(h2, n, sigma, k, A, B, D, r0_theta=0.0):
    """
    Analytic output field of HGB through ABCD system (Eq. 12).

    Parameters:
        h2: transverse output coordinate (array)
        n: order of HGB (can be non-integer)
        sigma: initial beam width parameter
        k: wave number
        A, B, D: ABCD matrix elements
        r0_theta: r0 * theta for phase factor

    Returns:
        Complex electric field E2(h2)
    """
    # Avoid division by zero when B=0
    if np.abs(B) < 1e-30:
        B = 1e-30

    alpha_param = 1.0 / sigma**2 - 1j * A * k / (2 * B)

    prefactor = (1.0 / np.sqrt(2 * np.pi)) * np.sqrt(k / (1j * B))
    prefactor *= np.exp(1j * D * k * h2**2 / (2 * B) + 1j * k * r0_theta)
    prefactor *= alpha_param**(-0.5 - n)
    prefactor *= (1.0 / sigma**2)**n
    prefactor *= gamma_func(0.5 + n)

    # Argument of confluent hypergeometric function
    z_arg = -h2**2 * k**2 * sigma**2 / (4 * B**2 - 2j * A * B * k * sigma**2)

    # Compute 1F1(1/2 + n, 1/2, z_arg)
    hf_values = np.zeros_like(h2, dtype=complex)
    for i, zi in enumerate(z_arg):
        hf_values[i] = hyp1f1(0.5 + n, 0.5, zi)

    return prefactor * hf_values


def intensity_hgb(h2, n, sigma, k, A, B, D, r0_theta=0.0):
    """Compute |E2|^2 for HGB propagation."""
    E = electric_field_hgb(h2, n, sigma, k, A, B, D, r0_theta)
    return np.abs(E)**2


def gaussian_beam_waist(sigma, k, A, B):
    """Gaussian beam (n=0) transverse waist width (Eq. 14)."""
    return sigma * np.sqrt(4 * B**2 / (k**2 * sigma**4) + A**2)


def divergence_coefficient(r, k, sigma):
    """Divergence coefficient gamma_d = 2r/(k*sigma^2) (below Eq. 15)."""
    return 2 * r / (k * sigma**2)


def second_order_moment_width(n, sigma, r0, r, theta, gamma_d):
    """Second-order moment beam width squared (Eq. 16)."""
    factor = (4 * n + 1) / 4.0
    sin2 = np.sin(r0 * theta / r)**2
    inner = (8 * n - 1) / (16 * n**2 - 1) * gamma_d**2 - 1
    return factor * sigma**2 * (1 + sin2 * inner)


def axial_intensity_integer_n(n, gamma_d, r0, r, theta):
    """
    On-axis intensity |E2|^2 at h2=0 for integer n (Eq. 18).
    Returns normalized intensity (without common prefactors).
    """
    sin_val = np.sin(r0 * theta / r)
    cot_val = np.cos(r0 * theta / r) / sin_val

    # Avoid division by zero
    mask = np.abs(sin_val) > 1e-15
    result = np.zeros_like(theta, dtype=float)

    numerator_coeff = ((math.factorial(2*n - 1) if n > 0 else 1) /
                       (2**(2*n - 1) * (math.factorial(n - 1) if n > 0 else 1)))**2 if n > 0 else 1.0

    if n == 0:
        # For n=0, Eq. 18 reduces to 1/|sin(r0*theta/r)|
        result[mask] = 1.0 / np.abs(sin_val[mask])
    else:
        result[mask] = (numerator_coeff *
                       gamma_d**(2*n) *
                       (gamma_d**2 + cot_val[mask]**2)**(-n - 0.5) /
                       np.abs(sin_val[mask]))

    return result


def axial_intensity_general_n(n, gamma_d, r0, r, theta):
    """
    On-axis intensity |E2|^2 at h2=0 for general n (Eq. 21).
    Includes cos^2(n*pi) factor for fractional orders.
    """
    sin_val = np.sin(r0 * theta / r)
    cot_val = np.cos(r0 * theta / r) / sin_val

    mask = np.abs(sin_val) > 1e-15
    result = np.zeros_like(theta, dtype=float)

    coeff = np.cos(n * np.pi)**2 / np.pi * gamma_func((1 + 2*n) / 2)**2

    result[mask] = (coeff * gamma_d**(2*n) *
                   (gamma_d**2 + cot_val[mask]**2)**(-n - 0.5) /
                   np.abs(sin_val[mask]))

    return result


def dark_spot_size_type_B(n, gamma_d, r0, r):
    """
    Longitudinal dark spot size for type B (Eq. 20).
    Valid when gamma_d > sqrt(2n+1).
    """
    arg = np.sqrt(2 * n / (gamma_d**2 - 1))
    if arg > 1:
        return np.nan
    return (r / r0) * np.arcsin(arg)


def compute_intensity_map(n, sigma, k, r0, r, theta_array, h_array):
    """
    Compute 2D intensity map I(theta, h) for HGB on CGCS.

    Parameters:
        n: HGB order
        sigma: initial beam width parameter
        k: wave number
        r0, r: CGCS parameters
        theta_array: 1D array of theta values
        h_array: 1D array of h values

    Returns:
        2D array of intensity values, shape (len(theta_array), len(h_array))
    """
    intensity_map = np.zeros((len(theta_array), len(h_array)))

    for i, theta in enumerate(theta_array):
        A, B, C, D = abcd_matrix(r0, r, theta)
        if np.abs(B) < 1e-30:
            continue
        intensity_map[i, :] = intensity_hgb(h_array, n, sigma, k, A, B, D, r0 * theta)

    return intensity_map


def compute_intensity_map_flat(n, sigma, k, z_array, h_array):
    """
    Compute 2D intensity map I(z, h) for HGB in flat space.

    Parameters:
        n: HGB order
        sigma: initial beam width parameter
        k: wave number
        z_array: 1D array of propagation distances
        h_array: 1D array of transverse positions

    Returns:
        2D array of intensity values
    """
    intensity_map = np.zeros((len(z_array), len(h_array)))

    for i, z in enumerate(z_array):
        A, B, C, D = abcd_free_space(z)
        if np.abs(B) < 1e-30:
            continue
        intensity_map[i, :] = intensity_hgb(h_array, n, sigma, k, A, B, D, 0.0)

    return intensity_map


def compute_skewness(n, sigma, k, r0, r, theta, h_array):
    """
    Compute skewness of intensity distribution at given theta.

    Skewness = E[(X - mu)^3] / std^3
    where X is the intensity-weighted position distribution.
    """
    A, B, C, D = abcd_matrix(r0, r, theta)
    if np.abs(B) < 1e-30:
        return 0.0

    I = intensity_hgb(h_array, n, sigma, k, A, B, D, r0 * theta)
    I_total = np.trapezoid(I, h_array)
    if I_total < 1e-30:
        return 0.0

    # Mean position
    mu = np.trapezoid(h_array * I, h_array) / I_total
    # Variance
    var = np.trapezoid((h_array - mu)**2 * I, h_array) / I_total
    if var < 1e-30:
        return 0.0
    std = np.sqrt(var)
    # Third central moment
    m3 = np.trapezoid((h_array - mu)**3 * I, h_array) / I_total

    return m3 / std**3
