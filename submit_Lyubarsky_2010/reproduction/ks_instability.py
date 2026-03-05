"""
Core module for Kruskal-Schwarzschild instability computations.

Based on: Lyubarsky (2010), "A new mechanism for dissipation of alternating
fields in Poynting dominated outflows"

Key equations:
- Eq. 26: Full dispersion relation for K-S instability
- Eq. 28: Short wavelength growth rate (kDelta >> 1)
- Eq. 29: Long wavelength growth rate (kDelta << 1)
- Eq. 43: Lorentz factor evolution gamma(r)
- Eq. 44: Dissipation scale r_diss
"""

import numpy as np


# --- Section 2: Kruskal-Schwarzschild instability ---

def enthalpy_ratio():
    """
    Eq. 27: Ratio (h0 + h0') / (h0 - h0') = 3
    for strongly magnetized domain (h0 = 2p0) and
    relativistically hot slab (h0' = 4p0).
    """
    return 3.0


def dispersion_relation_omega4(k, g, Delta):
    """
    Eq. 26: Full dispersion relation.
    omega^4 = k^2 g^2 * [1 - exp(-2 k Delta)] /
              [((h0+h0')/(h0-h0'))^2 - exp(-2 k Delta)]

    Returns omega^4 (positive means unstable).
    """
    R = enthalpy_ratio()  # = 3
    x = 2.0 * k * Delta
    exp_neg = np.exp(-x)
    numerator = 1.0 - exp_neg
    denominator = R**2 - exp_neg  # 9 - exp(-2kDelta)
    return k**2 * g**2 * numerator / denominator


def growth_rate_exact(k, g, Delta):
    """
    Exact growth rate eta = -Im(omega) from Eq. 26.
    eta = (omega^4)^(1/4) since omega^4 > 0 for instability.
    """
    w4 = dispersion_relation_omega4(k, g, Delta)
    return np.abs(w4)**0.25


def growth_rate_short_wavelength(k, g):
    """
    Eq. 28: Short wavelength limit (k*Delta >> 1).
    eta = sqrt(k*g / 3)
    """
    return np.sqrt(k * g / 3.0)


def growth_rate_long_wavelength(k, g, Delta):
    """
    Eq. 29: Long wavelength limit (k*Delta << 1).
    eta = (g/2)^(1/2) * Delta^(1/4) * k^(3/4)
    """
    return np.sqrt(g / 2.0) * Delta**0.25 * k**0.75


# --- Section 3: Self-consistent acceleration/dissipation ---

def lorentz_factor(r, zeta, gamma_max, l):
    """
    Eq. 43: Lorentz factor evolution.
    gamma = (9 * zeta^2 * gamma_max * r / (4 * l))^(1/3)

    Parameters:
        r: distance from source (array or scalar)
        zeta: numerical reconnection efficiency factor
        gamma_max: sigma_0 * gamma_0, maximum achievable Lorentz factor
        l: stripe width
    """
    return (9.0 * zeta**2 * gamma_max * r / (4.0 * l))**(1.0 / 3.0)


def current_sheet_width_ratio(gamma, gamma_0, sigma_0):
    """
    Eq. 41: Delta/l = gamma / (3 * sigma_0 * gamma_0)
    Valid when Delta/l >> 1/gamma_0^2.
    """
    return gamma / (3.0 * sigma_0 * gamma_0)


def dissipation_scale(gamma_max, zeta, l):
    """
    Eq. 44: Dissipation scale.
    r_diss = 12 * (gamma_max / zeta)^2 * l
    """
    return 12.0 * (gamma_max / zeta)**2 * l


def solve_gamma_ode(r_array, zeta, sigma_0, gamma_0, l):
    """
    Solve the ODE (Eq. 42):
    gamma^2 * dgamma/dr = 3 * zeta^2 * sigma_0 * gamma_0 / (4 * l)

    with gamma(r_0) = gamma_0.
    Returns gamma(r) array.
    """
    rhs = 3.0 * zeta**2 * sigma_0 * gamma_0 / (4.0 * l)
    # gamma^3 = gamma_0^3 + 3 * rhs * (r - r_0)
    r0 = r_array[0]
    gamma_cubed = gamma_0**3 + 3.0 * rhs * (r_array - r0)
    return np.cbrt(gamma_cubed)
