"""
Core plasma physics functions for reproducing
D. Zhou, Phys. Plasmas 23, 070701 (2016):
"The General Solution to Vlasov Equation and Linear Landau Damping"

All figure scripts import from this module.
"""

import numpy as np
from scipy.special import wofz
from scipy.constants import e as e_charge, m_e, epsilon_0, pi
from scipy.optimize import brentq
import warnings

warnings.filterwarnings("ignore", category=RuntimeWarning)


# ---------------------------------------------------------------------------
# Physical constants and default parameters
# ---------------------------------------------------------------------------

class PlasmaParams:
    """Parameters for a 1D electrostatic unmagnetized electron plasma."""

    def __init__(self, n0=1e18, T_eV=10.0):
        self.n0 = n0            # background density [m^-3]
        self.T_eV = T_eV        # electron temperature [eV]
        self.T_e = T_eV * e_charge  # temperature [J]

        # Thermal velocity v_th = sqrt(T_e / m_e)
        self.v_th = np.sqrt(self.T_e / m_e)

        # Plasma frequency
        self.omega_pe = np.sqrt(n0 * e_charge**2 / (epsilon_0 * m_e))

        # Debye length
        self.lambda_D = self.v_th / self.omega_pe


# ---------------------------------------------------------------------------
# Plasma Dispersion Function Z(zeta) and its derivative
# ---------------------------------------------------------------------------

def plasma_Z(zeta):
    """Plasma dispersion function Z(zeta) = i*sqrt(pi)*w(zeta)
    where w is the Faddeeva function.

    Valid for all complex zeta, uses the Landau contour.
    """
    return 1j * np.sqrt(pi) * wofz(zeta)


def plasma_Z_prime(zeta):
    """Derivative Z'(zeta) = -2[1 + zeta*Z(zeta)]."""
    return -2.0 * (1.0 + zeta * plasma_Z(zeta))


# ---------------------------------------------------------------------------
# Dielectric function
# ---------------------------------------------------------------------------

def dielectric_landau(omega, k, params):
    """Full Landau dielectric function (Eq. 23 of the paper with Landau contour).

    D(omega,k) = 1 + 1/(k^2 lambda_D^2) * [1 + zeta * Z(zeta)]

    where zeta = omega / (sqrt(2) * k * v_th).
    """
    zeta = omega / (np.sqrt(2) * k * params.v_th)
    return 1.0 + (1.0 / (k * params.lambda_D)**2) * (1.0 + zeta * plasma_Z(zeta))


def dielectric_vlasov(omega_r, k, params):
    """Vlasov dielectric function (Eq. 1 of the paper) for real omega.

    Uses principal-value integral only (no Landau correction S term).
    D_V = 1 + (omega_pe^2)/(n0*k^2) * P.V. int (df0/dv)/(v - omega/k) dv

    For a Maxwellian, this gives:
    D_V = 1 + 1/(k*lambda_D)^2 * [1 + zeta * Re(Z(zeta))]
    where zeta is real (omega is real).

    Note: Im part from principal value = 0 for real zeta.
    But Re(Z(zeta_r)) involves the principal value integral.
    """
    zeta = omega_r / (np.sqrt(2) * k * params.v_th)
    Z_val = plasma_Z(zeta + 0j)
    # Principal value contribution is the real part of Z for real zeta
    # The full Landau result adds i*pi*sign correction
    return 1.0 + (1.0 / (k * params.lambda_D)**2) * (1.0 + zeta * np.real(Z_val))


def dielectric_landau_real_omega(omega_r, k, params):
    """Landau dielectric function evaluated at real omega.

    For gamma=0, S = pi*i, and the dispersion relation (Eq. 23) becomes:
    D = 1 + 1/(k*lD)^2 * [P.V. integral + pi*i * (df0/dv)|_{v=w/k}]
    which equals the full Z-function result.
    """
    zeta = omega_r / (np.sqrt(2) * k * params.v_th)
    Z_val = plasma_Z(zeta + 0j)
    chi = (1.0 / (k * params.lambda_D)**2) * (1.0 + zeta * Z_val)
    return 1.0 + chi


# ---------------------------------------------------------------------------
# Langmuir wave dispersion relation
# ---------------------------------------------------------------------------

def langmuir_dispersion_approx(k, params):
    """Approximate Langmuir wave dispersion: omega^2 = omega_pe^2 (1 + 3 k^2 lambda_D^2).

    Returns (omega_r, gamma) where gamma is the Landau damping rate.
    """
    kld = k * params.lambda_D
    omega_r = params.omega_pe * np.sqrt(1.0 + 3.0 * kld**2)

    # Landau damping from Eq. (27): gamma = (pi * Omega^3) / (2 * n0 * k^2) * df0/dv|_{v=w/k}
    # For Maxwellian: df0/dv = -n0*v/(sqrt(2*pi)*v_th^3) * exp(-v^2/(2*v_th^2))
    # At v = omega_r/k:
    v_res = omega_r / k
    df0_dv = -params.n0 * v_res / (np.sqrt(2 * pi) * params.v_th**3) * np.exp(
        -v_res**2 / (2 * params.v_th**2))

    gamma = pi * omega_r**3 / (2 * params.n0 * k**2) * df0_dv
    return omega_r, gamma


def solve_langmuir_dispersion(k, params):
    """Solve the full Landau dispersion relation numerically.

    Finds complex omega = Omega + i*gamma satisfying D(omega,k) = 0.
    Uses Newton's method starting from approximate solution.
    """
    omega_r0, gamma0 = langmuir_dispersion_approx(k, params)
    omega = complex(omega_r0, gamma0)

    for _ in range(200):
        D_val = dielectric_landau(omega, k, params)
        dw = abs(omega) * 1e-8 + 1.0
        D_plus = dielectric_landau(omega + dw, k, params)
        D_minus = dielectric_landau(omega - dw, k, params)
        dD = (D_plus - D_minus) / (2 * dw)

        if abs(dD) < 1e-50:
            break

        delta = -D_val / dD

        # Limit step
        max_step = 0.1 * abs(omega)
        if abs(delta) > max_step:
            delta = max_step * delta / abs(delta)

        omega = omega + delta

        if abs(delta / omega) < 1e-12:
            break

    return omega.real, omega.imag


# ---------------------------------------------------------------------------
# Quasi-linear diffusion coefficient (Eq. 19)
# ---------------------------------------------------------------------------

def diffusion_vlasov(v, omega_r, gamma, k, Ek2=1.0):
    """Vlasov quasi-linear diffusion coefficient (Eq. 3).

    D propto gamma / ((Omega - k*v)^2 + gamma^2)

    Returns normalized D (without e^2/m^2 |Ek|^2 prefactor).
    """
    return gamma / ((omega_r - k * v)**2 + gamma**2)


def diffusion_general(v, omega_r, gamma, k, params):
    """General quasi-linear diffusion coefficient (Eq. 19).

    D = Im[1/(kv - omega) + S * Delta(v, omega/k) / k]

    For gamma > 0: S = 0, reduce to Vlasov result
    For gamma = 0: S = pi*i
    For gamma < 0: S = 2*pi*i

    The Delta function series terms contribute corrections that ensure
    continuity at gamma = 0.
    """
    omega = complex(omega_r, gamma)
    kv_minus_omega = k * v - omega

    # First term: Im(1/(kv - omega)) = gamma / ((kv - Omega)^2 + gamma^2)
    term1 = np.imag(1.0 / kv_minus_omega)

    if gamma > 1e-20:
        S = 0.0
    elif gamma < -1e-20:
        S = 2.0 * pi
    else:
        S = pi

    if abs(S) < 1e-30:
        return term1

    # Delta function contributions (Eq. 10)
    # For smooth test: approximate delta functions with Gaussians
    # In the limit, delta(v - Omega/k) -> narrow Gaussian
    # But for the diffusion coefficient in the continuous sense:
    #
    # The full expression Eq. (19) is:
    # D = (e/m)^2 |Ek|^2 { gamma/((Omega-kv)^2+gamma^2) + S/(ik) * [delta(v-Omega/k)
    #   - (gamma/k)^2 delta''(v-Omega/k) + ...]}
    #
    # For a smooth evaluation at finite gamma, we can use the series result which gives:
    # D_total = pi * delta(kv - Omega) for gamma -> 0+ and gamma -> 0-
    #
    # But for visualization we evaluate numerically.
    # The key insight is: term1 + S_correction -> pi * delta(kv-Omega) as gamma->0
    #
    # Using the explicit formula, the total is:
    # Im[1/(kv - omega)] + S * Re[Delta_smooth]
    # where Delta_smooth is the distributional kernel evaluated at finite gamma.
    #
    # For numerical evaluation with the series Delta:
    # Sum_n ((-i*gamma/k)^n / n!) * delta^(n)(v - Omega/k)
    # At finite gamma, approximate delta -> narrow Gaussian with width sigma = |gamma|/(2k)
    sigma = max(abs(gamma) / (2.0 * k), 1e-30 * abs(omega_r) / k)
    v_res = omega_r / k
    x = (v - v_res) / sigma

    # Delta series with Gaussian regularization
    gauss = np.exp(-0.5 * x**2) / (sigma * np.sqrt(2 * pi))

    # Leading term of S * Delta / k
    # S/(ik) * delta(v - Omega/k) -> S/(ik) * gauss
    # Im part: S/k * gauss (since S is real, i factor gives real part of 1/(ik))
    # Actually S/(ik) = -iS/k, so Im(-iS/k * gauss) = -S/k * gauss ... wait
    # Let me re-derive from Eq (19):
    # D = (e/m)^2 |Ek|^2 Im[1/(kv-w) + S*Delta(v,w/k)/k]
    # For gamma<0, S = 2*pi*i
    # S/k = 2*pi*i/k
    # The leading delta term: (2*pi*i/k) * delta(v - Omega/k)
    # Im[(2*pi*i/k) * delta(v-Omega/k)] = (2*pi/k) * delta(v-Omega/k)

    # Higher order terms are small for |gamma| << |Omega|
    S_contribution = (S / k) * gauss

    return term1 + S_contribution


def diffusion_coefficient_exact(v, omega_r, gamma, k):
    """Exact quasi-linear diffusion coefficient from the general solution.

    For all gamma, the full result using Eq (19) with S from Eq (17) gives:
    D_total = pi * delta(kv - Omega) in the limit gamma -> 0.

    For finite gamma:
    - gamma > 0: D = gamma / ((Omega - kv)^2 + gamma^2)  [Lorentzian]
    - gamma < 0: D = gamma / ((Omega - kv)^2 + gamma^2) + 2*pi/k * delta_reg(v - Omega/k)
    - gamma = 0: D = pi/k * delta(v - Omega/k)  [= pi * delta(kv - Omega)]

    Using Gaussian regularization for the delta function when gamma < 0.
    """
    Omega = omega_r
    lorentzian = gamma / ((Omega - k * v)**2 + gamma**2)

    if gamma > 0:
        return lorentzian
    elif abs(gamma) < 1e-30:
        # Limit: pi * delta(kv - Omega)
        # Represent as very narrow Gaussian
        sigma_kv = 1e-3 * abs(Omega)
        return pi * np.exp(-0.5 * ((k * v - Omega) / sigma_kv)**2) / (
            sigma_kv * np.sqrt(2 * pi))
    else:
        # gamma < 0: Lorentzian + correction
        # The correction makes the total continuous at gamma=0
        # Total should equal pi*delta(kv-Omega) as gamma->0-
        # Correction = 2*pi * |gamma| / ((Omega - kv)^2 + gamma^2) - lorentzian
        # = 2*pi * |gamma| / ((kv-Omega)^2 + gamma^2) + gamma/((Omega-kv)^2+gamma^2)
        # = (2*pi*|gamma| + gamma) / ((kv-Omega)^2 + gamma^2)
        # Since gamma < 0: |gamma| = -gamma, so 2*pi*(-gamma) + gamma = -(2*pi-1)*gamma
        # Hmm, that doesn't simplify nicely. Let me use the explicit form.

        # From Eq (19), the total D for gamma < 0 is:
        # Im[1/(kv-w)] + (2pi/k) * delta_series(v - Omega/k, gamma/k)
        # = gamma/((kv-Omega)^2 + gamma^2) + (2pi/k)*delta_reg(v-Omega/k)

        # The delta regularization width should be ~ |gamma/k| for consistency
        sigma_v = abs(gamma / k)
        v_res = Omega / k
        delta_reg = np.exp(-0.5 * ((v - v_res) / sigma_v)**2) / (
            sigma_v * np.sqrt(2 * pi))

        return lorentzian + (2 * pi / k) * delta_reg
