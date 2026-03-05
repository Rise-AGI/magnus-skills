"""
Shared plasma parameters and core functions for reproducing
Fulop et al., Phys. Plasmas 13, 062506 (2006).

All figures import from this module to avoid code duplication.
"""

import numpy as np
from scipy import special
from scipy.constants import e, m_e, epsilon_0, mu_0, c, pi, m_p
from scipy.integrate import quad
import warnings

warnings.filterwarnings("ignore", category=RuntimeWarning)


class PlasmaParams:
    """Plasma parameters for pure deuterium (m_i = 2m_p, Z=1)."""

    def __init__(self, n_e=5e19, B_T=2.0, T_eV=10.0, n_ratio=5e-3, Z=1,
                 E_par_norm=50, lnLambda=15):
        self.n_e = n_e
        self.B_T = B_T
        self.T_eV = T_eV
        self.n_ratio = n_ratio
        self.Z = Z
        self.E_par_norm = E_par_norm
        self.lnLambda = lnLambda

        self.n_r = n_ratio * n_e
        self.T_e = T_eV * e
        self.v_Te = np.sqrt(2 * self.T_e / m_e)
        self.m_i = 2 * m_p

        self.omega_pe = np.sqrt(n_e * e**2 / (m_e * epsilon_0))
        self.omega_pi = np.sqrt(n_e * e**2 / (self.m_i * epsilon_0))
        self.omega_ce = e * B_T / m_e
        self.omega_ci = e * B_T / self.m_i

        self.v_A = B_T / np.sqrt(mu_0 * n_e * self.m_i)

        self.tau = (4 * pi * epsilon_0**2 * m_e**2 * c**3) / (
            n_e * e**4 * lnLambda)

        self.c_Z = np.sqrt(3 * (Z + 5) / pi) * lnLambda

        self.alpha = (E_par_norm - 1) / (1 + Z)

        self.tau_ei = (3 * pi**1.5 * m_e**2 * self.v_Te**3 * epsilon_0**2 /
                       (n_e * Z**2 * e**4 * lnLambda))

        self.gamma_d = -1.5 / self.tau_ei

        self.omega_pr = np.sqrt(self.n_r * e**2 / (m_e * epsilon_0))


def omega0_dispersion(k, k_parallel, params):
    return k * params.v_A * np.sqrt(1 + k_parallel**2 * c**2 / params.omega_pi**2)


def omega_magnetosonic(k, params):
    return k * params.v_A


def omega_whistler(k, k_parallel, params):
    return k * k_parallel * params.v_A**2 / params.omega_ci


def stix_parameters(omega, params):
    wpe2 = params.omega_pe**2
    wpi2 = params.omega_pi**2
    wce = params.omega_ce
    wci = params.omega_ci
    w = omega
    w2 = w**2

    R = 1 - wpe2 / (w * (w + wce)) - wpi2 / (w * (w - wci))
    L = 1 - wpe2 / (w * (w - wce)) - wpi2 / (w * (w + wci))
    S = 0.5 * (R + L)
    D = 0.5 * (R - L)
    P = 1 - wpe2 / w2 - wpi2 / w2
    return S, D, P


def cold_dispersion_det(omega, k, theta, params):
    S, D_stix, P = stix_parameters(omega, params)
    k_par = k * np.cos(theta)
    n2 = (k * c / omega)**2
    npar2 = (k_par * c / omega)**2
    return (S - npar2) * (S - n2) - D_stix**2


def solve_cold_dispersion(k, theta, params, omega_guess=None):
    k_par = k * np.cos(theta)
    if omega_guess is None:
        omega_guess = omega0_dispersion(k, k_par, params)

    omega = omega_guess

    for _ in range(500):
        f = cold_dispersion_det(omega, k, theta, params)
        dw = omega * 1e-8 + 1.0
        fp = cold_dispersion_det(omega + dw, k, theta, params)
        fm = cold_dispersion_det(omega - dw, k, theta, params)
        df = (fp - fm) / (2 * dw)
        if abs(df) < 1e-50:
            break
        delta = -f / df

        max_step = 0.1 * abs(omega)
        if abs(delta) > max_step:
            delta = max_step * np.sign(delta)

        omega_new = omega + delta
        if omega_new <= params.omega_ci:
            omega_new = params.omega_ci * 1.1
        if omega_new > 0.5 * params.omega_ce:
            omega_new = 0.5 * params.omega_ce

        omega = omega_new
        if abs(delta / omega) < 1e-12:
            break
    return omega


def solve_cold_dispersion_sweep(k_array, theta, params):
    n = len(k_array)
    omega_out = np.zeros(n)

    for i in range(n):
        k = k_array[i]
        if i == 0:
            guess = omega0_dispersion(k, k * np.cos(theta), params)
        else:
            if i >= 2:
                dk = k_array[i] - k_array[i - 1]
                dk_prev = k_array[i - 1] - k_array[i - 2]
                if abs(dk_prev) > 0:
                    guess = omega_out[i - 1] + (omega_out[i - 1] - omega_out[i - 2]) * dk / dk_prev
                else:
                    guess = omega_out[i - 1]
            else:
                guess = omega_out[i - 1]

            if guess <= params.omega_ci:
                guess = omega0_dispersion(k, k * np.cos(theta), params)

        omega_out[i] = solve_cold_dispersion(k, theta, params, omega_guess=guess)
    return omega_out


def growth_rate_eq22(k, theta, params, omega=None):
    k_par = k * np.cos(theta)
    w0 = omega if omega is not None else omega0_dispersion(k, k_par, params)

    if abs(w0) < 1e-10:
        return 0.0

    if k_par * c <= w0:
        return 0.0

    exponent = -params.omega_ce / ((k_par * c - w0) * params.c_Z)
    if exponent < -500:
        return 0.0

    prefactor = (pi / (4 * params.c_Z)) * (params.omega_pr**2 / params.omega_pi**2)
    gamma_i = prefactor * (k**2 * params.v_A**2 / w0) * np.exp(exponent)
    return gamma_i


def growth_rate_eq20(k, theta, params, omega=None):
    k_par = k * np.cos(theta)
    k_perp = k * np.sin(theta)
    K_par = k_par * c / params.omega_ce
    K_perp = k_perp * c / params.omega_ce

    w0 = omega if omega is not None else omega0_dispersion(k, k_par, params)

    if abs(k_par * c) < 1e-10 or abs(w0) < 1e-10:
        return 0.0

    y = w0 / (k_par * c)

    if y >= 1.0:
        return 0.0

    a_n = params.alpha * K_par * (1 - y) / (-2)
    b_n = params.alpha * (1 - y) - K_par * (1 - y) - 1 / params.c_Z

    if abs(a_n) < 1e-30:
        return 0.0

    lam = K_perp**2 / (2 * a_n)

    exponent_C = -1.0 / (K_par * (1 - y) * params.c_Z)
    if exponent_C < -500:
        return 0.0

    C_hat = ((1 - y) * pi * params.alpha * params.omega_pr**2 /
             (2 * params.c_Z * params.omega_pi**2) *
             (k**2 * params.v_A**2 / w0**2) *
             (k_par**2 / k_perp**2) *
             np.exp(exponent_C))

    try:
        I0 = special.i0(lam)
        I1 = special.i1(lam)
    except (OverflowError, FloatingPointError):
        return 0.0

    denom = 4 * a_n**2
    term1 = K_perp**2 * K_par * (1 - y) * I0
    term2 = (2 * a_n * b_n + K_perp**2 * K_par * (1 - y)) * I1

    gamma_over_w0 = (C_hat / denom) * (term1 + term2) * np.exp(lam)
    gamma_i = gamma_over_w0 * w0

    return np.real(gamma_i)


def _eq19_integrand(p_perp, K_perp, K_par, y, a_n, b_n, params):
    z = K_perp * p_perp
    J1 = special.j1(z)
    exp_val = a_n * p_perp**2
    if exp_val < -500:
        return 0.0

    bracket = b_n + a_n * (1 - y) * K_par * p_perp**2 / (-1)
    return p_perp * J1**2 * np.exp(exp_val) * bracket


def growth_rate_eq19_numerical(k, theta, params, omega=None):
    k_par = k * np.cos(theta)
    k_perp = k * np.sin(theta)
    K_par = k_par * c / params.omega_ce
    K_perp = k_perp * c / params.omega_ce

    w0 = omega if omega is not None else omega0_dispersion(k, k_par, params)

    if abs(k_par * c) < 1e-10 or abs(w0) < 1e-10 or abs(k_perp) < 1e-10:
        return 0.0

    y = w0 / (k_par * c)

    if y >= 1.0:
        return 0.0

    a_n = params.alpha * K_par * (1 - y) / (-2)
    b_n = params.alpha * (1 - y) - K_par * (1 - y) - 1 / params.c_Z

    if abs(a_n) < 1e-30:
        return 0.0

    exponent_C = -1.0 / (K_par * (1 - y) * params.c_Z)
    if exponent_C < -500:
        return 0.0

    C_hat = ((1 - y) * pi * params.alpha * params.omega_pr**2 /
             (2 * params.c_Z * params.omega_pi**2) *
             (k**2 * params.v_A**2 / w0**2) *
             (k_par**2 / k_perp**2) *
             np.exp(exponent_C))

    if a_n < 0:
        p_perp_max = min(np.sqrt(10 / abs(a_n)) + 5, 200)
    else:
        p_perp_max = 50

    try:
        result, _ = quad(
            _eq19_integrand, 0, p_perp_max,
            args=(K_perp, K_par, y, a_n, b_n, params),
            limit=300, epsrel=1e-8, epsabs=1e-30
        )
    except Exception:
        return 0.0

    gamma_over_w0 = C_hat * result
    return gamma_over_w0 * w0


def growth_rate_numerical(k, theta, params):
    k_par = k * np.cos(theta)
    w0_anal = omega0_dispersion(k, k_par, params)
    omega_stix = solve_cold_dispersion(k, theta, params)

    if omega_stix > 2 * w0_anal or omega_stix < 0.2 * w0_anal:
        omega_stix = w0_anal

    return growth_rate_eq19_numerical(k, theta, params, omega=omega_stix)


def find_max_growth_rate(theta, params, k_min=10, k_max=2000, n_scan=200,
                         use_eq22=True):
    k_array = np.linspace(k_min, k_max, n_scan)
    gamma_array = np.zeros(n_scan)

    func = growth_rate_eq22 if use_eq22 else growth_rate_eq20

    for i, k in enumerate(k_array):
        gamma_array[i] = func(k, theta, params)

    idx = np.argmax(gamma_array)
    return k_array[idx], gamma_array[idx]


def find_max_growth_rate_numerical(theta, params, k_min=10, k_max=2000,
                                    n_scan=80):
    k_array = np.linspace(k_min, k_max, n_scan)
    gamma_array = np.zeros(n_scan)

    for i, k in enumerate(k_array):
        gamma_array[i] = growth_rate_numerical(k, theta, params)

    idx = np.argmax(gamma_array)
    return k_array[idx], gamma_array[idx]
