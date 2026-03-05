"""
Core physics module for runaway electron dynamics in tokamak fields.

Implements the relaxation (gyro-center) equations for momentum evolution
of runaway electrons with synchrotron radiation, following the physical
model in Wang, Qin & Liu, Phys. Plasmas 23, 062505 (2016).

All ODE integration uses normalized momentum ptilde = p / (m0*c) to
avoid numerical overflow from SI-unit momenta (~1e-21 kg*m/s).
"""

import numpy as np
from scipy.integrate import solve_ivp

# Physical constants (SI)
E_CHARGE = 1.602176634e-19      # electron charge [C]
M_E = 9.1093837015e-31          # electron rest mass [kg]
C_LIGHT = 2.99792458e8          # speed of light [m/s]
EPSILON_0 = 8.8541878128e-12    # vacuum permittivity [F/m]
PI = np.pi
M0C = M_E * C_LIGHT             # m0*c [kg*m/s]
M0C2 = M_E * C_LIGHT**2         # rest energy [J]
M0C2_MEV = 0.51099895           # rest energy [MeV]

# Synchrotron radiation: P_sync = C_SYNC * B^2 * p_perp_SI^2
C_SYNC = E_CHARGE**4 / (6 * PI * EPSILON_0 * M_E**4 * C_LIGHT**3)

# Curvature radiation: P_curv = e^2 * gamma^4 * v^4 / (6*pi*eps0 * c^3 * R^2)
#   In normalized form: P_curv = e^2 * c * ptilde_par^4 / (6*pi*eps0 * R^2)

# Drag coefficient D such that dp~/dt = alpha - D*ptilde_par:
#   D = P_total * gamma / (ptilde^2 * m0 * c^2)
#
# Sync contribution: D_sync = C_SYNC * B^2 * m0 * ptilde_perp^2 * gamma / ptilde^2
# Curv contribution: D_curv = [e^2/(6*pi*eps0*c*m0)] * ptilde_par^4 * gamma / (ptilde^2 * R^2)

TAU_S_COEFF = C_SYNC * M_E     # multiply by B0^2 to get tau_s [1/s]
TAU_C_COEFF = E_CHARGE**2 / (6 * PI * EPSILON_0 * C_LIGHT * M_E)  # divide by R_eff^2


def effective_curvature_radius(R0, q, r_orbit=0.0):
    """Compute effective orbit curvature radius."""
    kappa_tor_sq = 1.0 / R0**2
    if r_orbit > 0:
        kappa_hel = r_orbit / (r_orbit**2 + q**2 * R0**2)
        kappa_sq = kappa_tor_sq + kappa_hel**2
    else:
        kappa_sq = kappa_tor_sq
    return 1.0 / np.sqrt(kappa_sq)


def relaxation_rhs_normalized(t, y, alpha, tau_s, tau_c):
    """RHS of relaxation equations in normalized momentum units.

    y = [ptilde_par, ptilde_perp]

    dp~/dt_par  = alpha - D * ptilde_par
    dp~/dt_perp =       - D * ptilde_perp

    where D = (tau_s * ptilde_perp^2 + tau_c * ptilde_par^4) * gamma / ptilde^2
    """
    pp, pperp = y
    pperp = max(abs(pperp), 1e-30)

    p2 = pp**2 + pperp**2
    if p2 < 1e-60:
        return [alpha, 0.0]

    gamma = np.sqrt(1.0 + p2)
    D = (tau_s * pperp**2 + tau_c * pp**4) * gamma / p2

    return [alpha - D * pp, -D * pperp]


def solve_momentum_evolution(p_par0_m0c, p_perp0_m0c, E_l, B0, R0,
                             q=2.0, r_orbit=0.2, t_max=3.0, n_points=3000):
    """Solve relaxation equations for momentum evolution.

    Parameters
    ----------
    p_par0_m0c, p_perp0_m0c : float
        Initial momenta in units of m0*c.
    E_l : float
        Loop electric field [V/m].
    B0 : float
        Toroidal magnetic field [T].
    R0 : float
        Major radius [m].
    q : float
        Safety factor.
    r_orbit : float
        Typical orbit minor radius [m].
    t_max : float
        Maximum simulation time [s].
    n_points : int
        Number of output time points.

    Returns
    -------
    dict with keys: t, p_par, p_perp, gamma, energy_MeV, P_sync, P_curv
    """
    R_eff = effective_curvature_radius(R0, q, r_orbit)

    alpha = E_CHARGE * E_l / M0C              # [1/s]
    tau_s = TAU_S_COEFF * B0**2               # [1/s]
    tau_c = TAU_C_COEFF / R_eff**2            # [1/s]

    y0 = [p_par0_m0c, p_perp0_m0c]
    t_eval = np.linspace(0, t_max, n_points)

    sol = solve_ivp(
        relaxation_rhs_normalized, [0, t_max], y0,
        args=(alpha, tau_s, tau_c),
        t_eval=t_eval,
        method='RK45',
        rtol=1e-10, atol=1e-15,
        max_step=t_max / 100
    )

    p_par = sol.y[0]
    p_perp = np.abs(sol.y[1])
    gamma = np.sqrt(1 + p_par**2 + p_perp**2)
    energy_MeV = (gamma - 1) * M0C2_MEV

    # Radiation powers in Watts (for diagnostics)
    P_sync_arr = C_SYNC * B0**2 * (p_perp * M0C)**2
    P_curv_arr = (E_CHARGE**2 * C_LIGHT / (6 * PI * EPSILON_0)) * p_par**4 / R_eff**2

    return {
        't': sol.t,
        'p_par': p_par,
        'p_perp': p_perp,
        'gamma': gamma,
        'energy_MeV': energy_MeV,
        'P_sync': P_sync_arr,
        'P_curv': P_curv_arr,
    }


def find_energy_limit(E_l, B0, R0, q=2.0, r_orbit=0.2,
                      p_par0=5.0, p_perp0=1.0, t_max=5.0):
    """Find the synchrotron energy limit and balance time.

    Returns
    -------
    E_max_MeV : float
        Maximum kinetic energy [MeV].
    t_balance : float
        Time to reach 90% of energy limit [s].
    """
    result = solve_momentum_evolution(
        p_par0, p_perp0, E_l, B0, R0, q, r_orbit, t_max, n_points=5000)

    energy = result['energy_MeV']
    t = result['t']

    E_max = np.max(energy)

    threshold = 0.9 * E_max
    idx = np.searchsorted(energy, threshold)
    t_balance = t[min(idx, len(t) - 1)]

    return E_max, t_balance


def analytical_energy_limit(E_l, R_eff):
    """Analytical energy limit from curvature radiation balance.

    At equilibrium (p_perp -> 0, v -> c):
        e * E_l = e^2 * gamma^4 * c / (6*pi*eps0 * R_eff^2)
        gamma^4 = 6*pi*eps0 * E_l * R_eff^2 / e

    Returns gamma_max and E_max in MeV.
    """
    gamma4 = 6 * PI * EPSILON_0 * E_l * R_eff**2 / E_CHARGE
    gamma_max = gamma4**0.25
    E_max_MeV = (gamma_max - 1) * M0C2_MEV
    return gamma_max, E_max_MeV
