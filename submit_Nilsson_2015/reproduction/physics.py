"""
Core physics module for reproducing Nilsson et al. (2015)
"Kinetic modelling of runaway electron avalanches in tokamak plasmas"
Plasma Phys. Control. Fusion 57 (2015) 095006

Key equations implemented:
- Critical electric field Ec (Eq. 1)
- Relativistic collision time tau (Eq. 5)
- Analytic avalanche growth rate (Eq. 6, Rosenbluth)
- Dreicer generation rate (Eq. 20, Connor-Hastie)
- Toroidicity-dependent Dreicer fit (Fig. 7)
- Toroidicity-dependent avalanche rate (Eq. A.4)
"""

import numpy as np
from scipy.integrate import solve_ivp

# Physical constants (SI)
e = 1.602176634e-19       # elementary charge [C]
m_e = 9.1093837015e-31    # electron rest mass [kg]
c = 2.99792458e8          # speed of light [m/s]
eps0 = 8.8541878128e-12   # vacuum permittivity [F/m]
k_B = 1.380649e-23        # Boltzmann constant [J/K]


def critical_field(n_e, ln_Lambda):
    """Critical electric field Ec (Eq. 1)."""
    return n_e * e**3 * ln_Lambda / (4 * np.pi * eps0**2 * m_e * c**2)


def dreicer_field(n_e, T_eV, ln_Lambda):
    """Dreicer field ED = (c/vth)^2 * Ec."""
    T_J = T_eV * e
    v_th = np.sqrt(T_J / m_e)
    Ec = critical_field(n_e, ln_Lambda)
    return (c / v_th)**2 * Ec


def collision_time_rel(n_e, ln_Lambda):
    """Relativistic collision time tau (Eq. 5)."""
    return 4 * np.pi * eps0**2 * m_e**2 * c**3 / (n_e * e**4 * ln_Lambda)


def thermal_collision_freq(n_e, T_eV, ln_Lambda):
    """Thermal collision frequency nu_th = 1/tau(v_th).
    tau(v) = 4*pi*eps0^2*m_e^2*v^3 / (n_e*e^4*ln_Lambda)
    """
    T_J = T_eV * e
    v_th = np.sqrt(T_J / m_e)
    tau_th = 4 * np.pi * eps0**2 * m_e**2 * v_th**3 / (n_e * e**4 * ln_Lambda)
    return 1.0 / tau_th


def avalanche_growth_rate(E_over_Ec, n_e, ln_Lambda):
    """Analytic avalanche multiplication factor gamma_A_bar (from Eq. 6).
    (1/nr)(dnr/dt) = gamma_A_bar * nr  where
    gamma_A_bar = 1/(2*tau*ln_Lambda) * (E/Ec - 1)
    Returns rate in [1/s] (per unit nr).
    """
    tau = collision_time_rel(n_e, ln_Lambda)
    return (1.0 / (2 * tau * ln_Lambda)) * (E_over_Ec - 1)


def dreicer_rate(E_over_Ec, n_e, T_eV, ln_Lambda, Z_eff=1):
    """Dreicer generation rate gamma_D (Connor-Hastie, Ref [6]).
    Uses the standard parameterization matching Kulsrud numerical results.

    gamma_D = ne * nu_th * C_R * (ED/E)^h * exp(-ED/(4E) - sqrt((1+Z)*ED/(2E)))
    where h = 3(1+Z)/16

    Returns: total generation rate dnr/dt [m^-3 s^-1]
    """
    T_J = T_eV * e
    v_th = np.sqrt(T_J / m_e)
    Ec = critical_field(n_e, ln_Lambda)
    ED_over_Ec = (c / v_th)**2
    E_over_ED = E_over_Ec / ED_over_Ec
    ED_over_E = 1.0 / E_over_ED

    nu_th = thermal_collision_freq(n_e, T_eV, ln_Lambda)

    if E_over_ED <= 0:
        return 0.0

    h = 3.0 * (1 + Z_eff) / 16.0
    C_R = 1.0  # prefactor

    exponent = -ED_over_E / 4.0 - np.sqrt((1 + Z_eff) * ED_over_E / 2.0)
    if exponent < -500:
        return 0.0

    rate = C_R * n_e * nu_th * ED_over_E**h * np.exp(exponent)
    return rate


def dreicer_rate_normalized(E_over_Ec, T_eV, ln_Lambda, Z_eff=1):
    """Dreicer generation rate normalized to thermal collision frequency.
    Returns gamma_D / (ne * nu_th) (dimensionless).
    """
    T_J = T_eV * e
    v_th = np.sqrt(T_J / m_e)
    ED_over_Ec = (c / v_th)**2
    E_over_ED = E_over_Ec / ED_over_Ec
    ED_over_E = 1.0 / E_over_ED

    if E_over_ED <= 0:
        return 0.0

    h = 3.0 * (1 + Z_eff) / 16.0
    exponent = -ED_over_E / 4.0 - np.sqrt((1 + Z_eff) * ED_over_E / 2.0)
    if exponent < -500:
        return 0.0

    return ED_over_E**h * np.exp(exponent)


def avalanche_rate_normalized(E_over_Ec, T_eV, ln_Lambda):
    """Avalanche multiplication factor normalized to thermal collision freq.
    Returns gamma_A_bar / nu_th (dimensionless, per unit nr).
    gamma_A_bar = 1/(2*tau*ln_Lambda) * (E/Ec - 1)
    nu_th = 1/tau_th = (c/v_th)^3 / tau
    So gamma_A_bar / nu_th = (v_th/c)^3 / (2*ln_Lambda) * (E/Ec - 1)
    """
    T_J = T_eV * e
    v_th = np.sqrt(T_J / m_e)
    beta_th = v_th / c
    return beta_th**3 / (2 * ln_Lambda) * (E_over_Ec - 1)


def solve_runaway_evolution(E_over_Ec, T_eV, n_e, ln_Lambda,
                            t_max_tau_th, n_points=2000,
                            include_avalanche=True, epsilon=0.0):
    """Solve ODE for runaway electron fraction nr/ntot vs time.

    dnr/dt = (ntot - nr) * (gamma_D + nr * gamma_A_bar)

    Time is normalized to thermal collision time tau_th.
    epsilon: inverse aspect ratio (0 = cylindrical)

    Returns: t_normalized, nr_over_ntot
    """
    nu_th = thermal_collision_freq(n_e, T_eV, ln_Lambda)
    tau_th = 1.0 / nu_th

    # Dreicer rate (already total rate, divide by ne to get per-density)
    gamma_D_total = dreicer_rate(E_over_Ec, n_e, T_eV, ln_Lambda)
    gamma_D_per_ne = gamma_D_total / n_e if n_e > 0 else 0.0

    # Apply toroidicity correction to Dreicer
    if epsilon > 0:
        gamma_D_per_ne *= dreicer_toroidicity_factor(epsilon)

    # Avalanche multiplication factor
    gamma_A_bar = avalanche_growth_rate(E_over_Ec, n_e, ln_Lambda) if include_avalanche else 0.0

    # Apply toroidicity correction to avalanche
    if epsilon > 0 and include_avalanche:
        gamma_A_bar *= avalanche_toroidicity_factor(epsilon, E_over_Ec)

    # gamma_A_bar = (1/(2*tau*lnL)) * (E/Ec - 1) is the exponential growth rate
    # In the Eq. 10 formulation: (1/(ntot-nr))*dnr/dt = gamma_D + nr * gamma_A_bar / ne
    # So: dnr/dt = ne_local * gamma_D_per_ne + nr * gamma_A_bar * ne_local / ne
    gamma_A_normalized = gamma_A_bar / n_e if n_e > 0 else 0.0

    def rhs(t, y):
        nr = y[0]
        ne_local = n_e - nr
        if ne_local < 0:
            ne_local = 0
        dnr_dt = ne_local * (gamma_D_per_ne + nr * gamma_A_normalized)
        return [dnr_dt]

    # Initial condition: tiny seed
    nr0 = n_e * 1e-20
    t_span = (0, t_max_tau_th * tau_th)
    t_eval = np.linspace(0, t_max_tau_th * tau_th, n_points)

    sol = solve_ivp(rhs, t_span, [nr0], t_eval=t_eval,
                    method='RK45', rtol=1e-10, atol=1e-30, max_step=tau_th*t_max_tau_th/1000)

    t_norm = sol.t / tau_th
    nr_frac = sol.y[0] / n_e

    return t_norm, nr_frac


def dreicer_toroidicity_factor(epsilon):
    """Fit from Fig. 7: gamma_D/gamma_D,cyl = 1 - 1.2*sqrt(2*epsilon/(1+epsilon))."""
    val = 1.0 - 1.2 * np.sqrt(2 * epsilon / (1 + epsilon))
    return np.maximum(val, 0.0)


def avalanche_toroidicity_factor(epsilon, E_over_Ec):
    """Toroidicity-dependent avalanche growth rate factor (Eq. A.4).
    Returns ratio gamma_A(epsilon) / gamma_A(0).
    """
    if np.isscalar(epsilon):
        return _avalanche_toro_scalar(epsilon, E_over_Ec)
    return np.array([_avalanche_toro_scalar(eps, E_over_Ec) for eps in epsilon])


def _avalanche_toro_scalar(epsilon, E_over_Ec):
    if epsilon <= 0:
        return 1.0
    if epsilon >= 1:
        return 0.0

    # Check if trapping matters: epsilon*E/Ec*(1-epsilon)^(-2) > 1/4
    ratio = epsilon * E_over_Ec / (1 - epsilon)**2

    if ratio <= 0.25:
        # p_c is always the lower limit, no trapping effect
        return 1.0

    # theta_bound from Eq. A.2
    cos_val = (1 - epsilon)**2 / (2 * epsilon * E_over_Ec) - 1
    cos_val = np.clip(cos_val, -1, 1)
    theta_b = np.arccos(cos_val)

    if theta_b <= 0:
        return 1.0
    if theta_b >= np.pi:
        # Full trapping limit
        return (1 - epsilon)**2 / (np.pi * np.sqrt(epsilon * E_over_Ec))

    sin_tb = np.sin(theta_b)

    # From Eq. A.4:
    # factor = (1 - theta_b/pi - epsilon*sin(theta_b)/pi) +
    #          (1-epsilon)^2 * Ec/(2*epsilon*pi*E) * (sqrt(1-epsilon)*sqrt(4*epsilon*E/Ec - (1-epsilon)) + epsilon*theta_b)
    term1 = 1 - theta_b / np.pi - epsilon * sin_tb / np.pi

    inner_sqrt_arg = 4 * epsilon * E_over_Ec - (1 - epsilon)**2
    if inner_sqrt_arg < 0:
        inner_sqrt = 0
    else:
        inner_sqrt = np.sqrt(inner_sqrt_arg)

    # tan(theta_b/2) term contribution
    term2 = ((1 - epsilon)**2 / (2 * epsilon * np.pi * E_over_Ec)) * \
            (np.sqrt(1 - epsilon) * inner_sqrt + epsilon * theta_b)

    factor = term1 + term2

    # Normalize: at epsilon=0, factor should be 1 (=E/Ec / (E/Ec))
    # The factor multiplies (E/Ec) in the full growth rate expression
    # Since at epsilon=0, the growth rate is proportional to (E/Ec - 1) ~ E/Ec for large E/Ec
    # The returned factor is already the ratio gamma_A(eps)/gamma_A(0)
    return np.clip(factor, 0, 1)


def growth_rate_ratio_analytic(E_over_Ec, T_eV, ln_Lambda, nr_over_ne=0.01, Z_eff=1):
    """Ratio gamma_A/gamma_D from Eq. 21 (modified for consistency).

    Uses: gamma_A = nr * (1/(2*tau*lnL)) * (E/Ec - 1)
          gamma_D = ne * nu_th * (ED/E)^h * exp(-ED/(4E) - sqrt((1+Z)*ED/(2E)))
    """
    T_J = T_eV * e
    v_th = np.sqrt(T_J / m_e)
    beta_th = v_th / c
    ED_over_Ec = (c / v_th)**2
    E_over_ED = E_over_Ec / ED_over_Ec
    ED_over_E = 1.0 / E_over_ED

    if E_over_ED <= 0 or E_over_Ec <= 1:
        return 0.0

    h = 3.0 * (1 + Z_eff) / 16.0

    # gamma_A_bar / ne = (1/(2*tau*lnL)) * (E/Ec - 1) / ne
    # = beta_th^3 / (2*lnL) * nu_th * (E/Ec - 1)
    # gamma_D / ne = nu_th * (ED/E)^h * exp(...)
    # Ratio = nr * gamma_A_bar/ne / (gamma_D/ne)
    #       = nr * beta_th^3 * (E/Ec-1) / (2*lnL) / ((ED/E)^h * exp(...))

    exponent_D = -ED_over_E / 4.0 - np.sqrt((1 + Z_eff) * ED_over_E / 2.0)
    if exponent_D < -500:
        return np.inf  # Dreicer negligible, avalanche dominates

    dreicer_norm = ED_over_E**h * np.exp(exponent_D)
    if dreicer_norm <= 0:
        return np.inf

    avalanche_norm = beta_th**3 / (2 * ln_Lambda) * (E_over_Ec - 1) * nr_over_ne

    return avalanche_norm / dreicer_norm


def compute_avalanche_fraction(E_over_Ec, T_eV, n_e, ln_Lambda,
                               target_frac=0.01):
    """Compute fraction of runaways from knock-on collisions when
    nr/ntot reaches target_frac.

    Returns: n_A/n_r (fraction from avalanche)
    """
    nu_th = thermal_collision_freq(n_e, T_eV, ln_Lambda)
    tau_th = 1.0 / nu_th

    gamma_D_total = dreicer_rate(E_over_Ec, n_e, T_eV, ln_Lambda)
    gamma_D_per_ne = gamma_D_total / n_e if n_e > 0 else 0.0
    gamma_A_bar = avalanche_growth_rate(E_over_Ec, n_e, ln_Lambda)
    gamma_A_normalized = gamma_A_bar / n_e if n_e > 0 else 0.0

    nr_target = target_frac * n_e

    # Track both Dreicer and avalanche contributions
    nr = n_e * 1e-20
    nr_dreicer = 0.0
    nr_avalanche = 0.0
    dt = tau_th * 0.01
    max_steps = int(1e8)

    for _ in range(max_steps):
        ne_local = n_e - nr
        if ne_local < 0:
            ne_local = 0

        d_dreicer = ne_local * gamma_D_per_ne * dt
        d_avalanche = ne_local * nr * gamma_A_normalized * dt

        nr_dreicer += d_dreicer
        nr_avalanche += d_avalanche
        nr += d_dreicer + d_avalanche

        if nr >= nr_target:
            break

        # Adaptive timestep
        total_rate = ne_local * (gamma_D_per_ne + nr * gamma_A_normalized)
        if total_rate > 0:
            dt = min(0.01 * nr / total_rate, tau_th * 1.0)
            dt = max(dt, tau_th * 1e-6)

    if nr_dreicer + nr_avalanche > 0:
        return nr_avalanche / (nr_dreicer + nr_avalanche)
    return 0.0


def time_to_fraction(E_over_Ec, T_eV, n_e, ln_Lambda,
                     target_frac=0.01, include_avalanche=True,
                     max_tau_th=1e12):
    """Compute time (in tau_th units) for nr/ntot to reach target_frac.
    Returns np.inf if not reached within max_tau_th.
    """
    nu_th = thermal_collision_freq(n_e, T_eV, ln_Lambda)
    tau_th = 1.0 / nu_th

    gamma_D_total = dreicer_rate(E_over_Ec, n_e, T_eV, ln_Lambda)
    gamma_D_per_ne = gamma_D_total / n_e if n_e > 0 else 0.0
    gamma_A_bar = avalanche_growth_rate(E_over_Ec, n_e, ln_Lambda) if include_avalanche else 0.0
    gamma_A_normalized = gamma_A_bar / n_e if n_e > 0 else 0.0

    nr_target = target_frac * n_e
    nr = n_e * 1e-20
    t = 0.0
    dt = tau_th * 0.01
    t_max = max_tau_th * tau_th

    while t < t_max:
        ne_local = n_e - nr
        if ne_local < 0:
            ne_local = 0

        total_rate = ne_local * (gamma_D_per_ne + nr * gamma_A_normalized)

        if total_rate > 0:
            dt = min(0.005 * max(nr, n_e * 1e-15) / total_rate, tau_th * 10.0)
            dt = max(dt, tau_th * 1e-8)

        dnr = total_rate * dt
        nr += dnr
        t += dt

        if nr >= nr_target:
            return t / tau_th

    return np.inf
