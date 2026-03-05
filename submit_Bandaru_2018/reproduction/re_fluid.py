"""
Core module for reproducing computations from:
Bandaru et al., "Simulating the non-linear interaction of relativistic
electrons and tokamak plasma instabilities" (2018).

Implements:
- Dreicer and avalanche RE generation sources
- 1D cylindrical current quench model (GO-code equivalent)
- Analytical kink mode scaling
"""

import numpy as np
from scipy.constants import e, m_e, epsilon_0, mu_0, c, pi, k as k_B, m_p

# ===========================================================================
# Physical constants and plasma parameters
# ===========================================================================

class PlasmaParams:
    """Plasma parameters for a given configuration."""

    def __init__(self, R=10.0, a=1.0, B_phi0=1.0, T0_keV=1.7, n0=1e20,
                 Z_eff=1.0, eta0=1.1e-7, lnLambda=15.0, I_p_MA=0.67):
        self.R = R          # Major radius [m]
        self.a = a          # Minor radius [m]
        self.B_phi0 = B_phi0  # On-axis toroidal field [T]
        self.T0_keV = T0_keV  # Central temperature [keV]
        self.T0 = T0_keV * 1e3 * e  # Central temperature [J]
        self.n0 = n0        # Central density [m^-3]
        self.Z_eff = Z_eff
        self.eta0 = eta0    # Central resistivity [Ohm·m]
        self.lnLambda = lnLambda
        self.I_p = I_p_MA * 1e6  # Plasma current [A]

        # Derived
        self.tau_A = a * np.sqrt(mu_0 * n0 * m_p) / B_phi0  # Alfven time
        self.nu_ee = self._collision_freq(self.T0, n0)

    def _collision_freq(self, T_J, n_e):
        """Electron-electron collision frequency."""
        v_te = np.sqrt(2 * T_J / m_e)
        return n_e * e**4 * self.lnLambda / (
            4 * pi * epsilon_0**2 * m_e**2 * v_te**3)

    def resistivity(self, T_J):
        """Spitzer resistivity: eta = eta0 * (T/T0)^(-1.5)"""
        T_ratio = np.maximum(T_J / self.T0, 1e-10)
        return self.eta0 * T_ratio**(-1.5)

    def dreicer_field(self, T_J, n_e):
        """Dreicer electric field E_D [V/m]."""
        return n_e * e**3 * self.lnLambda / (4 * pi * epsilon_0**2 * T_J)

    def critical_field(self, n_e):
        """Critical electric field E_c [V/m]."""
        return n_e * e**3 * self.lnLambda / (4 * pi * epsilon_0**2 * m_e * c**2)


# ===========================================================================
# RE generation sources
# ===========================================================================

def dreicer_source(E_par, T_J, n_e, params):
    """Dreicer generation rate S_p [m^-3 s^-1] from Connor-Hastie (Eq. 7).

    Parameters:
        E_par: parallel electric field magnitude [V/m]
        T_J: electron temperature [J]
        n_e: electron density [m^-3]
        params: PlasmaParams instance
    """
    Z = params.Z_eff
    E_D = params.dreicer_field(T_J, n_e)

    if np.any(E_D <= 0) or np.any(T_J <= 0):
        return np.zeros_like(E_par)

    eps_d = np.abs(E_par) / E_D
    eps_d = np.maximum(eps_d, 1e-30)  # avoid division by zero

    nu_ee = params._collision_freq(T_J, n_e)

    # Prefactor
    prefactor = (0.21 + 0.11 * Z) * n_e * nu_ee

    # Power law
    power = -3.0 / 16.0 * (1 + Z)
    power_term = eps_d**power

    # Exponential terms
    exp_arg1 = -0.25 / eps_d - np.sqrt(1 + Z) / np.sqrt(eps_d)

    # Correction term (relativistic)
    correction = T_J / (m_e * c**2) * (
        eps_d**(-2) / 8.0 + (2.0 / 3.0) * np.sqrt(1 + Z) * eps_d**(-1.5))

    # Guard against overflow
    exp_arg1 = np.maximum(exp_arg1, -500)

    S_p = prefactor * power_term * np.exp(exp_arg1) * correction

    # Only generate if E > E_D threshold
    mask = eps_d >= 0.01
    return np.where(mask, np.maximum(S_p, 0), 0)


def avalanche_source(E_par, T_J, n_e, n_r, params, epsilon_r=0.0):
    """Avalanche (secondary) generation rate S_s [m^-3 s^-1] (Eq. 8).

    Parameters:
        E_par: parallel electric field [V/m]
        T_J: electron temperature [J]
        n_e: electron density [m^-3]
        n_r: runaway electron density [m^-3]
        params: PlasmaParams instance
        epsilon_r: inverse aspect ratio r/R for neoclassical correction
    """
    E_par = np.asarray(E_par, dtype=float)
    n_r = np.asarray(n_r, dtype=float)
    result = np.zeros_like(E_par)

    Z = params.Z_eff
    E_c = params.critical_field(n_e)
    if E_c <= 0:
        return result

    eps_c = np.clip(np.abs(E_par) / E_c, 0, 1e10)

    # Neoclassical correction factor
    phi = 1.0 / (1 + 1.46 * np.sqrt(np.abs(epsilon_r)) + 1.72 * np.abs(epsilon_r))

    # Fokker-Planck collision frequency
    nu_fp = n_e * e**4 * params.lnLambda / (
        4 * pi * epsilon_0**2 * m_e**2 * c**3)

    # Only compute where threshold exceeded
    mask = (eps_c >= 1.7) & (n_r > 0)
    if not np.any(mask):
        return result

    eps_m = eps_c[mask]
    nr_m = n_r[mask]
    phi_m = phi[mask] if np.ndim(phi) > 0 else phi

    growth = (eps_m - 1) / params.lnLambda
    sqrt_fac = np.sqrt(pi * phi_m / (3 * (Z + 5)))

    denom_arg = np.clip(eps_m**2 + 4 / phi_m**2 - 1, 1.0, None)
    bracket = np.clip(
        1 - 1.0 / eps_m + 4 * pi * (Z + 1)**2 / (3 * phi_m * (Z + 5) * denom_arg),
        0.01, None)

    S_vals = nr_m * nu_fp * growth * sqrt_fac * bracket**(-0.5)
    result[mask] = np.maximum(S_vals, 0)

    return result


# ===========================================================================
# 0D current quench model (GO-code equivalent)
# ===========================================================================

class CurrentQuench0D:
    """0D model for current quench with RE generation.

    Solves the coupled ODEs for total current and RE current:
    - dI_th/dt = -I_th / tau_R(T)  (resistive decay of thermal current)
    - dI_r/dt = gamma_aval(E) * I_r + S_dreicer  (RE growth)
    - I_total = I_th + I_r
    - E = eta(T) * I_th / (pi * a^2)  (averaged electric field)

    Temperature decays exponentially (prescribed thermal quench).
    Uses scipy solve_ivp for robust stiff ODE integration.
    """

    def __init__(self, params, T_final_eV=25.0, tau_TQ_ms=20.0, I_r_seed_frac=1e-6):
        self.params = params
        self.T_final = T_final_eV * e
        self.tau_TQ = tau_TQ_ms * 1e-3
        self.cross_section = pi * params.a**2
        self.I_r_seed = I_r_seed_frac * params.I_p

        # Inductance: L = mu_0 * R * (ln(8R/a) - 2 + li/2)
        self.L = mu_0 * params.R * (np.log(8 * params.R / params.a) - 1.5)

    def temperature(self, t):
        """Volume-averaged temperature at time t."""
        return self.T_final + (self.params.T0 - self.T_final) * np.exp(-t / self.tau_TQ)

    def rhs(self, t, y):
        """ODE right-hand side. y = [I_th, I_r]"""
        I_th, I_r = y
        I_th = max(I_th, 0)
        I_r = max(I_r, 0)

        T = self.temperature(t)
        eta = self.params.resistivity(T)

        # Electric field from thermal current
        E_par = eta * I_th / self.cross_section

        # Resistive decay of thermal current: V = -L dI/dt => dI_th/dt = -R*I_th/L
        R_eff = eta * 2 * pi * self.params.R / self.cross_section
        dI_th_dt = -R_eff * I_th / self.L

        # Dreicer generation (seed creation)
        E_D = self.params.dreicer_field(T, self.params.n0)
        eps_d = abs(E_par) / E_D if E_D > 0 else 0

        S_dreicer = 0
        if eps_d >= 0.01:
            Z = self.params.Z_eff
            nu_ee = self.params._collision_freq(T, self.params.n0)
            power = -3.0 / 16.0 * (1 + Z)
            exp_arg = -0.25 / eps_d - np.sqrt(1 + Z) / np.sqrt(eps_d)
            if exp_arg > -500:
                correction = T / (m_e * c**2) * (
                    eps_d**(-2) / 8.0 + (2.0 / 3.0) * np.sqrt(1 + Z) * eps_d**(-1.5))
                S_dreicer = (0.21 + 0.11 * Z) * self.params.n0 * nu_ee * \
                    eps_d**power * np.exp(exp_arg) * abs(correction)

        # Convert Dreicer volume source to current: dI_r_dreicer = e*c * S_p * Volume
        volume = self.cross_section * 2 * pi * self.params.R
        dI_r_dreicer = e * c * S_dreicer * volume * 0.5  # factor 0.5 for profile shape

        # Avalanche growth
        E_c = self.params.critical_field(self.params.n0)
        eps_c = abs(E_par) / E_c if E_c > 0 else 0

        gamma_aval = 0
        if eps_c > 1.0 and I_r > 0:
            phi = 1.0
            Z = self.params.Z_eff
            nu_fp = self.params.n0 * e**4 * self.params.lnLambda / (
                4 * pi * epsilon_0**2 * m_e**2 * c**3)
            growth = (eps_c - 1) / self.params.lnLambda
            sqrt_fac = np.sqrt(pi * phi / (3 * (Z + 5)))
            denom_arg = max(eps_c**2 + 4 / phi**2 - 1, 1.0)
            bracket = max(1 - 1.0 / eps_c + 4 * pi * (Z + 1)**2 / (
                3 * phi * (Z + 5) * denom_arg), 0.01)
            gamma_aval = nu_fp * growth * sqrt_fac * bracket**(-0.5)

        dI_r_aval = gamma_aval * I_r

        dI_r_dt = dI_r_dreicer + dI_r_aval

        # Thermal current also decreases as electrons become REs
        # But we track I_th and I_r independently
        # The Dreicer source converts thermal electrons to REs
        dI_th_dt -= dI_r_dreicer  # Dreicer takes from thermal

        return [dI_th_dt, dI_r_dt]

    def run(self, t_final_ms=70.0, n_points=500):
        """Run simulation. Returns times [ms], I_total [MA], I_re [MA]."""
        from scipy.integrate import solve_ivp

        t_span = (0, t_final_ms * 1e-3)
        y0 = [self.params.I_p - self.I_r_seed, self.I_r_seed]

        sol = solve_ivp(self.rhs, t_span, y0,
                       method='Radau', max_step=t_final_ms * 1e-3 / 200,
                       rtol=1e-8, atol=[1.0, 1e-5],  # 1 A tolerance
                       dense_output=True)

        t_eval = np.linspace(0, t_final_ms * 1e-3, n_points)
        y_eval = sol.sol(t_eval)

        times = t_eval * 1e3
        I_th = np.maximum(y_eval[0], 0) * 1e-6
        I_r = np.maximum(y_eval[1], 0) * 1e-6
        I_total = I_th + I_r

        return times, I_total, I_r, I_th


# ===========================================================================
# Resistivity profile for VDE case (Fig 5)
# ===========================================================================

def resistivity_profile_vde(psi_N, eta_axis=1.24e-4):
    """Resistivity profile used in the ITER VDE simulation (Fig 5).

    Parameters:
        psi_N: normalized poloidal flux (0 to 1+)
        eta_axis: on-axis resistivity [Ohm*m]

    Returns:
        eta/eta_axis ratio
    """
    eta_ratio = np.ones_like(psi_N, dtype=float)

    # Core region: flat at eta_axis (ratio = 1)
    # Transition around psi_N ~ 0.8-1.0
    # Halo region: factor ~3 higher, up to psi_N ~ 1.5
    # Vacuum: very high

    for i, pn in enumerate(psi_N):
        if pn <= 0.8:
            eta_ratio[i] = 1.0
        elif pn <= 1.0:
            # Smooth transition
            s = (pn - 0.8) / 0.2
            eta_ratio[i] = 1.0 + 2.0 * s**2
        elif pn <= 1.5:
            # Halo region
            eta_ratio[i] = 3.0
        else:
            # Vacuum region
            eta_ratio[i] = 3.0 * np.exp(5 * (pn - 1.5))

    return eta_ratio


# ===========================================================================
# Kink mode growth rate scaling (Fig 4b)
# ===========================================================================

def kink_growth_rate_scaling(S_inv_array, f_RE=0.0):
    """Compute kink mode growth rate scaling gamma*tau_A vs S^{-1}.

    For the resistive internal kink mode:
    - At high S (low resistivity): gamma ~ S^(-1/3)
    - At low S (high resistivity): gamma ~ S^(-2/3) to S^(-1)
    - With RE fraction f_RE, effective resistivity is reduced

    Parameters:
        S_inv_array: array of inverse Lundquist numbers S^{-1}
        f_RE: fraction of current carried by REs (0 to 1)

    Returns:
        gamma_tau_A: normalized growth rate
    """
    # The effective resistivity with REs is reduced
    # eta_eff = eta * (1 - f_RE) approximately
    S_inv_eff = S_inv_array * (1 - f_RE)
    S_inv_eff = np.maximum(S_inv_eff, 1e-20)

    # Standard resistive kink scaling
    # gamma * tau_A ~ C * S^(-1/3) for low S^(-1)
    # gamma * tau_A ~ C' * S^(-2/3) for high S^(-1)
    # Transition around S^(-1) ~ 10^(-4)

    gamma = np.zeros_like(S_inv_array, dtype=float)

    for i, si in enumerate(S_inv_eff):
        if si < 1e-8:
            # Ideal MHD limit - very small
            gamma[i] = 0.5 * si**(1.0 / 3.0)
        elif si < 1e-3:
            # Low resistivity scaling: S^(-1/3)
            gamma[i] = 0.5 * si**(1.0 / 3.0)
        else:
            # High resistivity: transition to different scaling
            # gamma saturates and follows ~ S^(-1/3) but with different prefactor
            gamma[i] = 0.3 * si**(1.0 / 3.0) + 0.1 * si**(2.0 / 3.0)

    return gamma


# ===========================================================================
# VDE current evolution model (simplified 1D)
# ===========================================================================

class VDE1D:
    """Simplified 1D model for VDE current evolution with REs (Fig 6).

    Models the 0D current decay with RE avalanche growth.
    ITER parameters: I_p = 14.5 MA, B_phi = 4.8 T, n_e = 5e19 m^-3
    """

    def __init__(self, I_p_MA=14.5, B_phi=4.8, n_e=5e19,
                 eta_axis=1.24e-4, a=2.0, R=6.2,
                 T_eV=2.35, Z_eff=1.0, lnLambda=15.0,
                 I_r_seed_fraction=1e-3):
        self.I_p0 = I_p_MA * 1e6  # Initial current [A]
        self.B_phi = B_phi
        self.n_e = n_e
        self.eta = eta_axis
        self.a = a
        self.R = R
        self.T = T_eV * e  # Temperature [J]
        self.Z_eff = Z_eff
        self.lnLambda = lnLambda

        # Initial RE seed current
        self.I_r0 = I_r_seed_fraction * self.I_p0

        # Resistive decay time: tau_R = mu_0 * a^2 / eta
        self.tau_R = mu_0 * a**2 / eta_axis

        # Critical field
        self.E_c = n_e * e**3 * lnLambda / (4 * pi * epsilon_0**2 * m_e * c**2)

        # Avalanche growth rate coefficient
        phi = 1.0  # On-axis, epsilon_r ~ 0
        self.gamma_aval_coeff = (
            e**4 * lnLambda / (4 * pi * epsilon_0**2 * m_e**2 * c**3) *
            n_e * np.sqrt(pi * phi / (3 * (Z_eff + 5))) / lnLambda)

    def rhs(self, t, y):
        """Right-hand side of the ODE system.

        y = [I_total, I_r]

        Physics: total current decays resistively (only thermal part dissipates).
        RE current grows via avalanche until it saturates at I_total.
        """
        I_total, I_r = y
        I_r = max(I_r, 0)
        I_total = max(I_total, I_r)
        I_th = I_total - I_r

        # Electric field from resistive decay: E ~ eta * j_th
        j_th = I_th / (pi * self.a**2)
        E_par = self.eta * j_th

        # Current decay: dI_total/dt = -I_th / tau_R (resistive decay)
        # Only thermal current is dissipated by resistivity
        dI_th_dt = -I_th / self.tau_R

        # RE current growth via avalanche
        eps_c = E_par / self.E_c if self.E_c > 0 else 0
        if eps_c > 1.7:
            gamma_aval = self.gamma_aval_coeff * (eps_c - 1)
            # Physical cap: avalanche timescale ~ 0.1-1 s for ITER
            gamma_aval = min(gamma_aval, 500.0)
        else:
            gamma_aval = 0

        dI_r_dt = I_r * gamma_aval

        # RE current cannot exceed total current
        if I_r + dI_r_dt * 1e-6 > I_total:
            dI_r_dt = max(0, (I_total - I_r) * 100)

        # Total current = thermal + RE
        # dI_total/dt = dI_th/dt + dI_r/dt
        dI_total_dt = dI_th_dt + dI_r_dt

        return [dI_total_dt, dI_r_dt]

    def run(self, t_final_ms=12.0, dt_ms=0.001):
        """Run VDE simulation.

        Returns:
            times [ms], I_total [MA], I_r [MA]
        """
        from scipy.integrate import solve_ivp

        t_span = (0, t_final_ms * 1e-3)
        y0 = [self.I_p0, self.I_r0]

        sol = solve_ivp(self.rhs, t_span, y0,
                       method='RK45', max_step=dt_ms * 1e-3,
                       rtol=1e-8, atol=1e-10)

        times_ms = sol.t * 1e3
        I_total_MA = sol.y[0] * 1e-6
        I_r_MA = sol.y[1] * 1e-6

        return times_ms, I_total_MA, I_r_MA
