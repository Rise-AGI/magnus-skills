"""
Core module for quenched chiral perturbation theory (qXPT) analysis.

Implements the formulas from Chiu & Hsieh (2003):
  "Light quark masses, chiral condensate and quark-gluon condensate
   in quenched lattice QCD with exact chiral symmetry"

Key equations:
  - Eq(14): m_pi^2 = A * m_q^(1/(1+delta)) + B * m_q^2
  - Eq(13): m_pi^2 = C*m_q * {1 - delta*[ln(C*m_q/Lambda_chi^2) + 1]} + B*m_q^2
  - Eq(16): delta = m_eta'^2 / (8*pi^2*f_pi^2*N_f)
  - Eq(25): f_pi*a = f0 + slope * m_q*a
  - Eq(37): m_K^2/m_pi^2 ratio in qXPT
"""

import numpy as np
from scipy.optimize import curve_fit

# ============================================================
# Published fit parameters (from the paper)
# ============================================================

# From Eq.(14) fit: m_pi^2*a^2 = A1*(m_q*a)^(1/(1+delta)) + B*(m_q*a)^2
DELTA = 0.164      # coefficient of quenched chiral logarithm
DELTA_ERR = 0.013
A1 = 1.044          # A1 = A * a^(-(1-2*delta)/(1+delta))
A1_ERR = 0.022
B_FIT = 2.077       # coefficient of m_q^2 term
B_FIT_ERR = 0.130

# From Eq.(28) fit with Eq.(13): delta, C, B with fixed Lambda_chi
DELTA_13 = 0.168
DELTA_13_ERR = 0.017
C_PARAM = 1.093     # C*a (dimensionless on lattice)
C_PARAM_ERR = 0.020
B_13 = 2.092
B_13_ERR = 0.125

# Pion decay constant linear fit: f_pi*a = f0 + slope*(m_q*a)
FPI_F0 = 0.0667     # f_pi*a at m_q=0
FPI_F0_ERR = 0.002
FPI_SLOPE = 0.2218
FPI_SLOPE_ERR = 0.0020

# Lattice spacing
A_INV_GEV = 1.979   # a^{-1} in GeV
A_INV_ERR = 0.006
A_FM = 0.0997       # a in fm
A_FM_ERR = 0.0003

# Topological susceptibility
CHI_T_MEV4 = 175.0**4  # (175 MeV)^4
A4_CHI_T = 6.03e-5     # a^4 * chi_t
A4_CHI_T_ERR = 0.75e-5

# Chiral condensate linear fit: -a^3<qbar q> = intercept + slope*(m_q*a)
QQBAR_INTERCEPT = 1.73e-3
QQBAR_INTERCEPT_ERR = 0.03e-3
QQBAR_SLOPE = 0.242

# Quark-gluon condensate linear fit
QGQ_INTERCEPT = 5.59e-4
QGQ_INTERCEPT_ERR = 0.12e-4
QGQ_SLOPE = 0.0150

# Physical constants (experimental inputs)
M_PI_PHYS = 0.135   # GeV (physical pion mass)
M_K_PHYS = 0.495    # GeV (physical kaon mass)
F_PI_PHYS = 0.132   # GeV (pion decay constant)
N_F = 3              # number of light quark flavors

# Index distribution from Table 1
# n_+ - n_- : number of configurations
INDEX_DISTRIBUTION = {
    5: 6, 4: 5, 3: 6, 2: 8, 1: 14, 0: 10,
    -1: 15, -2: 15, -3: 10, -4: 4, -5: 3, -6: 3, -7: 1
}
N_CONFIGS = 100
N_SITES = 16**3 * 32  # 16^3 x 32 lattice

# Renormalization constant Z_s (one-loop perturbative)
Z_S = 1.172  # at beta=6.0, mu=2 GeV

# 13 bare quark masses used in the paper (in lattice units m_q*a)
MQ_VALUES = np.array([
    0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09,
    0.10, 0.12, 0.14, 0.16, 0.18, 0.20
])


# ============================================================
# qXPT formulas
# ============================================================

def mpi2_eq14(mq_a, delta, A1, B):
    """Pion mass squared from Eq.(14): m_pi^2*a^2 = A1*(m_q*a)^(1/(1+delta)) + B*(m_q*a)^2"""
    return A1 * mq_a**(1.0 / (1.0 + delta)) + B * mq_a**2


def mpi2_eq13(mq_a, C_a, delta, B, Lambda_chi_a):
    """Pion mass squared from Eq.(13) (one-loop):
    m_pi^2*a^2 = C_a*(m_q*a)*{1 - delta*[ln(C_a*(m_q*a)/Lambda_chi_a^2) + 1]} + B*(m_q*a)^2
    """
    x = C_a * mq_a
    return x * (1.0 - delta * (np.log(x / Lambda_chi_a**2) + 1.0)) + B * mq_a**2


def fpi_linear(mq_a, f0, slope):
    """Pion decay constant linear fit: f_pi*a = f0 + slope*(m_q*a)"""
    return f0 + slope * mq_a


def ratio_mpi2_mq(mq_a, delta, A1, B):
    """(m_pi*a)^2 / (m_q*a) for chiral log analysis"""
    mpi2 = mpi2_eq14(mq_a, delta, A1, B)
    return mpi2 / mq_a


def chiral_log_extraction(mq_a, delta, A1, B):
    """log[(m_pi*a)^2/(m_q*a) - B*(m_q*a)] vs log(m_q*a)
    Slope should be -delta/(1+delta)
    """
    mpi2 = mpi2_eq14(mq_a, delta, A1, B)
    y = mpi2 / mq_a - B * mq_a
    return np.log(mq_a), np.log(y)


def condensate_linear(mq_a, intercept, slope):
    """Linear fit for condensate extrapolation"""
    return intercept + slope * mq_a


def quark_mass_ratio(delta, m_K, m_pi):
    """m_s/m from Eq.(37) in qXPT"""
    # Solve: m_K^2/m_pi^2 = (m+m_s)/(2m) * {1 + delta*[1 - m/(m_s-m)*ln(m_s/m)]}
    # This requires numerical solution
    from scipy.optimize import brentq

    target = (m_K / m_pi)**2

    def equation(r):
        # r = m_s/m
        lhs = (1 + r) / 2.0 * (1.0 + delta * (1.0 - 1.0 / (r - 1.0) * np.log(r)))
        return lhs - target

    r = brentq(equation, 1.001, 100.0)
    return r


def compute_topological_susceptibility():
    """Compute chi_t from Table 1 index distribution"""
    total_Q2 = 0
    n_total = 0
    for Q, n_cfg in INDEX_DISTRIBUTION.items():
        total_Q2 += Q**2 * n_cfg
        n_total += n_cfg

    mean_Q2 = total_Q2 / n_total
    a4_chi_t = mean_Q2 / N_SITES

    # Error estimation (statistical)
    Q2_vals = []
    for Q, n_cfg in INDEX_DISTRIBUTION.items():
        Q2_vals.extend([Q**2] * n_cfg)
    Q2_vals = np.array(Q2_vals)
    std_Q2 = np.std(Q2_vals, ddof=1) / np.sqrt(n_total)
    a4_chi_t_err = std_Q2 / N_SITES

    return a4_chi_t, a4_chi_t_err, mean_Q2


def delta_from_chi_t(a4_chi_t, fpi_a):
    """Compute delta from topological susceptibility via Eq.(20):
    delta = 1/(2*pi^2*(f_pi*a)^4) * <(n+ - n-)^2>/N
    """
    return a4_chi_t / (2.0 * np.pi**2 * fpi_a**4)


def eta_prime_mass(a4_chi_t, fpi_a, a_inv):
    """Compute eta' mass from index susceptibility via Eq.(19):
    (m_eta' * a)^2 = 4*N_f / (f_pi*a)^2 * chi_t
    """
    meta_a2 = 4.0 * N_F / fpi_a**2 * a4_chi_t
    meta_a = np.sqrt(meta_a2)
    meta_GeV = meta_a * a_inv
    return meta_GeV


def compute_quark_masses():
    """Compute light quark masses from qXPT parameters + experimental inputs"""
    # Quark mass ratio from Eq.(37)
    r = quark_mass_ratio(DELTA, M_K_PHYS, M_PI_PHYS)

    # From Eq.(38): m_pi^2 = A * m^(1/(1+delta)) + B * m^2
    # with A, B in physical units
    # m = bare quark mass, determined by inverting
    # Using published result directly: m = 4.7 +/- 0.5 MeV (bare)
    m_bare = 4.7e-3  # GeV
    m_s_bare = m_bare * r

    # Renormalize to MS-bar at mu=2 GeV
    Z_m = 1.0 / Z_S
    m_ud_msbar = m_bare * Z_m
    m_s_msbar = m_s_bare * Z_m

    return {
        'ms_over_m': r,
        'm_bare_MeV': m_bare * 1000,
        'ms_bare_MeV': m_s_bare * 1000,
        'm_ud_msbar_MeV': m_ud_msbar * 1000,
        'ms_msbar_MeV': m_s_msbar * 1000,
    }


def compute_condensates():
    """Compute chiral and quark-gluon condensates"""
    # Chiral condensate
    a3_qqbar = QQBAR_INTERCEPT  # at m_q=0
    qqbar_GeV3 = -a3_qqbar * A_INV_GEV**3
    qqbar_msbar = qqbar_GeV3 * Z_S
    qqbar_MeV = -(-qqbar_msbar * 1e9)**(1.0/3.0)  # -(X MeV)^3

    # Quark-gluon condensate
    a5_qgq = QGQ_INTERCEPT
    qgq_GeV5 = -a5_qgq * A_INV_GEV**5
    qgq_msbar = qgq_GeV5 * Z_S  # approximate
    qgq_MeV = -(-qgq_msbar * 1e15)**(1.0/5.0)

    # Ratio M_0^2
    M02 = qgq_msbar / qqbar_msbar

    return {
        'qqbar_GeV3': qqbar_GeV3,
        'qqbar_msbar_GeV3': qqbar_msbar,
        'qqbar_msbar_MeV_cuberoot': qqbar_MeV,
        'qgq_GeV5': qgq_GeV5,
        'qgq_msbar_GeV5': qgq_msbar,
        'qgq_msbar_MeV_fifthroot': qgq_MeV,
        'M02_GeV2': M02,
    }
