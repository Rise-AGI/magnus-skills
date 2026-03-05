"""
Data tables and derived quantities from:
Bernard et al. (MILC Collaboration),
"QCD Thermodynamics with an Improved Lattice Action",
Phys. Rev. D 56, 5584 (1997). hep-lat/9612025.

All figures import from this module to avoid code duplication.
"""

import numpy as np


# ============================================================================
# Table 1: Hadron masses on 8^3 x 16 lattice (improved Wilson fermions)
# Points marked with on_crossover=True lie on the Nt=4 thermal crossover.
# ============================================================================

TABLE1 = [
    # (beta, kappa, u0, n_props, aMPS, aMPS_err, aMV, aMV_err, aMN, aMN_err,
    #  MPS_MV, MPS_MV_err, MN_MV, MN_MV_err, on_crossover)
    (6.40, 0.145,  0.826, 81,  0.931, 0.004, 1.351, 0.018, 2.14,  0.03,
     0.689, 0.010, 1.58, 0.03, False),
    (6.40, 0.1475, 0.828, 30,  0.664, 0.008, 1.26,  0.06,  1.63,  0.11,
     0.527, 0.026, 1.29, 0.11, True),
    (6.60, 0.140,  0.834, 64,  1.173, 0.004, 1.481, 0.009, 2.30,  0.03,
     0.792, 0.006, 1.55, 0.02, False),
    (6.60, 0.143,  0.841, 144, 0.927, 0.004, 1.280, 0.008, 1.958, 0.018,
     0.724, 0.006, 1.530, 0.015, True),
    (6.60, 0.146,  0.855, 40,  0.468, 0.015, 1.04,  0.13,  1.34,  0.05,
     0.45,  0.06,  1.29, 0.17, False),
    (6.80, 0.1325, 0.842, 79,  1.494, 0.003, 1.700, 0.007, 2.651, 0.013,
     0.879, 0.004, 1.56, 0.01, False),
    (6.80, 0.137,  0.849, 120, 1.187, 0.003, 1.421, 0.006, 2.190, 0.010,
     0.835, 0.004, 1.541, 0.010, True),
    (6.80, 0.140,  0.857, 43,  0.885, 0.008, 1.182, 0.016, 1.75,  0.04,
     0.749, 0.012, 1.48, 0.04, False),
    (7.20, 0.118,  0.864, 143, 1.915, 0.003, 1.994, 0.003, 3.110, 0.005,
     0.960, 0.002, 1.560, 0.003, True),
    (7.30, 0.114,  0.8695, 30, 2.043, 0.004, 2.106, 0.005, 3.297, 0.011,
     0.970, 0.003, 1.614, 0.006, True),
]


# ============================================================================
# Table 2: Heavy quark potential fits along Nt=4 crossover (improved Wilson)
# Columns: beta, kappa, n_configs, rmin_rmax, t, aV0, aV0_err,
#           a2sigma, a2sigma_err, b(=e), b_err, f, f_err, r0_over_a, r0_over_a_err
# ============================================================================

TABLE2 = [
    (6.40, 0.1475, 30,  "1.41-4.47", 2, 1.0,  0.3,  0.41,  0.08,  0.7,  0.2,
     5.6, 0.6, 1.52, 0.04),
    (6.60, 0.1430, 108, "1.41-6.93", 2, 0.65, 0.09, 0.42,  0.03,  0.34, 0.08,
     3.21, 0.24, 1.77, 0.02),
    (6.80, 0.1370, 95,  "1.41-6.93", 2, 0.70, 0.06, 0.346, 0.015, 0.38, 0.05,
     2.45, 0.18, 1.913, 0.013),
    (7.20, 0.1180, 117, "1.41-5.66", 2, 0.65, 0.02, 0.253, 0.006, 0.33, 0.02,
     1.06, 0.08, 2.287, 0.013),
]


# ============================================================================
# Table 3: Heavy quark potential fits along crossovers (Kogut-Susskind fermions)
# The last entry (*) is Nt=6 crossover; the first three are Nt=4.
# ============================================================================

TABLE3 = [
    # (beta, amq, n_configs, rmin_rmax, t, aV0, aV0_err,
    #  a2sigma, a2sigma_err, b, b_err, f, f_err, r0_over_a, r0_over_a_err, Nt)
    (5.2875, 0.025, 55,  "1.41-6.93", 2, 0.80, 0.10, 0.30,  0.03,  0.46, 0.10,
     1.46, 0.20, 1.99, 0.04, 4),
    (5.3200, 0.050, 67,  "2.24-6.93", 2, 0.68, 0.22, 0.29,  0.04,  0.2,  0.3,
     5.7, 1.1, 2.17, 0.11, 4),
    (5.3750, 0.100, 90,  "1.00-5.66", 3, 0.62, 0.08, 0.288, 0.023, 0.26, 0.07,
     0.56, 0.12, 2.20, 0.04, 4),
    (5.415,  0.0125, 280, "2.24-6.71", 3, 0.76, 0.02, 0.130, 0.005, 0.36, 0.02,
     1.0, 0.2, 3.14, 0.05, 6),
]


# ============================================================================
# KS meson mass ratios (from Blum et al., PRD 51, 5153 (1995), Table 1)
# Approximate values along Nt=4 crossover for 2-flavor KS
# ============================================================================

KS_MESON_DATA = [
    # (beta, amq, MPS_MV, MPS_MV_err)
    (5.2875, 0.025, 0.325, 0.035),
    (5.3200, 0.050, 0.545, 0.025),
    (5.3750, 0.100, 0.700, 0.015),
    (5.415,  0.0125, 0.35, 0.04),  # Nt=6
]


# ============================================================================
# Derived quantity computations
# ============================================================================

def get_crossover_data():
    """Return only the crossover points from Table 1 (on_crossover=True)."""
    return [row for row in TABLE1 if row[14]]


def compute_tc_over_mv(aMV, aMV_err, Nt=4):
    """
    Compute Tc/MV = 1/(Nt * aMV) with error propagation.
    Tc = 1/(Nt * a), so Tc/MV = 1/(Nt * aMV).
    """
    val = 1.0 / (Nt * aMV)
    err = aMV_err / (Nt * aMV**2)
    return val, err


def compute_tc_over_sqrtsigma(a2sigma, a2sigma_err, Nt=4):
    """
    Compute Tc/sqrt(sigma) = 1/(Nt * sqrt(a^2 sigma)).
    """
    sqrts = np.sqrt(a2sigma)
    val = 1.0 / (Nt * sqrts)
    err = a2sigma_err / (2 * Nt * a2sigma * sqrts)
    return val, err


def compute_r0_tc(r0_over_a, r0_over_a_err, Nt=4):
    """
    Compute r0*Tc = (r0/a) / Nt.
    """
    val = r0_over_a / Nt
    err = r0_over_a_err / Nt
    return val, err


def compute_a_sqrtsigma(a2sigma, a2sigma_err):
    """
    Compute a*sqrt(sigma) = sqrt(a^2 sigma).
    """
    val = np.sqrt(a2sigma)
    err = a2sigma_err / (2 * val)
    return val, err


def compute_r0_sqrtsigma(r0_over_a, r0_over_a_err, a2sigma, a2sigma_err):
    """
    Compute r0*sqrt(sigma) = (r0/a) * sqrt(a^2 sigma).
    """
    sqrts = np.sqrt(a2sigma)
    val = r0_over_a * sqrts
    err = val * np.sqrt((r0_over_a_err / r0_over_a)**2 +
                         (a2sigma_err / (2 * a2sigma))**2)
    return val, err


def compute_a_over_r0(r0_over_a, r0_over_a_err):
    """
    Compute a/r0 = 1/(r0/a).
    """
    val = 1.0 / r0_over_a
    err = r0_over_a_err / r0_over_a**2
    return val, err


def compute_mv_r0(aMV, aMV_err, r0_over_a, r0_over_a_err):
    """
    Compute MV*r0 = aMV * (r0/a).
    """
    val = aMV * r0_over_a
    err = val * np.sqrt((aMV_err / aMV)**2 + (r0_over_a_err / r0_over_a)**2)
    return val, err


def strong_coupling_prediction(beta):
    """
    Strong coupling prediction for a^2*K (string tension).
    a^2 K ~ -ln(beta/4) for SU(2).
    For SU(3): a^2 K ~ -ln(beta/6) approximately.
    """
    return -np.log(beta / 6.0)


def asymptotic_freedom_prediction(beta):
    """
    Weak coupling (asymptotic freedom) prediction for SU(3).
    a^2 K ~ C * exp(-6*pi^2/(33-2*Nf) * beta)
    For Nf=2: a^2 K ~ C * exp(-6*pi^2/29 * beta)
    The overall constant is not fixed by perturbation theory.
    """
    Nf = 2
    b0 = (33 - 2 * Nf) / (48 * np.pi**2)
    return np.exp(-1.0 / (2 * b0 * 6.0 / beta))


# ============================================================================
# Phase diagram data: kappa_T(beta) from Polyakov loop crossover
# and kappa_c(beta) from vanishing pion mass
# ============================================================================

# Crossover kappa_T(beta) values (from * entries in Table 1)
KAPPA_T = {
    6.40: 0.1475,
    6.60: 0.143,
    6.80: 0.137,
    7.20: 0.118,
    7.30: 0.114,
}

# Critical kappa_c(beta) estimates where MPS -> 0
# Extrapolated from spectroscopy data (linear in 1/kappa vs (aMPS)^2)
# These are approximate estimates from the paper's Figure 9
KAPPA_C = {
    6.40: 0.1495,
    6.60: 0.1480,
    6.80: 0.1430,
    7.20: 0.1195,
    7.30: 0.1155,
}


def estimate_kappa_c(beta_val, kappa_vals, amps_vals):
    """
    Estimate kappa_c by linear extrapolation of (aMPS)^2 vs 1/kappa to zero.
    """
    inv_kappa = 1.0 / np.array(kappa_vals)
    amps2 = np.array(amps_vals)**2
    # Linear fit: amps2 = a * inv_kappa + b
    coeffs = np.polyfit(inv_kappa, amps2, 1)
    # At amps2 = 0: inv_kappa_c = -b/a
    inv_kappa_c = -coeffs[1] / coeffs[0]
    return 1.0 / inv_kappa_c


# Quenched SU(3) Nt=4 Tc/sqrt(sigma) from Boyd et al. (Ref [36])
QUENCHED_TC_SQRTSIGMA_NT4 = 0.630
