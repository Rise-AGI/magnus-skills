"""
Perturbative QCD at NNLO for quenched lattice QCD.

Based on: Giusti, Rapuano, Talevi, Vladikas,
"The QCD Chiral Condensate from the Lattice" (1998)
hep-lat/9807014, Nucl. Phys. B538 (1999) 249-277.

Implements:
- Beta-function coefficients (Eq. B.4)
- Quark mass anomalous dimension (Eq. B.5)
- Running coupling alpha_s at NNLO (Eq. B.7)
- RI/MSbar matching coefficient (Eq. B.1-B.3)
- RG evolution coefficient c_S^MSbar (Eq. B.6)
- Chiral condensate determination from lattice data (Tables 1-4)
"""

import numpy as np


# ============================================================
# Physical constants
# ============================================================
ZETA_3 = 1.2020569031595942  # Riemann zeta(3)

# QCD parameters for quenched approximation
N_C = 3    # Number of colors
N_F = 0    # Number of flavors (quenched)
LAMBDA_QCD_MSBAR = 0.251  # GeV, quenched, from ref [37]
LAMBDA_QCD_MSBAR_ERR = 0.021  # GeV

# Experimental inputs
F_CHI = 0.1282      # GeV, pion decay constant in chiral limit
M_KSTAR = 0.8931    # GeV, K* meson mass for scale setting


# ============================================================
# Beta-function coefficients (Eq. B.4)
# ============================================================
def beta_coefficients(n_c=N_C, n_f=N_F):
    """Compute beta-function coefficients beta_0, beta_1, beta_2^MSbar."""
    beta_0 = (11.0 / 3.0) * n_c - (2.0 / 3.0) * n_f

    beta_1 = (34.0 / 3.0) * n_c**2 - (10.0 / 3.0) * n_c * n_f \
             - ((n_c**2 - 1.0) / n_c) * n_f

    beta_2 = (2857.0 / 54.0) * n_c**3 \
             + ((n_c**2 - 1.0)**2 / (4.0 * n_c**2)) * n_f \
             - (205.0 / 36.0) * (n_c**2 - 1.0) * n_f \
             - (1415.0 / 54.0) * n_c**2 * n_f \
             + (11.0 / 18.0) * ((n_c**2 - 1.0) / n_c) * n_f**2 \
             + (79.0 / 54.0) * n_c * n_f**2

    return beta_0, beta_1, beta_2


# ============================================================
# Anomalous dimension coefficients (Eq. B.5)
# ============================================================
def gamma_m_coefficients(n_c=N_C, n_f=N_F):
    """Compute quark mass anomalous dimension coefficients gamma_m^(i)."""
    cf = (n_c**2 - 1.0) / n_c  # Casimir C_F

    gamma_0 = 3.0 * cf

    gamma_1 = cf / n_c * (-3.0 / 4.0 + (203.0 / 12.0) * n_c**2
                          - (5.0 / 3.0) * n_c * n_f)

    gamma_2 = cf / n_c**2 * (
        129.0 / 8.0 - (129.0 / 8.0) * n_c**2
        + (11413.0 / 108.0) * n_c**4
        + n_f * (23.0 / 2.0 * n_c - (1177.0 / 54.0) * n_c**3
                 - 12.0 * n_c * ZETA_3 - 12.0 * n_c**3 * ZETA_3)
        - (35.0 / 27.0) * n_c**2 * n_f**2
    )

    return gamma_0, gamma_1, gamma_2


# ============================================================
# Running coupling alpha_s at NNLO (Eq. B.7)
# ============================================================
def alpha_s_nnlo(mu, lambda_qcd=LAMBDA_QCD_MSBAR, n_c=N_C, n_f=N_F):
    """
    Compute alpha_s^MSbar(mu) at NNLO.

    Parameters
    ----------
    mu : float or array
        Renormalization scale in GeV.
    lambda_qcd : float
        Lambda_QCD^MSbar in GeV.

    Returns
    -------
    alpha_s / (4*pi)
    """
    beta_0, beta_1, beta_2 = beta_coefficients(n_c, n_f)

    q2 = (mu / lambda_qcd)**2
    L = np.log(q2)
    LL = np.log(L)

    # Leading order
    result = 1.0 / (beta_0 * L)

    # NLO correction
    result -= beta_1 / beta_0**3 * LL / L**2

    # NNLO correction
    result += 1.0 / (beta_0**5 * L**3) * (
        beta_1**2 * LL**2
        - beta_1**2 * LL
        + beta_2 * beta_0
        - beta_1**2
    )

    return result


def alpha_s_value(mu, lambda_qcd=LAMBDA_QCD_MSBAR, n_c=N_C, n_f=N_F):
    """Return alpha_s (not alpha_s/(4*pi))."""
    return 4.0 * np.pi * alpha_s_nnlo(mu, lambda_qcd, n_c, n_f)


# ============================================================
# RI/MSbar matching coefficient (Eq. B.1-B.3)
# ============================================================
def ri_msbar_matching_coefficients(n_c=N_C, n_f=N_F):
    """Compute C^(1) and C^(2) for RI/MSbar matching."""
    c1 = 8.0 * (n_c**2 - 1.0) / (4.0 * n_c)

    c2 = (n_c**2 - 1.0) / (96.0 * n_c**2) * (
        -309.0 + 3029.0 * n_c**2
        - 288.0 * ZETA_3 - 576.0 * n_c**2 * ZETA_3
        - 356.0 * n_c * n_f
    )

    return c1, c2


def delta_z_ri_msbar(mu, lambda_qcd=LAMBDA_QCD_MSBAR, n_c=N_C, n_f=N_F):
    """
    Compute RI/MSbar matching coefficient Delta_Z^{RI/MSbar}.
    Eq. B.1.
    """
    c1, c2 = ri_msbar_matching_coefficients(n_c, n_f)
    a_s_over_4pi = alpha_s_nnlo(mu, lambda_qcd, n_c, n_f)
    a_s = 4.0 * np.pi * a_s_over_4pi

    dz = 1.0 + (a_s / (4.0 * np.pi)) * c1 \
         + (a_s / (4.0 * np.pi))**2 * c2

    return dz


# ============================================================
# Wave function renormalization correction (Eq. A.6-A.7)
# ============================================================
def delta_q_correction(n_c=N_C, n_f=N_F):
    """Compute Delta_q^(2) for Z_q/Z_q' matching. Eq. A.7."""
    return (n_c**2 - 1.0) / (16.0 * n_c**2) * (3.0 + 22.0 * n_c**2 - 4.0 * n_c * n_f)


# ============================================================
# RG evolution coefficient c_S^MSbar (Eq. B.6)
# ============================================================
def evolution_coefficient(mu, lambda_qcd=LAMBDA_QCD_MSBAR, n_c=N_C, n_f=N_F):
    """
    Compute the RG evolution coefficient c_S^MSbar(mu) at NNLO.
    Eq. B.6.
    """
    beta_0, beta_1, beta_2 = beta_coefficients(n_c, n_f)
    gamma_0, gamma_1, gamma_2 = gamma_m_coefficients(n_c, n_f)

    # Barred quantities
    beta_1_bar = beta_1 / beta_0
    beta_2_bar = beta_2 / beta_0

    # gamma_S = -gamma_m (scalar condensate anomalous dimension)
    gamma_s0_bar = -gamma_0 / (2.0 * beta_0)
    gamma_s1_bar = -gamma_1 / (2.0 * beta_0)
    gamma_s2_bar = -gamma_2 / (2.0 * beta_0)

    a_s = alpha_s_value(mu, lambda_qcd, n_c, n_f)
    a_s_4pi = a_s / (4.0 * np.pi)

    # Leading power
    c_s = a_s**gamma_s0_bar

    # NLO correction
    nlo_term = 1.0 + a_s_4pi * (gamma_s1_bar - beta_1_bar * gamma_s0_bar)

    # NNLO correction
    bracket = gamma_s1_bar - beta_1_bar * gamma_s0_bar
    nnlo_term = 0.5 * a_s_4pi**2 * (
        bracket**2
        + gamma_s2_bar
        + beta_1_bar**2 * gamma_s0_bar
        - beta_1_bar * gamma_s1_bar
        - beta_2_bar * gamma_s0_bar
    )

    c_s *= (nlo_term + nnlo_term)

    return c_s


def run_condensate(condensate_mu, mu_from, mu_to,
                   lambda_qcd=LAMBDA_QCD_MSBAR, n_c=N_C, n_f=N_F):
    """
    Run the chiral condensate from scale mu_from to mu_to using NNLO RG.
    Eq. 13.
    """
    c_from = evolution_coefficient(mu_from, lambda_qcd, n_c, n_f)
    c_to = evolution_coefficient(mu_to, lambda_qcd, n_c, n_f)
    return condensate_mu * c_to / c_from


# ============================================================
# Lattice data from Tables 1-4
# ============================================================

# Table 1: Lattice parameters
LATTICE_RUNS = {
    'C60a': {'beta': 6.0, 'action': 'Clover', 'volume': '18^3x64',
             'n_configs': 490, 'a_inv': 2.12, 'a_inv_err': 0.06,
             'kappa': [0.1425, 0.1432, 0.1440]},
    'C60b': {'beta': 6.0, 'action': 'Clover', 'volume': '24^3x40',
             'n_configs': 600, 'a_inv': 2.16, 'a_inv_err': 0.04,
             'kappa': [0.1425, 0.1432, 0.1440]},
    'C62':  {'beta': 6.2, 'action': 'Clover', 'volume': '24^3x64',
             'n_configs': 250, 'a_inv': 2.7, 'a_inv_err': 0.1,
             'kappa': [0.14144, 0.14184, 0.14224, 0.14264]},
    'C64':  {'beta': 6.4, 'action': 'Clover', 'volume': '24^3x64',
             'n_configs': 400, 'a_inv': 4.0, 'a_inv_err': 0.2,
             'kappa': [0.1400, 0.1403, 0.1406, 0.1409]},
    'W60':  {'beta': 6.0, 'action': 'Wilson', 'volume': '18^3x64',
             'n_configs': 320, 'a_inv': 2.26, 'a_inv_err': 0.05,
             'kappa': [0.1530, 0.1540, 0.1550]},
    'W62':  {'beta': 6.2, 'action': 'Wilson', 'volume': '24^3x64',
             'n_configs': 250, 'a_inv': 3.00, 'a_inv_err': 0.09,
             'kappa': [0.1510, 0.1515, 0.1520, 0.1526]},
    'W64':  {'beta': 6.4, 'action': 'Wilson', 'volume': '24^3x64',
             'n_configs': 400, 'a_inv': 4.1, 'a_inv_err': 0.2,
             'kappa': [0.1488, 0.1492, 0.1496, 0.1500]},
}

# Table 2: Renormalization constants (NP method, RI scheme, chiral limit)
RENORM_CONSTANTS = {
    'C60': {'Z_S': (0.83, 0.02), 'Z_P': (0.41, 0.06),
            'Z_P_over_Z_S': (0.49, 0.06), 'Z_A': (1.05, 0.03)},
    'W60': {'Z_S': (0.68, 0.01), 'Z_P': (0.45, 0.06),
            'Z_P_over_Z_S': (0.66, 0.08), 'Z_A': (0.81, 0.01)},
    'C62': {'Z_S': (0.85, 0.02), 'Z_P': (0.47, 0.05),
            'Z_P_over_Z_S': (0.56, 0.05), 'Z_A': (1.02, 0.02)},
    'W62': {'Z_S': (0.72, 0.01), 'Z_P': (0.50, 0.05),
            'Z_P_over_Z_S': (0.69, 0.07), 'Z_A': (0.81, 0.01)},
    'C64': {'Z_S': (0.85, 0.02), 'Z_P': (0.55, 0.03),
            'Z_P_over_Z_S': (0.65, 0.02), 'Z_A': (1.01, 0.01)},
    'W64': {'Z_S': (0.74, 0.01), 'Z_P': (0.57, 0.04),
            'Z_P_over_Z_S': (0.77, 0.05), 'Z_A': (0.82, 0.01)},
}

# Table 3: Subtracted chiral condensate in lattice units
# Values are (value, statistical_error); overall minus sign omitted
CONDENSATE_LATTICE = {
    'C60a': {'chi1_over_ZA': (0.0040, 0.0003), 'chi2_over_ZA': (0.0040, 0.0003)},
    'C60b': {'chi1_over_ZA': (0.0039, 0.0002), 'chi2_over_ZA': (0.0039, 0.0002)},
    'W60':  {'chi1_over_ZA': (0.0045, 0.0003), 'chi2_over_ZA': (0.0051, 0.0004)},
    'C62':  {'chi1_over_ZA': (0.0018, 0.0002), 'chi2_over_ZA': (0.0015, 0.0002)},
    'W62':  {'chi1_over_ZA': (0.0019, 0.0001), 'chi2_over_ZA': (0.0016, 0.0001)},
    'C64':  {'chi1_over_ZA': (0.0007, 0.0001), 'chi2_over_ZA': (0.0007, 0.0001)},
    'W64':  {'chi1_over_ZA': (0.0008, 0.0001), 'chi2_over_ZA': (0.0007, 0.0001)},
}

# Table 4: C^HS, C^AWI slopes and physical condensate (GeV^3)
CONDENSATE_PHYSICAL = {
    'C60a': {'C_HS': (2.98, 0.08), 'C_AWI': (3.9, 0.1),
             'psi1': (None, None, None), 'psi2': (0.0141, 0.0005, 0.0021)},
    'C60b': {'C_HS': (3.04, 0.07), 'C_AWI': (4.1, 0.1),
             'psi1': (None, None, None), 'psi2': (0.0146, 0.0003, 0.0021)},
    'W60':  {'C_HS': (2.40, 0.05), 'C_AWI': (3.01, 0.07),
             'psi1': (0.0150, 0.0003, 0.0002), 'psi2': (0.0152, 0.0004, 0.0020)},
    'C62':  {'C_HS': (2.9, 0.1), 'C_AWI': (3.7, 0.2),
             'psi1': (None, None, None), 'psi2': (0.0147, 0.0008, 0.0016)},
    'W62':  {'C_HS': (2.52, 0.08), 'C_AWI': (2.98, 0.09),
             'psi1': (0.0156, 0.0005, 0.0002), 'psi2': (0.0158, 0.0005, 0.0016)},
    'C64':  {'C_HS': (3.5, 0.2), 'C_AWI': (4.2, 0.2),
             'psi1': (None, None, None), 'psi2': (0.0187, 0.0009, 0.0010)},
    'W64':  {'C_HS': (2.9, 0.1), 'C_AWI': (3.2, 0.1),
             'psi1': (0.0172, 0.0007, 0.0002), 'psi2': (0.0179, 0.0008, 0.0013)},
}


def compute_condensate_method1(run_key, mu_target=2.0):
    """
    Compute chiral condensate using method 1 (Eq. 41):
    <psibar psi>_1 = -0.5 * a^{-1} * f_chi^2 * Z_S * C^HS

    Only available for Wilson action.
    """
    run = LATTICE_RUNS[run_key]
    beta_key = run['action'][0] + str(int(run['beta'] * 10))
    rc = RENORM_CONSTANTS[beta_key]

    a_inv = run['a_inv']  # GeV
    z_s = rc['Z_S'][0]
    c_hs = CONDENSATE_PHYSICAL[run_key]['C_HS'][0]

    # Bare condensate in RI scheme at scale mu ~ a^{-1}
    condensate_ri = 0.5 * a_inv * F_CHI**2 * z_s * c_hs

    # Convert RI -> MSbar
    dz = delta_z_ri_msbar(a_inv)
    condensate_msbar = dz * condensate_ri

    # Run to target scale
    condensate_final = run_condensate(condensate_msbar, a_inv, mu_target)

    return condensate_final


def compute_condensate_method2(run_key, mu_target=2.0):
    """
    Compute chiral condensate using method 2 (Eq. 42):
    <psibar psi>_2 = -0.5 * a^{-1} * f_chi^2 * (Z_P/Z_A) * C^AWI
    """
    run = LATTICE_RUNS[run_key]
    beta_str = str(int(run['beta'] * 10))
    beta_key = run['action'][0] + beta_str
    rc = RENORM_CONSTANTS[beta_key]

    a_inv = run['a_inv']  # GeV
    z_p = rc['Z_P'][0]
    z_a = rc['Z_A'][0]
    c_awi = CONDENSATE_PHYSICAL[run_key]['C_AWI'][0]

    # Bare condensate in RI scheme at scale mu ~ a^{-1}
    condensate_ri = 0.5 * a_inv * F_CHI**2 * (z_p / z_a) * c_awi

    # Convert RI -> MSbar
    dz = delta_z_ri_msbar(a_inv)
    condensate_msbar = dz * condensate_ri

    # Run to target scale
    condensate_final = run_condensate(condensate_msbar, a_inv, mu_target)

    return condensate_final
