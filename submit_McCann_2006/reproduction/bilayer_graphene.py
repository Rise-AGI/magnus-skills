"""
Core module for reproducing McCann, Phys. Rev. B 74, 161403(R) (2006).
Asymmetry gap in the electronic band structure of bilayer graphene.

Implements:
- 4x4 tight-binding Hamiltonian (Eq. 3)
- Band energies (Eq. 4)
- Self-consistent Hartree calculation of layer asymmetry Delta(n)
- Analytic approximations (Eqs. 1, 6)
- Cyclotron mass (Eq. 7)
- Landau level spectrum (Eqs. 8, 9)
"""

import numpy as np
from scipy.constants import e as e_charge, epsilon_0, hbar, m_e

# ---------------------------------------------------------------------------
# Physical parameters from Ref [11] of the paper
# ---------------------------------------------------------------------------
GAMMA1_EV = 0.39        # interlayer coupling [eV]
GAMMA1 = GAMMA1_EV * e_charge   # [J]
C0 = 3.35e-10           # interlayer separation [m]
VF = 8.0e5              # in-plane velocity [m/s]

# Electron mass for cyclotron mass normalization
M_E = m_e

# Effective mass at Delta=0: m = gamma1/(2v^2)
M_EFF = GAMMA1 / (2 * VF**2)


def screening_lambda(eps_r=1.0):
    """Screening parameter Lambda = e^2 gamma1 c0 / (2pi hbar^2 v^2 eps_r eps_0)."""
    return (e_charge**2 * GAMMA1 * C0) / (
        2 * np.pi * hbar**2 * VF**2 * eps_r * epsilon_0)


def cap_factor(eps_r=1.0):
    """Capacitance prefactor: e^2 c0 / (eps_r eps_0) [J*m^2]."""
    return e_charge**2 * C0 / (eps_r * epsilon_0)


# ---------------------------------------------------------------------------
# Band structure
# ---------------------------------------------------------------------------

def band_energies(p_arr, delta):
    """Compute 4 band energies analytically from Eq. (4).

    Parameters
    ----------
    p_arr : array_like
        Momentum magnitudes [kg*m/s]
    delta : float
        Layer asymmetry [J]

    Returns
    -------
    energies : ndarray, shape (N, 4)
        Sorted band energies [J]: [-E2, -E1, +E1, +E2]
    """
    p_arr = np.atleast_1d(p_arr)
    vp2 = (VF * p_arr)**2
    g1 = GAMMA1

    A = g1**2 / 2 + delta**2 / 4 + vp2
    B = np.sqrt(g1**4 / 4 + vp2 * (g1**2 + delta**2))

    e1 = np.sqrt(np.maximum(A - B, 0))
    e2 = np.sqrt(np.maximum(A + B, 0))

    return np.column_stack([-e2, -e1, e1, e2])


def diag_batch(p_arr, delta):
    """Diagonalize 4x4 Hamiltonian for array of momenta.

    Basis: (psi_A1, psi_B2, psi_A2, psi_B1) for valley K (xi=+1).

    Returns eigenvalues (N,4) sorted ascending, eigenvectors (N,4,4).
    """
    p_arr = np.atleast_1d(p_arr).astype(float)
    N = len(p_arr)
    H = np.zeros((N, 4, 4))
    vp = VF * p_arr

    H[:, 0, 0] = -delta / 2
    H[:, 1, 1] = delta / 2
    H[:, 2, 2] = delta / 2
    H[:, 3, 3] = -delta / 2
    H[:, 0, 3] = vp
    H[:, 3, 0] = vp
    H[:, 1, 2] = vp
    H[:, 2, 1] = vp
    H[:, 2, 3] = GAMMA1
    H[:, 3, 2] = GAMMA1

    return np.linalg.eigh(H)


def layer2_weights(evecs):
    """Layer-2 weight (B2 idx=1, A2 idx=2) for each eigenstate.

    Returns shape (N, 4) array.
    """
    return np.abs(evecs[:, 1, :])**2 + np.abs(evecs[:, 2, :])**2


def layer1_weights(evecs):
    """Layer-1 weight (A1 idx=0, B1 idx=3) for each eigenstate."""
    return np.abs(evecs[:, 0, :])**2 + np.abs(evecs[:, 3, :])**2


# ---------------------------------------------------------------------------
# Density computation
# ---------------------------------------------------------------------------

def pF_from_density(n):
    """Fermi momentum from 2D density with 4-fold degeneracy.

    n = p_F^2 / (pi hbar^2)
    """
    return np.sqrt(np.pi * hbar**2 * abs(n))


def compute_n2(n_density, delta, eps_r=1.0, Np=600):
    """Compute excess density on layer 2.

    Parameters
    ----------
    n_density : float
        Total excess electron density [m^-2]. Positive = electron doping.
    delta : float
        Layer asymmetry gap [J]
    eps_r : float
        Dielectric constant of bilayer
    Np : int
        Momentum grid points

    Returns
    -------
    n2 : float
        Excess density on layer 2 [m^-2]
    """
    prefactor = 2.0 / (np.pi * hbar**2)

    # --- Valence band redistribution ---
    p_max = 15 * GAMMA1 / VF
    p_val = np.linspace(0, p_max, Np + 1)[1:]  # skip p=0

    _, evecs_d = diag_batch(p_val, delta)
    _, evecs_0 = diag_batch(p_val, 0.0)

    w2_d = layer2_weights(evecs_d)
    w2_0 = layer2_weights(evecs_0)

    # Change in layer-2 weight for valence bands (indices 0, 1)
    dw2 = (w2_d[:, 0] + w2_d[:, 1]) - (w2_0[:, 0] + w2_0[:, 1])
    dn2_val = np.trapezoid(p_val * dw2, p_val) * prefactor

    # --- Partial band filling/emptying ---
    n2_extra = 0.0
    if abs(n_density) > 1e8:
        p_F = pF_from_density(n_density)
        Nc = max(Np // 3, 80)
        p_cond = np.linspace(0, p_F, Nc + 1)[1:]

        _, evecs_c = diag_batch(p_cond, delta)
        w2_c = layer2_weights(evecs_c)

        band_idx = 2 if n_density > 0 else 1
        n2_extra = np.sign(n_density) * np.trapezoid(
            p_cond * w2_c[:, band_idx], p_cond) * prefactor

    return dn2_val + n2_extra


def compute_layer_densities(n_density, delta, eps_r=1.0, Np=600):
    """Compute excess densities on both layers.

    Returns (n1, n2) in m^-2.
    """
    n2 = compute_n2(n_density, delta, eps_r, Np)
    n1 = n_density - n2
    return n1, n2


# ---------------------------------------------------------------------------
# Self-consistent solver
# ---------------------------------------------------------------------------

def self_consistent_delta(n_density, delta0=0.0, eps_r=1.0,
                          max_iter=200, tol=1e-5, mix=0.3):
    """Solve for self-consistent asymmetry gap Delta(n).

    Parameters
    ----------
    n_density : float
        Total excess density [m^-2]
    delta0 : float
        Bare asymmetry [J]
    eps_r : float
        Dielectric constant

    Returns
    -------
    delta : float
        Self-consistent gap [J]
    """
    cf = cap_factor(eps_r)

    # Initial guess
    delta = delta0 + cf * n_density / 2
    if abs(delta) < 1e-26:
        delta = delta0 if abs(delta0) > 1e-26 else 1e-26

    for _ in range(max_iter):
        n2 = compute_n2(n_density, delta, eps_r)
        delta_new = delta0 + cf * n2

        err = abs(delta_new - delta)
        if err < tol * (abs(delta_new) + 1e-26):
            return delta_new

        delta = mix * delta_new + (1 - mix) * delta

    return delta


def self_consistent_delta_at_zero(delta0, eps_r=1.0,
                                  max_iter=200, tol=1e-5, mix=0.3):
    """Self-consistent Delta(0) at zero excess density.

    Simpler version: only valence band redistribution matters.
    """
    return self_consistent_delta(0.0, delta0, eps_r, max_iter, tol, mix)


# ---------------------------------------------------------------------------
# Analytic approximations
# ---------------------------------------------------------------------------

def delta_analytic_eq1(n_density):
    """Analytic Delta(n) from Eq. (1) for Delta0=0, low density limit.

    Valid for Delta << gamma1 and 4pi hbar^2 v^2 |n| << gamma1^2.
    """
    if abs(n_density) < 1e8:
        return 0.0

    cf = cap_factor(1.0)
    lam = screening_lambda(1.0)

    x = hbar * VF * np.sqrt(np.pi * abs(n_density)) / GAMMA1

    numerator = cf * n_density / 2
    denominator = 1 + lam * (x**2 * np.pi * abs(n_density) * hbar**2 * VF**2 / GAMMA1**2
                             - np.log(x))

    # Simpler: from Eq.(1) directly
    # denominator = 1 + Lambda * (hbar^2 v^2 pi|n| / gamma1^2 - ln(hbar*v*sqrt(pi|n|)/gamma1))
    # = 1 + Lambda * (x^2 - ln(x)) where x = hbar*v*sqrt(pi|n|)/gamma1
    denominator = 1 + lam * (x**2 - np.log(x))

    return numerator / denominator


def delta_analytic_eq6(n_density, delta0=0.0, eps_r=1.0):
    """Analytic Delta(n) from Eq. (6), more general than Eq. (1).

    Valid for |Delta|, |Delta0| << epsilon_F < gamma1.
    Breaks down at n=0 for Delta0 != 0.
    """
    if abs(n_density) < 1e8:
        return delta0  # approximate; exact for delta0=0

    cf = cap_factor(eps_r)
    lam = screening_lambda(eps_r)

    p_F = pF_from_density(n_density)
    vp_F = VF * p_F

    # epsilon_0 from Eq. (6) - not the permittivity!
    eps0 = (GAMMA1 / 2) * (np.sqrt(1 + 4 * vp_F**2 / GAMMA1**2) - 1)
    x = eps0 / GAMMA1

    screening = lam * (x + x**2 - 0.5 * np.log(max(x, 1e-30)))

    return (delta0 + cf * n_density / 2) / (1 + screening)


# ---------------------------------------------------------------------------
# Cyclotron mass
# ---------------------------------------------------------------------------

def cyclotron_mass(p, delta, alpha_band=1):
    """Cyclotron mass from Eq. (7).

    alpha_band=1 for low-energy bands, alpha_band=2 for high-energy bands.
    Returns mass in kg.
    """
    vp = VF * p
    g1 = GAMMA1

    eps_sq = (g1**2 / 2 + delta**2 / 4 + vp**2
              + (-1)**alpha_band * np.sqrt(g1**4 / 4 + vp**2 * (g1**2 + delta**2)))
    eps = np.sqrt(max(eps_sq, 0))

    sqrt_term = np.sqrt(vp**2 * (g1**2 + delta**2) + g1**4 / 4)
    bracket = 1 + (-1)**alpha_band * (g1**2 + delta**2) / (2 * sqrt_term)

    if abs(bracket) < 1e-30:
        return np.inf

    return eps / (VF**2 * bracket)


def cyclotron_mass_array(n_arr, delta_arr):
    """Cyclotron mass for arrays of density and corresponding Delta.

    For the low-energy band (alpha=1).
    Returns mass in units of free electron mass m_e.
    """
    mc = np.zeros_like(n_arr, dtype=float)
    for i, (n, d) in enumerate(zip(n_arr, delta_arr)):
        if abs(n) < 1e8:
            # At zero density, mc = gamma1/(2v^2) for Delta=0
            mc[i] = cyclotron_mass(0, d, alpha_band=1) / m_e
        else:
            p_F = pF_from_density(n)
            mc[i] = cyclotron_mass(p_F, d, alpha_band=1) / m_e
    return mc


# ---------------------------------------------------------------------------
# Landau levels
# ---------------------------------------------------------------------------

def landau_levels_low(B, delta, N_max=6, alpha=0, beta=0):
    """Low-energy Landau level spectrum from Eqs. (8, 9).

    Parameters
    ----------
    B : float
        Magnetic field [T]
    delta : float
        Layer asymmetry [J]
    N_max : int
        Maximum LL index
    alpha, beta : float
        Dimensionless asymmetry parameters (typically ~ 0)

    Returns
    -------
    levels : dict
        Keys are (N, xi, sign) tuples, values are energies [J]
    """
    m = GAMMA1 / (2 * VF**2)
    omega_c = e_charge * B / m
    hw = hbar * omega_c

    levels = {}

    for xi in [+1, -1]:
        # N=0 level (Eq. 9)
        e0 = -0.5 * xi * delta - 0.5 * (alpha - beta) * hw
        levels[(0, xi)] = e0

        # N=1 level (Eq. 9)
        e1 = (-0.5 * xi * delta + xi * delta * hw / GAMMA1
              - 0.5 * (3 * alpha - beta) * hw)
        levels[(1, xi)] = e1

        # N >= 2 levels (Eq. 8)
        for N in range(2, N_max + 1):
            e_shift = xi * delta * hw / (2 * GAMMA1)
            e_diag = -alpha * hw * (N - 0.5)
            e_off = hw * np.sqrt(N * (N - 1))

            levels[(N, xi, +1)] = e_shift + e_diag + e_off
            levels[(N, xi, -1)] = e_shift + e_diag - e_off

    return levels
