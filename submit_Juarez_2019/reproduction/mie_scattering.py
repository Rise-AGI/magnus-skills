"""
Core module for Mie scattering and multiple scattering calculations.

Reproduces the scattering approach from Section 2.2 of:
  Juarez, Mendoza, Mochan - "Mie Scattering in the Macroscopic Response
  and the Photonic Bands of Metamaterials"

All figures use the multiple scattering method for an array of dielectric
cylinders in vacuum. Key equations: Eqs. (24)-(33).
"""

import numpy as np
from scipy.special import jv, jvp, hankel1, h1vp


def mie_scattering_coeff(l, n, q, R):
    """
    Compute Mie scattering coefficient s_l (Eq. 26).

    Parameters
    ----------
    l : int
        Angular momentum quantum number.
    n : float
        Refractive index of cylinder.
    q : float
        Free-space wavevector omega/c.
    R : float
        Cylinder radius.

    Returns
    -------
    s_l : complex
        Scattering coefficient b_l/a_l.
    """
    nqR = n * q * R
    qR = q * R
    Jl_nqR = jv(l, nqR)
    Jlp_nqR = jvp(l, nqR)
    Jl_qR = jv(l, qR)
    Jlp_qR = jvp(l, qR)
    Hl_qR = hankel1(l, qR)
    Hlp_qR = h1vp(l, qR)

    num = Jlp_nqR * Jl_qR - n * Jlp_qR * Jl_nqR
    den = n * Jl_nqR * Hlp_qR - Hl_qR * Jlp_nqR
    return num / den


def single_cylinder_coefficients(q, R, n, l_max):
    """
    For a single cylinder illuminated by a plane wave along x,
    the scattering coefficients are b_l = s_l * i^l.

    Returns dict mapping l -> b_l for l in [-l_max, l_max].
    """
    coeffs = {}
    for l in range(-l_max, l_max + 1):
        s_l = mie_scattering_coeff(l, n, q, R)
        a_l = 1j ** l
        coeffs[l] = s_l * a_l
    return coeffs


def lattice_sum(m, q, k_x, N, a, damping=0.0):
    """
    Compute lattice sum S_m for a square lattice of NxN cylinders.

    S_m = sum_{n' != 0} exp(i*m*phi_{0n'}) * H_m(q*R_{0n'}) * exp(i*k*X_n')
                         * exp(-damping * R_{0n'})

    Parameters
    ----------
    m : int
        Order of the Hankel function.
    q : float
        Free-space wavevector.
    k_x : float
        Bloch wavevector along x.
    N : int
        Number of cylinders per side (NxN array).
    a : float
        Lattice constant.
    damping : float
        Exponential damping factor for convergence.

    Returns
    -------
    S_m : complex
        Lattice sum.
    """
    half_N = N // 2
    S = 0.0 + 0.0j
    for n1 in range(-half_N, half_N + 1):
        for n2 in range(-half_N, half_N + 1):
            if n1 == 0 and n2 == 0:
                continue
            X_n = n1 * a
            Y_n = n2 * a
            R_n = np.sqrt(X_n**2 + Y_n**2)
            phi_n = np.arctan2(Y_n, X_n)
            qR = q * R_n
            if qR < 1e-15:
                continue
            H_m = hankel1(m, qR)
            phase = np.exp(1j * k_x * X_n)
            damp = np.exp(-damping * R_n)
            S += np.exp(1j * m * phi_n) * H_m * phase * damp
    return S


def lattice_sums_vectorized(q, k_x, N, a, l_max, damping=0.0):
    """
    Compute all needed lattice sums S_m for m in [-(2*l_max), 2*l_max].
    Uses vectorized numpy operations over the lattice.

    Returns dict mapping m -> S_m.
    """
    half_N = N // 2
    n1 = np.arange(-half_N, half_N + 1)
    n2 = np.arange(-half_N, half_N + 1)
    N1, N2 = np.meshgrid(n1, n2)
    N1 = N1.ravel()
    N2 = N2.ravel()

    # Remove origin
    mask = ~((N1 == 0) & (N2 == 0))
    N1 = N1[mask]
    N2 = N2[mask]

    X = N1 * a
    Y = N2 * a
    R_arr = np.sqrt(X**2 + Y**2)
    phi_arr = np.arctan2(Y, X)
    qR_arr = q * R_arr

    phase_arr = np.exp(1j * k_x * X)
    damp_arr = np.exp(-damping * R_arr)

    sums = {}
    for m in range(-2 * l_max, 2 * l_max + 1):
        H_m = hankel1(m, qR_arr)
        angular = np.exp(1j * m * phi_arr)
        S_m = np.sum(angular * H_m * phase_arr * damp_arr)
        sums[m] = S_m
    return sums


def solve_bloch_scattering(q, k_x, R, n_refr, a, N, l_max, damping=0.0):
    """
    Solve the multiple scattering problem for an NxN array of cylinders
    using Bloch's theorem (all cylinders scatter identically up to a phase).

    Solves: beta_l - s_l * sum_l' S_{l-l'} * beta_l' = s_l * i^l

    Parameters
    ----------
    q : float
        Free-space wavevector omega/c.
    k_x : float
        Bloch wavevector along x.
    R : float
        Cylinder radius.
    n_refr : float
        Refractive index.
    a : float
        Lattice constant.
    N : int
        Array size (NxN).
    l_max : int
        Max angular momentum.
    damping : float
        Damping factor for lattice sums.

    Returns
    -------
    beta : dict
        Scattering coefficients beta_l for l in [-l_max, l_max].
    """
    ls = list(range(-l_max, l_max + 1))
    n_l = len(ls)

    # Compute scattering coefficients
    s = {}
    for l in ls:
        s[l] = mie_scattering_coeff(l, n_refr, q, R)

    # Compute lattice sums
    lsums = lattice_sums_vectorized(q, k_x, N, a, l_max, damping)

    # Build system: (I - S_mat) beta = rhs
    # T_{ll'} = delta_{ll'} - s_l * S_{l-l'}  (for l != 0 in lattice, but here S already excludes origin)
    # Actually: T beta = rhs where T_{ll'} = delta_{ll'} - s_l * S_{l-l'}
    T_mat = np.eye(n_l, dtype=complex)
    rhs = np.zeros(n_l, dtype=complex)

    for i, l in enumerate(ls):
        rhs[i] = s[l] * (1j ** l)
        for j, lp in enumerate(ls):
            m = l - lp
            if m in lsums:
                T_mat[i, j] -= s[l] * lsums[m]

    beta_arr = np.linalg.solve(T_mat, rhs)
    beta = {}
    for i, l in enumerate(ls):
        beta[l] = beta_arr[i]
    return beta


def compute_scattered_coefficients(nqR_array, R_over_a, n_refr, N_array, l_max,
                                    k_factor=1.01, damping_per_a=0.0):
    """
    Compute scattering coefficients beta_l for a range of frequencies.

    Parameters
    ----------
    nqR_array : array
        Values of n*q*R to scan.
    R_over_a : float
        R/a ratio.
    n_refr : float
        Refractive index.
    N_array : int
        Array size.
    l_max : int
        Max angular momentum.
    k_factor : float
        k/q ratio (k = k_factor * q).
    damping_per_a : float
        Damping factor per lattice constant.

    Returns
    -------
    results : dict
        For each l, array of beta_l values.
    """
    a = 1.0  # Lattice constant (normalized)
    R = R_over_a * a

    results = {l: np.zeros(len(nqR_array), dtype=complex) for l in range(-l_max, l_max + 1)}

    for idx, nqR in enumerate(nqR_array):
        q = nqR / (n_refr * R)
        k_x = k_factor * q
        damping = damping_per_a / a

        try:
            beta = solve_bloch_scattering(q, k_x, R, n_refr, a, N_array, l_max, damping)
            for l in range(-l_max, l_max + 1):
                results[l][idx] = beta[l]
        except (np.linalg.LinAlgError, ValueError):
            pass

    return results
