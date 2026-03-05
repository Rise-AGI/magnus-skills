"""
Core physics module for the Lattice Schwinger Model.

Implements:
- One-particle energy spectra for various fermion derivatives
- Anomalous commutator W(k) computation
- Coupled dispersion relation eigenvalues (modified SLAC)
- Perfect Wilson derivative spectrum

Reference: Melnikov & Weinstein, "The Lattice Schwinger Model:
Confinement, Anomalies, Chiral Fermions and All That" (2000)
hep-lat/0006029
"""

import numpy as np


# =============================================================
# Fermion derivative definitions
# =============================================================

def naive_derivative(xi):
    """Naive fermion derivative: Z = sin(xi), X = 0."""
    return np.sin(xi), np.zeros_like(xi)


def wilson_derivative(xi, r=1.0):
    """Wilson fermion derivative: Z = sin(xi), X = r(1 - cos(xi))."""
    return np.sin(xi), r * (1.0 - np.cos(xi))


def slac_derivative(xi):
    """SLAC derivative: Z = xi, X = 0 (for xi in [0, pi])."""
    return xi.copy(), np.zeros_like(xi)


def modified_slac_derivative(xi, mu):
    """
    Modified SLAC derivative (Eq. 86).
    Z_k = k * theta(mu*pi - k) + mu/(mu-1) * (pi - k) * theta(k - mu*pi)
    X_k = 0
    Here xi = k*a is in [0, pi].
    """
    Z = np.where(
        xi < mu * np.pi,
        xi,
        (mu / (mu - 1.0)) * (np.pi - xi)
    )
    return Z, np.zeros_like(xi)


def perfect_wilson_derivative(xi):
    """
    Perfect Wilson derivative (Eq. 100).
    For k > 0 (xi = ka in [0, pi]):
    Z = xi for xi < pi/2
    Z = (pi/2) * sin(pi - xi) for xi >= pi/2
    X = 0 for xi < pi/2
    X = (pi/2) * cos(pi - xi) for xi >= pi/2
    """
    Z = np.where(
        xi < np.pi / 2,
        xi,
        (np.pi / 2) * np.sin(np.pi - xi)
    )
    X = np.where(
        xi < np.pi / 2,
        0.0,
        (np.pi / 2) * np.cos(np.pi - xi)
    )
    return Z, X


# =============================================================
# Energy spectrum
# =============================================================

def energy_spectrum(Z, X):
    """One-particle energy E_k = sqrt(Z_k^2 + X_k^2)."""
    return np.sqrt(Z**2 + X**2)


# =============================================================
# Anomalous commutator (Eq. 83)
# =============================================================

def anomalous_commutator_continuum_limit(derivative_func, *args, npts=10000):
    """
    Compute the continuum limit (ak -> 0) of the anomalous commutator W.

    W = k/pi * integral_0^pi [d^2Z/dxi^2 * cos(theta) + d^2X/dxi^2 * sin(theta)] dxi

    Returns W / (k/pi) so the result is the dimensionless prefactor.
    For the correct continuum Schwinger model, this should be -1.

    NOTE: Only works correctly for smooth derivative functions (naive, Wilson).
    For piecewise-linear derivatives (modified SLAC, perfect Wilson),
    use anomalous_commutator_analytical() instead.
    """
    xi = np.linspace(0, np.pi, npts)
    dxi = xi[1] - xi[0]

    Z, X = derivative_func(xi, *args)
    E = energy_spectrum(Z, X)

    E = np.where(E < 1e-15, 1e-15, E)
    cos_theta = Z / E
    sin_theta = X / E

    d2Z = np.gradient(np.gradient(Z, dxi), dxi)
    d2X = np.gradient(np.gradient(X, dxi), dxi)

    integrand = d2Z * cos_theta + d2X * sin_theta
    result = np.trapezoid(integrand, xi)

    return result


def anomalous_commutator_analytical(derivative_name, **kwargs):
    """
    Analytical results for the anomalous commutator W/(k/pi) from the paper.

    Returns the dimensionless coefficient such that W = coefficient * k/pi.
    """
    if derivative_name == "naive":
        return -2.0
    elif derivative_name == "wilson":
        return -2.0  # r-independent in a->0 limit (Eq. 85)
    elif derivative_name == "slac":
        return -1.0  # Correct continuum value
    elif derivative_name == "modified_slac":
        mu = kwargs.get("mu", 0.5)
        return -(1.0 + mu / (1.0 - mu))  # Eq. (88): -(1 + c) where c=mu/(1-mu)
    elif derivative_name == "perfect_wilson":
        return -2.0  # Differs from continuum (discussed after Eq. 100)
    else:
        raise ValueError(f"Unknown derivative: {derivative_name}")


def anomalous_commutator_full(derivative_func, k_vals, a=1.0, *args, npts=2000):
    """
    Full anomalous commutator W(k) for finite lattice spacing a.

    From Eq. (82):
    W(k) = (a * e^{-ika/2}) / (2 sin(ka/2)) * integral

    Returns W(k) for each k value.
    """
    results = []
    xi_pts = np.linspace(-np.pi, np.pi, npts, endpoint=False)
    dxi = xi_pts[1] - xi_pts[0]

    Z_full, X_full = derivative_func(np.abs(xi_pts), *args)
    # Make Z odd and X even for full range
    Z_full = np.sign(xi_pts) * Z_full

    E_full = energy_spectrum(Z_full, X_full)
    E_full = np.where(E_full < 1e-15, 1e-15, E_full)
    cos_theta = Z_full / E_full
    sin_theta = X_full / E_full

    for k in k_vals:
        ka = k * a
        if abs(np.sin(ka / 2)) < 1e-10:
            results.append(0.0)
            continue

        # Integrand from Eq. (82)
        integrand = (
            (2 * Z_full - np.roll(Z_full, int(k * npts / (2 * np.pi / a)))
             - np.roll(Z_full, -int(k * npts / (2 * np.pi / a)))) * cos_theta
            + (2 * X_full - np.roll(X_full, int(k * npts / (2 * np.pi / a)))
               - np.roll(X_full, -int(k * npts / (2 * np.pi / a)))) * sin_theta
        )

        prefactor = a * np.exp(-1j * ka / 2) / (2 * np.sin(ka / 2))
        W = prefactor * np.trapezoid(integrand, xi_pts) / (2 * np.pi)
        results.append(np.real(W))

    return np.array(results)


# =============================================================
# Coupled dispersion relation (Eq. 93)
# =============================================================

def coupled_dispersion_eigenvalues(k_vals, c, e2_over_pi=1.0):
    """
    Eigenvalues of the coupled system matrix M (Eq. 93):

    M = [[k^2 + e^2/pi,      sqrt(c)*e^2/pi],
         [sqrt(c)*e^2/pi,  c^2*k^2 + c*e^2/pi]]

    Returns (E1, E2) for each k value where E = sqrt(eigenvalue).
    """
    E1_list = []
    E2_list = []

    for k in k_vals:
        M = np.array([
            [k**2 + e2_over_pi, np.sqrt(c) * e2_over_pi],
            [np.sqrt(c) * e2_over_pi, c**2 * k**2 + c * e2_over_pi]
        ])
        eigenvalues = np.sort(np.linalg.eigvalsh(M))
        E1_list.append(np.sqrt(max(eigenvalues[0], 0)))
        E2_list.append(np.sqrt(max(eigenvalues[1], 0)))

    return np.array(E1_list), np.array(E2_list)


def uncoupled_dispersion(k_vals, c, e2_over_pi=1.0):
    """
    Uncoupled dispersion relations (Eq. 94), valid for c*k^2 >> e^2/pi.
    E1 = sqrt(k^2 + e^2/pi)
    E2 = sqrt(c^2*k^2 + c*e^2/pi)
    """
    E1 = np.sqrt(k_vals**2 + e2_over_pi)
    E2 = np.sqrt(c**2 * k_vals**2 + c * e2_over_pi)
    return E1, E2


def small_k_dispersion(k_vals, c, e2_over_pi=1.0):
    """
    Small momentum limit (Eq. 92), valid for c*k^2 << e^2/pi.
    E1 = sqrt(c) * |k|  (massless Goldstone mode)
    E2 = sqrt(1+c) * sqrt(e^2/pi)  (massive mode)
    """
    E1 = np.sqrt(c) * np.abs(k_vals)
    E2 = np.sqrt((1 + c) * e2_over_pi) * np.ones_like(k_vals)
    return E1, E2
