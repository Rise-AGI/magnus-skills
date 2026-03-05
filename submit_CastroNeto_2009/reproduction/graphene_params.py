"""
Core module for graphene electronic properties calculations.

Implements tight-binding models, density of states, Klein tunneling,
and nanoribbon spectra from Castro Neto et al., Rev. Mod. Phys. 81, 109 (2009).
"""

import numpy as np
from scipy.special import ellipk


# ============================================================
# Physical constants and default parameters
# ============================================================

A_CC = 1.42e-10      # carbon-carbon distance in meters (1.42 Angstrom)
A_LAT = A_CC * np.sqrt(3)  # lattice constant ~ 2.46 Angstrom
T_HOP = 2.7          # nearest-neighbor hopping energy in eV
T_PRIME = 0.2 * T_HOP  # next-nearest-neighbor hopping in eV
GAMMA1 = 0.4         # bilayer interlayer hopping gamma1 in eV
VF = 3 * T_HOP * A_CC / (2 * 1.0545718e-34) * 1.602e-19  # Fermi velocity in m/s


# ============================================================
# Lattice vectors
# ============================================================

def lattice_vectors():
    """Return lattice vectors a1, a2 and reciprocal vectors b1, b2."""
    a = A_CC
    a1 = a / 2 * np.array([3.0, np.sqrt(3)])
    a2 = a / 2 * np.array([3.0, -np.sqrt(3)])
    b1 = 2 * np.pi / (3 * a) * np.array([1.0, np.sqrt(3)])
    b2 = 2 * np.pi / (3 * a) * np.array([1.0, -np.sqrt(3)])
    return a1, a2, b1, b2


def dirac_points():
    """Return K and K' points."""
    a = A_CC
    K = np.array([2 * np.pi / (3 * a), 2 * np.pi / (3 * np.sqrt(3) * a)])
    Kp = np.array([2 * np.pi / (3 * a), -2 * np.pi / (3 * np.sqrt(3) * a)])
    return K, Kp


# ============================================================
# Single-layer tight-binding dispersion
# ============================================================

def f_function(kx, ky, a=A_CC):
    """Compute f(k) = 2cos(sqrt(3)*ky*a) + 4cos(sqrt(3)/2*ky*a)*cos(3/2*kx*a)."""
    return (2 * np.cos(np.sqrt(3) * ky * a)
            + 4 * np.cos(np.sqrt(3) / 2 * ky * a) * np.cos(3.0 / 2 * kx * a))


def energy_dispersion(kx, ky, t=T_HOP, tp=T_PRIME, a=A_CC):
    """
    Tight-binding energy dispersion E_pm(k), Eq. (6).

    Returns (E_plus, E_minus) arrays.
    """
    fk = f_function(kx, ky, a)
    sqrt_term = np.sqrt(3 + fk)
    E_plus = t * sqrt_term - tp * fk
    E_minus = -t * sqrt_term - tp * fk
    return E_plus, E_minus


def energy_dispersion_full(kx, ky, t=T_HOP, tp=0.0, a=A_CC):
    """
    Full tight-binding dispersion, Eq. (6), for arbitrary t'.
    Returns (E_plus, E_minus).
    """
    fk = f_function(kx, ky, a)
    # Clamp to avoid sqrt of negative due to numerical noise
    arg = np.clip(3 + fk, 0, None)
    sqrt_term = np.sqrt(arg)
    E_plus = t * sqrt_term - tp * fk
    E_minus = -t * sqrt_term - tp * fk
    return E_plus, E_minus


# ============================================================
# Density of states (analytical, Eq. 14)
# ============================================================

def dos_analytical(E_array, t=T_HOP):
    """
    Analytical density of states per unit cell for t'=0, Eq. (14).
    Uses complete elliptic integral of the first kind.
    """
    rho = np.zeros_like(E_array, dtype=float)

    for i, E in enumerate(E_array):
        eps = abs(E) / t
        if eps < 1e-12 or eps > 3.0 - 1e-12:
            rho[i] = 0.0
            continue

        if eps <= 1.0:
            Z0 = (1 + eps)**2 - ((eps**2 - 1)**2) / 4
            Z1 = 4 * eps
        else:  # 1 < eps < 3
            Z0 = 4 * eps
            Z1 = (1 + eps)**2 - ((eps**2 - 1)**2) / 4

        if Z0 <= 0 or Z1 / Z0 < 0 or Z1 / Z0 > 1:
            rho[i] = 0.0
            continue

        k_arg = np.sqrt(Z1 / Z0)
        K_val = ellipk(k_arg**2)  # ellipk takes m = k^2
        rho[i] = (4 / (np.pi**2)) * (abs(E) / t**2) * (1 / np.sqrt(Z0)) * K_val

    return rho


def dos_numerical(kx_grid, ky_grid, t=T_HOP, tp=0.0, n_bins=500, E_range=(-3.5, 3.5)):
    """
    Compute DOS numerically by histogramming the energy spectrum.
    """
    Ep, Em = energy_dispersion_full(kx_grid, ky_grid, t, tp)
    energies = np.concatenate([Ep.ravel(), Em.ravel()])

    bins = np.linspace(E_range[0] * t, E_range[1] * t, n_bins + 1)
    hist, bin_edges = np.histogram(energies, bins=bins, density=True)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    # Normalize: DOS per unit cell, factor of 4 for spin and valley degeneracy already in spectrum
    return bin_centers / t, hist * t  # return in units of t


# ============================================================
# Klein tunneling, Eq. (29)
# ============================================================

def klein_tunneling_T(phi_array, E, V0, D, vF_eff=1e6):
    """
    Transmission probability T(phi) for Dirac electrons through a square barrier.
    Eq. (29) of the paper.

    Parameters:
        phi_array: incidence angles in radians
        E: electron energy in eV
        V0: barrier height in eV
        D: barrier width in meters
        vF_eff: Fermi velocity in m/s

    Returns:
        T(phi) array
    """
    hbar = 1.0545718e-34  # J*s
    eV_to_J = 1.602e-19

    E_J = E * eV_to_J
    V0_J = V0 * eV_to_J

    kF = abs(E_J) / (hbar * vF_eff)

    s = np.sign(E)
    sp = np.sign(E - V0)

    T_arr = np.zeros_like(phi_array, dtype=float)

    for i, phi in enumerate(phi_array):
        kx = kF * np.cos(phi)
        ky = kF * np.sin(phi)

        qx_sq = ((V0_J - E_J) / (hbar * vF_eff))**2 - ky**2

        if qx_sq < 0:
            # Evanescent wave - very low transmission
            T_arr[i] = 0.0
            continue

        qx = np.sqrt(qx_sq)
        theta = np.arctan2(ky, qx)

        cos_phi = np.cos(phi)
        cos_theta = np.cos(theta)
        sin_phi = np.sin(phi)
        sin_theta = np.sin(theta)

        cos_Dqx = np.cos(D * qx)
        sin_Dqx = np.sin(D * qx)

        numerator = cos_theta**2 * cos_phi**2
        denominator = (cos_Dqx * cos_phi * cos_theta)**2 + sin_Dqx**2 * (1 - s * sp * sin_phi * sin_theta)**2

        if denominator < 1e-30:
            T_arr[i] = 1.0
        else:
            T_arr[i] = numerator / denominator

    return T_arr


# ============================================================
# Nanoribbon spectra (tight-binding)
# ============================================================

def zigzag_nanoribbon_spectrum(N, k_array):
    """
    Compute energy spectrum of a zigzag graphene nanoribbon.

    Parameters:
        N: number of zigzag chains (ribbon width)
        k_array: momentum along ribbon direction (in units of 1/(sqrt(3)*a))

    Returns:
        energies: array of shape (len(k_array), 2*N)
    """
    t = T_HOP
    energies = np.zeros((len(k_array), 2 * N))

    for ik, k in enumerate(k_array):
        # Build Hamiltonian for zigzag ribbon
        H = np.zeros((2 * N, 2 * N))

        for n in range(N):
            # A site index: 2*n, B site index: 2*n+1
            # Intra-cell A-B coupling
            H[2*n, 2*n+1] = -t * (1 + np.exp(-1j * k * np.sqrt(3)))
            H[2*n+1, 2*n] = -t * (1 + np.exp(1j * k * np.sqrt(3)))

            # Inter-cell B-A coupling (B_n to A_{n+1})
            if n < N - 1:
                H[2*n+1, 2*(n+1)] = -t
                H[2*(n+1), 2*n+1] = -t

        # Ensure Hermiticity
        H = (H + H.conj().T) / 2

        evals = np.linalg.eigvalsh(H)
        energies[ik, :] = np.sort(evals)

    return energies


def armchair_nanoribbon_spectrum(N, k_array):
    """
    Compute energy spectrum of an armchair graphene nanoribbon.

    Parameters:
        N: number of dimer lines across ribbon
        k_array: momentum along ribbon direction (in units of 1/a)

    Returns:
        energies: array of shape (len(k_array), N)
    """
    t = T_HOP
    energies = np.zeros((len(k_array), 2 * N))

    for ik, k in enumerate(k_array):
        H = np.zeros((2 * N, 2 * N), dtype=complex)

        for n in range(N):
            # A site: 2*n, B site: 2*n+1
            # A-B hopping within unit cell
            H[2*n, 2*n+1] = -t * (1 + np.exp(-1j * k))
            H[2*n+1, 2*n] = -t * (1 + np.exp(1j * k))

            # A_n - B_{n-1} and B_n - A_{n+1} coupling
            if n > 0:
                H[2*n, 2*(n-1)+1] = -t
                H[2*(n-1)+1, 2*n] = -t

        H = (H + H.conj().T) / 2
        evals = np.linalg.eigvalsh(H)
        energies[ik, :] = np.sort(evals)

    return energies


# ============================================================
# Bilayer graphene band structure
# ============================================================

def bilayer_bands(k_array, V=0.0, gamma1=GAMMA1, gamma3=0.0, t=T_HOP):
    """
    Compute bilayer graphene band structure.
    4x4 Hamiltonian, Eq. (38).

    Parameters:
        k_array: 1D array of |k| values (in 1/m or arbitrary units matching vF)
        V: interlayer potential asymmetry (half shift)
        gamma1: interlayer hopping (eV)
        gamma3: trigonal warping hopping (eV)
        t: intralayer hopping (eV)

    Returns:
        bands: array of shape (len(k_array), 4)
    """
    a = A_CC
    vF = 3 * t * a / 2  # in eV*m (natural units)

    bands = np.zeros((len(k_array), 4))

    for i, kval in enumerate(k_array):
        # Use k along kx direction
        kc = kval + 0j  # complex k = kx + i*ky

        H = np.array([
            [-V, vF * kc, 0, 3 * gamma3 * a * np.conj(kc)],
            [vF * np.conj(kc), -V, gamma1, 0],
            [0, gamma1, V, vF * kc],
            [3 * gamma3 * a * kc, 0, vF * np.conj(kc), V]
        ], dtype=complex)

        evals = np.linalg.eigvalsh(H)
        bands[i, :] = np.sort(evals.real)

    return bands
