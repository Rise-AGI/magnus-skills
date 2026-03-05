"""
Core module for Anderson localization computations.

Implements the Anderson tight-binding Hamiltonian on d-dimensional lattices
and provides routines for:
  - Hamiltonian construction (orthogonal, unitary, symplectic)
  - Exact diagonalization (eigenvalues + eigenvectors)
  - Transfer matrix method for quasi-1D conductance
  - Density of states, inverse participation ratio, level statistics
"""

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh
from scipy.linalg import eigh


# ---------------------------------------------------------------------------
# Hamiltonian construction
# ---------------------------------------------------------------------------

def anderson_hamiltonian_3d(L, W, disorder="box", V=1.0, periodic=True, seed=None):
    """
    Build the 3D Anderson Hamiltonian on an L x L x L cubic lattice.

    Parameters
    ----------
    L : int  — linear system size
    W : float — disorder strength
    disorder : str — "box" or "gaussian"
    V : float — hopping amplitude (uniform in all directions)
    periodic : bool — periodic boundary conditions
    seed : int or None

    Returns
    -------
    H : ndarray of shape (N, N), N = L**3
    """
    rng = np.random.default_rng(seed)
    N = L ** 3

    # On-site energies
    if disorder == "box":
        diag = W * rng.uniform(-0.5, 0.5, size=N)
    else:  # gaussian
        diag = W * rng.normal(0, 1, size=N)

    # Build sparse Hamiltonian
    row, col, val = [], [], []

    def idx(x, y, z):
        return (x % L) * L * L + (y % L) * L + (z % L)

    for x in range(L):
        for y in range(L):
            for z in range(L):
                i = idx(x, y, z)
                # Neighbors in +x, +y, +z
                for dx, dy, dz in [(1,0,0),(0,1,0),(0,0,1)]:
                    nx_, ny_, nz_ = x+dx, y+dy, z+dz
                    if periodic or (0 <= nx_ < L and 0 <= ny_ < L and 0 <= nz_ < L):
                        j = idx(nx_, ny_, nz_)
                        row.extend([i, j])
                        col.extend([j, i])
                        val.extend([V, V])

    H = sparse.csr_matrix((val, (row, col)), shape=(N, N), dtype=np.float64)
    H.setdiag(diag)
    return H


def anderson_hamiltonian_2d(L, W, disorder="box", V=1.0, periodic=True, seed=None):
    """
    Build the 2D Anderson Hamiltonian on an L x L square lattice.
    """
    rng = np.random.default_rng(seed)
    N = L * L

    if disorder == "box":
        diag = W * rng.uniform(-0.5, 0.5, size=N)
    else:
        diag = W * rng.normal(0, 1, size=N)

    row, col, val = [], [], []

    def idx(x, y):
        return (x % L) * L + (y % L)

    for x in range(L):
        for y in range(L):
            i = idx(x, y)
            for dx, dy in [(1,0),(0,1)]:
                nx_, ny_ = x+dx, y+dy
                if periodic or (0 <= nx_ < L and 0 <= ny_ < L):
                    j = idx(nx_, ny_)
                    row.extend([i, j])
                    col.extend([j, i])
                    val.extend([V, V])

    H = sparse.csr_matrix((val, (row, col)), shape=(N, N), dtype=np.float64)
    H.setdiag(diag)
    return H


# ---------------------------------------------------------------------------
# Density of states
# ---------------------------------------------------------------------------

def density_of_states(eigenvalues, bins=200, energy_range=None):
    """Compute density of states histogram from eigenvalues."""
    if energy_range is None:
        energy_range = (eigenvalues.min() - 0.5, eigenvalues.max() + 0.5)
    counts, edges = np.histogram(eigenvalues, bins=bins, range=energy_range, density=True)
    centers = 0.5 * (edges[:-1] + edges[1:])
    return centers, counts


# ---------------------------------------------------------------------------
# Level spacing statistics
# ---------------------------------------------------------------------------

def level_spacing_ratios(eigenvalues, energy_window=None):
    """
    Compute normalized nearest-neighbor level spacings.

    Parameters
    ----------
    eigenvalues : sorted array
    energy_window : tuple (E_min, E_max) or None

    Returns
    -------
    spacings : normalized spacings s with <s> = 1
    """
    E = np.sort(eigenvalues)
    if energy_window is not None:
        mask = (E >= energy_window[0]) & (E <= energy_window[1])
        E = E[mask]
    s = np.diff(E)
    s = s / np.mean(s)  # normalize so <s> = 1
    return s


def wigner_surmise(s, beta=1):
    """Wigner surmise for GOE (beta=1), GUE (beta=2), GSE (beta=4)."""
    if beta == 1:
        return (np.pi / 2) * s * np.exp(-np.pi * s**2 / 4)
    elif beta == 2:
        return (32 / np.pi**2) * s**2 * np.exp(-4 * s**2 / np.pi)
    elif beta == 4:
        return (2**18 / (3**6 * np.pi**3)) * s**4 * np.exp(-64 * s**2 / (9 * np.pi))


# ---------------------------------------------------------------------------
# Transfer matrix for quasi-1D conductance
# ---------------------------------------------------------------------------

def transfer_matrix_1d(Lz, W, E, disorder="box", seed=None):
    """
    Compute the transfer matrix product for a 1D Anderson chain.

    Returns x = arccosh(1/sqrt(g)) parameter where g = 1/cosh^2(x/2).
    Actually returns x such that g = 1/cosh^2(x/2).
    """
    rng = np.random.default_rng(seed)
    if disorder == "box":
        eps = W * rng.uniform(-0.5, 0.5, size=Lz)
    else:
        eps = W * rng.normal(0, 1, size=Lz)

    # Transfer matrix product: T_n = [[E-eps_n, -1], [1, 0]]
    # We track T = prod T_n
    T = np.eye(2)
    for n in range(Lz):
        Tn = np.array([[E - eps[n], -1.0], [1.0, 0.0]])
        T = Tn @ T

    # Conductance from transfer matrix:
    # g = 1 / |T[0,0]|^2  (for single channel)
    # Actually: g = 4 / (Tr(T^dagger T) + 2) for real T
    # More precisely: eigenvalues of T^T T give exp(2*gamma*L)
    # Use: g = 1/cosh^2(x/2) where cosh(x/2) = |T[0,0]|/sqrt(det T)
    # For 1D: det(T_n) = 1, so det(T) = 1
    # Then cosh^2(x/2) = (Tr(T^T T) + 2) / 4

    TtT = T.T @ T
    trace_val = np.trace(TtT)
    # cosh^2(x/2) = (trace_val + 2) / 4
    cosh2 = (trace_val + 2) / 4.0
    if cosh2 < 1.0:
        cosh2 = 1.0
    x = 2 * np.arccosh(np.sqrt(cosh2))
    g = 1.0 / cosh2
    return x, g


def transfer_matrix_quasi1d(Lz, L, W, E, disorder="box", V=1.0, t=1.0, seed=None):
    """
    Transfer matrix method for quasi-1D strip (2D: L x Lz) or bar (3D: L x L x Lz).
    Uses the recursive Green's function / transfer matrix approach.

    For a 2D strip of width L and length Lz:
    Returns the Landauer conductance g = Tr(t^dag t) in units of e^2/h.
    """
    rng = np.random.default_rng(seed)
    N = L  # number of channels (for 2D strip)

    # Hopping matrix between slices
    T_hop = t * np.eye(N)

    # Build inter-slice hopping (periodic in transverse direction)
    H_perp = np.zeros((N, N))
    for i in range(N):
        j = (i + 1) % N
        H_perp[i, j] = V
        H_perp[j, i] = V

    # Green's function recursive method
    # Initialize
    g_ret = np.zeros((N, N), dtype=complex)

    # Lead self-energy (semi-infinite lead)
    # For simplicity, use iterative method
    # Actually, let's use the transfer matrix approach directly

    # Transfer matrix: M_total = prod_n M_n
    # M_n = [[(E - H_slice_n) * T_hop^-1, -T_hop^T], [T_hop^-1, 0]]
    # This is N x N block matrix, total 2N x 2N

    T_hop_inv = np.linalg.inv(T_hop)
    M = np.eye(2 * N)

    for n in range(Lz):
        if disorder == "box":
            eps = W * rng.uniform(-0.5, 0.5, size=N)
        else:
            eps = W * rng.normal(0, 1, size=N)

        H_slice = H_perp.copy()
        H_slice[np.diag_indices(N)] += eps

        A = (E * np.eye(N) - H_slice) @ T_hop_inv
        B = -T_hop.T
        C = T_hop_inv
        D = np.zeros((N, N))

        Mn = np.block([[A, B], [C, D]])
        M = Mn @ M

        # Stabilize via QR periodically
        if (n + 1) % 10 == 0 and n < Lz - 1:
            Q, R = np.linalg.qr(M)
            M = Q  # We lose the R factor, but this is for Lyapunov exponents

    # Extract Lyapunov exponents from M
    # This simple approach is numerically unstable for large Lz
    # For proper implementation, use QR decomposition at each step
    # For now, return trace-based conductance estimate
    svd_vals = np.linalg.svd(M[:N, :N], compute_uv=False)
    g = np.sum(1.0 / (svd_vals**2 + 1e-30))
    return g


def transfer_matrix_1d_stable(Lz, W, E, disorder="box", seed=None):
    """
    Stable 1D transfer matrix using the Ricatti variable.
    Returns x = 2*Lz/lambda (Lyapunov exponent times length).
    """
    rng = np.random.default_rng(seed)
    if disorder == "box":
        eps = W * rng.uniform(-0.5, 0.5, size=Lz)
    else:
        eps = W * rng.normal(0, 1, size=Lz)

    # Use log of transfer matrix norm for stability
    log_norm = 0.0
    v = np.array([1.0, 0.0])  # initial vector

    for n in range(Lz):
        Tn = np.array([[E - eps[n], -1.0], [1.0, 0.0]])
        v = Tn @ v
        norm_v = np.linalg.norm(v)
        log_norm += np.log(norm_v)
        v /= norm_v

    # The Lyapunov exponent gamma = log_norm / Lz
    # x = 2 * gamma * Lz = 2 * log_norm
    # But x is defined via g = 1/cosh^2(x/2)
    # For large systems, gamma*Lz >> 1, so x ~ 2*gamma*Lz
    x = 2 * log_norm / Lz  # This is 2*gamma, multiply by Lz to get x
    return log_norm, log_norm  # return (gamma*Lz, gamma*Lz) for consistency


def conductance_1d_ensemble(Lz, W, E, n_samples, disorder="box", seed=42):
    """
    Compute conductance statistics for 1D Anderson chain.
    Returns arrays of x values and g values.
    """
    rng = np.random.default_rng(seed)
    x_vals = np.zeros(n_samples)
    g_vals = np.zeros(n_samples)

    for i in range(n_samples):
        x, g = transfer_matrix_1d(Lz, W, E, disorder=disorder, seed=rng.integers(0, 2**31))
        x_vals[i] = x
        g_vals[i] = g

    return x_vals, g_vals


def mean_x_vs_length(lengths, W, E, n_samples=1000, disorder="box", seed=42):
    """
    Compute mean <x> as a function of system length for 1D chain.
    """
    rng = np.random.default_rng(seed)
    mean_x = np.zeros(len(lengths))

    for il, Lz in enumerate(lengths):
        x_sum = 0.0
        for i in range(n_samples):
            x, g = transfer_matrix_1d(Lz, W, E, disorder=disorder,
                                       seed=rng.integers(0, 2**31))
            x_sum += x
        mean_x[il] = x_sum / n_samples

    return mean_x


# ---------------------------------------------------------------------------
# 2D conductance via transfer matrix (proper implementation)
# ---------------------------------------------------------------------------

def conductance_2d_strip(L, Lz, W, E, disorder="box", seed=None, V=1.0):
    """
    Compute conductance of a 2D strip (width L, length Lz) using
    the transfer matrix method with QR stabilization.

    Returns g in units of e^2/h.
    """
    rng = np.random.default_rng(seed)
    N = L  # number of transverse channels

    # Transverse hopping matrix (periodic BC)
    H_perp = np.zeros((N, N))
    for i in range(N):
        j = (i + 1) % N
        H_perp[i, j] = V
        H_perp[j, i] = V

    # Transfer matrix approach with QR decomposition for stability
    # We track the product of transfer matrices via QR
    Q = np.eye(2 * N)
    log_R_diag = np.zeros(2 * N)

    for n in range(Lz):
        if disorder == "box":
            eps = W * rng.uniform(-0.5, 0.5, size=N)
        else:
            eps = W * rng.normal(0, 1, size=N)

        H_slice = H_perp + np.diag(eps)
        A = E * np.eye(N) - H_slice
        I_N = np.eye(N)

        # Transfer matrix for this slice:
        # [psi_{n+1}]   [A  -I] [psi_n  ]
        # [psi_n    ] = [I   0] [psi_{n-1}]
        Mn = np.block([[A, -I_N], [I_N, np.zeros((N, N))]])

        Q = Mn @ Q
        # QR stabilization every step
        Q, R = np.linalg.qr(Q)
        log_R_diag += np.log(np.abs(np.diag(R)) + 1e-300)

    # Lyapunov exponents = log_R_diag / Lz (sorted)
    lyap = np.sort(log_R_diag / Lz)[::-1]

    # For the conductance, we need the transmission eigenvalues
    # The N smallest Lyapunov exponents (by magnitude) give the
    # localization lengths: lambda_a = 1 / gamma_a
    # g = sum_a 1/cosh^2(gamma_a * Lz)

    # The 2N Lyapunov exponents come in +/- pairs
    # Take the N positive ones (smallest N from the sorted list)
    gamma_pos = lyap[:N]
    g = np.sum(1.0 / np.cosh(gamma_pos)**2)

    return g


# ---------------------------------------------------------------------------
# Inverse participation ratio
# ---------------------------------------------------------------------------

def inverse_participation_ratio(psi):
    """Compute IPR = sum |psi|^4 for a normalized eigenvector."""
    return np.sum(np.abs(psi)**4)


# ---------------------------------------------------------------------------
# Analytical predictions
# ---------------------------------------------------------------------------

def dos_3d_clean(E, V=1.0, n_points=500):
    """
    Compute density of states for clean 3D cubic lattice
    by sampling k-space.
    """
    k = np.linspace(-np.pi, np.pi, n_points)
    kx, ky, kz = np.meshgrid(k, k, k, indexing='ij')
    energies = 2*V*(np.cos(kx) + np.cos(ky) + np.cos(kz))
    return energies.ravel()
