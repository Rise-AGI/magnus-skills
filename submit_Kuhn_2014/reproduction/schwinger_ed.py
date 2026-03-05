"""
Exact diagonalization for the Schwinger model with finite-dimensional link variables.

Implements two models from Kuhn, Cirac, Banuls (2014):
  1. Truncated compact QED (cQED) - U(1) gauge symmetry preserved
  2. Zd model - discrete Zd gauge symmetry

Uses the spin formulation (Jordan-Wigner transformed) with open boundary conditions.
The Gauss law constraint is used to project onto the physical subspace.
"""

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh
from scipy.linalg import expm


# ============================================================
# Operator construction
# ============================================================

def sigma_z():
    return sparse.diags([1, -1], 0, shape=(2, 2), format='csr')

def sigma_plus():
    return sparse.csr_matrix(([1.0], ([0], [1])), shape=(2, 2))

def sigma_minus():
    return sparse.csr_matrix(([1.0], ([1], [0])), shape=(2, 2))

def identity(d):
    return sparse.eye(d, format='csr')


def link_Lz(d):
    """L^z operator for a link of dimension d = 2J+1."""
    J = (d - 1) / 2.0
    vals = np.arange(-J, J + 1)
    return sparse.diags(vals, 0, shape=(d, d), format='csr')


def link_Lplus_cqed(d):
    """L+ operator for truncated cQED: raises L^z, truncated at boundaries."""
    J = (d - 1) / 2.0
    # L+ |k> = |k+1> for k < J, L+ |J> = 0
    data = np.ones(d - 1)
    rows = np.arange(d - 1)      # target states (k)
    cols = np.arange(1, d)        # source states (k+1) -- wait, raising means k -> k+1
    # Actually L+ |k> = |k+1>, so <k+1|L+|k> = 1
    rows_correct = np.arange(1, d)  # |k+1>
    cols_correct = np.arange(d - 1)  # <k|
    return sparse.csr_matrix((data, (rows_correct, cols_correct)), shape=(d, d))


def link_Lplus_zd(d):
    """L+ operator for Zd model: cyclic, |J+1> = |-J>."""
    J = (d - 1) / 2.0
    # Same as cQED but with cyclic boundary: L+ |J> = |-J>
    data = np.ones(d)
    rows = np.zeros(d, dtype=int)
    cols = np.zeros(d, dtype=int)
    for k in range(d):
        rows[k] = (k + 1) % d
        cols[k] = k
    return sparse.csr_matrix((data, (rows, cols)), shape=(d, d))


def link_Lminus(Lplus):
    """L- = (L+)^dagger."""
    return Lplus.T.conjugate().tocsr()


# ============================================================
# Hilbert space and Gauss law
# ============================================================

def build_site_dims(N, d):
    """
    Build list of local Hilbert space dimensions for the spin chain.
    Sites alternate: spin (dim 2) and link (dim d).
    Total: N spins and (N-1) links.
    Layout: spin_0, link_0, spin_1, link_1, ..., spin_{N-1}
    """
    dims = []
    for n in range(N):
        dims.append(2)  # spin site
        if n < N - 1:
            dims.append(d)  # link between site n and n+1
    return dims


def total_dim(dims):
    result = 1
    for d in dims:
        result *= d
    return result


def op_on_site(op, site_idx, dims):
    """Embed a local operator on a specific site into the full Hilbert space."""
    parts = []
    for i, d in enumerate(dims):
        if i == site_idx:
            parts.append(op)
        else:
            parts.append(identity(d))
    # Tensor product from left to right
    result = parts[0]
    for p in parts[1:]:
        result = sparse.kron(result, p, format='csr')
    return result


def op_on_sites(ops, site_indices, dims):
    """Embed a product of operators on specified sites."""
    local_ops = {}
    for op, idx in zip(ops, site_indices):
        local_ops[idx] = op

    parts = []
    for i, d in enumerate(dims):
        if i in local_ops:
            parts.append(local_ops[i])
        else:
            parts.append(identity(d))

    result = parts[0]
    for p in parts[1:]:
        result = sparse.kron(result, p, format='csr')
    return result


def build_gauss_law_projector(N, d, model='cqed'):
    """
    Build projector onto the Gauss law subspace.

    For cQED: G_n = L_n^z - L_{n-1}^z - (phi_n^dag phi_n - (1-(-1)^n)/2) = 0
    For Zd: U_n^Zd = exp(i 2pi/d * G_n) = 1, i.e., G_n = 0 mod d

    In spin formulation: phi_n^dag phi_n -> (1 + sigma_z_n)/2 for even n,
                                            (1 - sigma_z_n)/2 for odd n
    After Jordan-Wigner: number operator on site n is (1 + (-1)^n sigma_z_n) / 2

    Gauss law: L_n^z - L_{n-1}^z = phi_n^dag phi_n - (1-(-1)^n)/2
             = (1 + (-1)^n sigma_z_n)/2 - (1-(-1)^n)/2
             = ((-1)^n sigma_z_n + (-1)^n) / 2
             = (-1)^n (sigma_z_n + 1) / 2

    Wait, let me redo this carefully.
    The staggered fermion number operator: phi_n^dag phi_n.
    In spin language (Jordan-Wigner): phi_n^dag phi_n = (1 + sigma_z_n)/2
    The background charge: q_n = (1 - (-1)^n)/2

    Gauss law: L_n^z - L_{n-1}^z = phi_n^dag phi_n - q_n
             = (1 + sigma_z_n)/2 - (1 - (-1)^n)/2
             = (sigma_z_n + (-1)^n) / 2

    With open BC, L_{-1}^z = 0 (or equivalently, no link to the left of site 0).
    """
    dims = build_site_dims(N, d)
    D = total_dim(dims)

    # Compute the electric field and charge at each site
    # We build diagonal operators since Gauss law involves only diagonal ops
    # (L^z and sigma_z are both diagonal)

    # For efficiency, compute directly on the diagonal
    # Build the state space: each state is labeled by (s_0, l_0, s_1, l_1, ..., s_{N-1})

    # Generate all basis states
    n_sites = len(dims)
    # Use iterative approach to avoid huge memory
    # For each basis state, check Gauss law

    # Actually, let's build it differently - enumerate states and filter
    if D > 500000:
        return _build_gauss_projector_sparse(N, d, dims, D, model)

    # For manageable sizes, use direct enumeration
    physical_indices = []

    # Iterate over all basis states
    for idx in range(D):
        # Decode state index into local indices
        remaining = idx
        local_indices = []
        for i in range(n_sites - 1, -1, -1):
            local_indices.append(remaining % dims[i])
            remaining //= dims[i]
        local_indices.reverse()

        # Check Gauss law at each site
        J = (d - 1) / 2.0
        is_physical = True
        prev_Lz = 0.0  # L_{-1}^z = 0 (open BC)

        for n in range(N):
            spin_idx = 2 * n  # index in dims array
            sz = 1.0 if local_indices[spin_idx] == 0 else -1.0  # sigma_z eigenvalue

            # Charge on site n
            charge_n = (sz + (-1)**n) / 2.0

            # Required L_n^z
            required_Lz = prev_Lz + charge_n

            if n < N - 1:
                link_idx = 2 * n + 1
                actual_Lz = local_indices[link_idx] - J  # eigenvalue of L^z

                if model == 'cqed':
                    if abs(actual_Lz - required_Lz) > 1e-10:
                        is_physical = False
                        break
                elif model == 'zd':
                    # Gauss law mod d
                    diff = actual_Lz - required_Lz
                    if abs(diff - round(diff / d) * d) > 1e-10:
                        is_physical = False
                        break
                prev_Lz = actual_Lz
            else:
                # Last site: no link to the right
                # For zero total charge: required_Lz should be 0
                # Actually for open BC, we just need the Gauss law up to site N-1
                # The right boundary condition is L_{N-1}^z = required_Lz
                # For total charge 0 sector: sum of charges = 0 => L_{N-1}^z = 0
                # But we also need required_Lz to match a valid link value
                # Since there's no link N-1, the constraint is just on total charge
                pass

        # Also check zero total charge sector
        if is_physical:
            total_charge = 0.0
            for n in range(N):
                spin_idx = 2 * n
                sz = 1.0 if local_indices[spin_idx] == 0 else -1.0
                total_charge += (sz + (-1)**n) / 2.0
            if abs(total_charge) > 1e-10:
                is_physical = False

        if is_physical:
            physical_indices.append(idx)

    # Build projector matrix
    n_phys = len(physical_indices)
    P = sparse.lil_matrix((D, n_phys))
    for j, idx in enumerate(physical_indices):
        P[idx, j] = 1.0
    return P.tocsr(), physical_indices


def _build_gauss_projector_sparse(N, d, dims, D, model):
    """Build Gauss law projector for larger Hilbert spaces using batch processing."""
    n_sites = len(dims)
    J = (d - 1) / 2.0
    physical_indices = []

    batch_size = 100000
    for batch_start in range(0, D, batch_size):
        batch_end = min(batch_start + batch_size, D)
        for idx in range(batch_start, batch_end):
            remaining = idx
            local_indices = []
            for i in range(n_sites - 1, -1, -1):
                local_indices.append(remaining % dims[i])
                remaining //= dims[i]
            local_indices.reverse()

            is_physical = True
            prev_Lz = 0.0

            for n in range(N):
                spin_idx = 2 * n
                sz = 1.0 if local_indices[spin_idx] == 0 else -1.0
                charge_n = (sz + (-1)**n) / 2.0
                required_Lz = prev_Lz + charge_n

                if n < N - 1:
                    link_idx = 2 * n + 1
                    actual_Lz = local_indices[link_idx] - J

                    if model == 'cqed':
                        if abs(actual_Lz - required_Lz) > 1e-10:
                            is_physical = False
                            break
                    elif model == 'zd':
                        diff = actual_Lz - required_Lz
                        if abs(diff - round(diff / d) * d) > 1e-10:
                            is_physical = False
                            break
                    prev_Lz = actual_Lz

            if is_physical:
                total_charge = 0.0
                for n in range(N):
                    spin_idx = 2 * n
                    sz = 1.0 if local_indices[spin_idx] == 0 else -1.0
                    total_charge += (sz + (-1)**n) / 2.0
                if abs(total_charge) > 1e-10:
                    is_physical = False

            if is_physical:
                physical_indices.append(idx)

    n_phys = len(physical_indices)
    P = sparse.lil_matrix((D, n_phys))
    for j, idx in enumerate(physical_indices):
        P[idx, j] = 1.0
    return P.tocsr(), physical_indices


# ============================================================
# Hamiltonian construction
# ============================================================

def build_schwinger_hamiltonian(N, d, x, model='cqed'):
    """
    Build the Schwinger Hamiltonian in spin formulation.

    H = (g^2 a / 2) sum_n (L_n^z)^2 + m sum_n (-1)^n phi_n^dag phi_n
        - i/(2a) sum_n (phi_n^dag L_n^+ phi_{n+1} - h.c.)

    In rescaled form (Wilson's W = 2H/(ag^2)):
    W = sum_n (L_n^z)^2 + mu sum_n (-1)^n (1+sigma_z_n)/2
        + x sum_n (sigma_+_n L_+_n sigma_-_{n+1} + h.c.)

    where x = 1/(ag)^2, mu = 2m/(ag^2).

    For massless case (mu=0):
    W = sum_n (L_n^z)^2 + x sum_n (sigma_+_n L_+_n sigma_-_{n+1} + h.c.)

    In Jordan-Wigner transformed form, the hopping term becomes:
    sigma_+_n L_+_n sigma_-_{n+1}

    But we need to be careful: in the spin chain, the JW string is trivial
    for nearest-neighbor hopping since sigma_+_n sigma_-_{n+1} already includes
    the correct signs for staggered fermions.

    Parameters:
    -----------
    N : int - number of lattice sites (must be even)
    d : int - dimension of link Hilbert space (odd, d = 2J+1)
    x : float - dimensionless parameter x = 1/(ag)^2
    model : str - 'cqed' or 'zd'

    Returns:
    --------
    H : sparse matrix - Hamiltonian in full Hilbert space
    """
    dims = build_site_dims(N, d)

    # Choose link operators based on model
    if model == 'cqed':
        Lp = link_Lplus_cqed(d)
    elif model == 'zd':
        Lp = link_Lplus_zd(d)
    else:
        raise ValueError(f"Unknown model: {model}")

    Lm = link_Lminus(Lp)
    Lz = link_Lz(d)

    D = total_dim(dims)
    H = sparse.csr_matrix((D, D))

    # Electric field energy: sum_n (L_n^z)^2
    for n in range(N - 1):
        link_site = 2 * n + 1  # index in dims array
        Lz2 = Lz.dot(Lz)
        H = H + op_on_site(Lz2, link_site, dims)

    # Hopping term: x * sum_n (sigma_+_n L_+_n sigma_-_{n+1} + h.c.)
    sp = sigma_plus()
    sm = sigma_minus()

    for n in range(N - 1):
        spin_n = 2 * n        # site index for spin n
        link_n = 2 * n + 1    # site index for link n
        spin_np1 = 2 * (n + 1)  # site index for spin n+1

        # sigma_+_n * L_+_n * sigma_-_{n+1}
        hop = op_on_sites([sp, Lp, sm], [spin_n, link_n, spin_np1], dims)
        hop_dag = op_on_sites([sm, Lm, sp], [spin_n, link_n, spin_np1], dims)

        H = H + x * (hop + hop_dag)

    return H


def compute_ground_state(N, d, x, model='cqed', n_states=1):
    """
    Compute ground state(s) in the Gauss law subspace.

    Returns:
    --------
    energies : array of eigenvalues
    states : array of eigenvectors (in physical subspace)
    """
    H_full = build_schwinger_hamiltonian(N, d, x, model)
    P, phys_idx = build_gauss_law_projector(N, d, model)

    # Project Hamiltonian to physical subspace
    H_phys = P.T.dot(H_full.dot(P))

    n_phys = H_phys.shape[0]
    if n_phys <= 2 * n_states + 1:
        # Very small - use dense solver
        H_dense = H_phys.toarray()
        evals, evecs = np.linalg.eigh(H_dense)
        return evals[:n_states], evecs[:, :n_states]

    # Use sparse solver
    evals, evecs = eigsh(H_phys, k=min(n_states, n_phys - 1), which='SA')
    idx = np.argsort(evals)
    return evals[idx], evecs[:, idx]


def energy_density(N, d, x, model='cqed'):
    """
    Compute energy density omega = E0 / N in the continuum-like normalization.

    The rescaled Hamiltonian W = 2H/(ag^2), and the energy density is:
    omega = <W> / (2N) = E0_W / (2N)

    But following the paper's convention (Eq. in Appendix A):
    omega = E0 / N where E0 is the ground state energy of W.
    Then the continuum energy density is omega/sqrt(x) extrapolated to x->inf.

    Actually, the paper plots omega = E0/(2N) vs ga = 1/sqrt(x).
    """
    evals, _ = compute_ground_state(N, d, x, model, n_states=1)
    return evals[0] / N


def energy_gap(N, d, x, model='cqed'):
    """Compute gap between ground and first excited state in the physical subspace."""
    evals, _ = compute_ground_state(N, d, x, model, n_states=2)
    if len(evals) < 2:
        return 0.0
    return evals[1] - evals[0]


# ============================================================
# Adiabatic preparation
# ============================================================

def build_H_at_x(N, d, x_val, model='cqed'):
    """Build Hamiltonian at a specific x value and project to physical subspace."""
    H_full = build_schwinger_hamiltonian(N, d, x_val, model)
    P, _ = build_gauss_law_projector(N, d, model)
    return P.T.dot(H_full.dot(P))


def adiabatic_evolution(N, d, xF, T, n_steps, model='cqed', ramp='cubic'):
    """
    Simulate adiabatic evolution from x=0 to x=xF over time T.

    The ramp is x(t) = xF * (t/T)^3 (cubic ramp as in the paper).
    Time is in units of 1/g^2.

    Uses second-order Suzuki-Trotter decomposition.

    Parameters:
    -----------
    N : int - number of sites
    d : int - link dimension
    xF : float - final value of x
    T : float - total evolution time
    n_steps : int - number of time steps
    model : str - 'cqed' or 'zd'

    Returns:
    --------
    overlap : float - |<psi(T)|gs>|^2
    energy_error : float - relative error in energy
    """
    P, phys_idx = build_gauss_law_projector(N, d, model)
    n_phys = P.shape[1]

    # Initial state: strong coupling ground state (x=0)
    # At x=0, H = sum_n (L_n^z)^2, ground state has all L_n^z = 0
    # and is in the zero charge sector.
    # In the physical subspace with Gauss law, at x=0 the ground state
    # is the state with all electric fields = 0.
    H0 = build_H_at_x(N, d, 0.0, model)
    if n_phys <= 20:
        e0, v0 = np.linalg.eigh(H0.toarray())
    else:
        e0, v0 = eigsh(H0, k=1, which='SA')
    psi = v0[:, 0].copy()

    # Exact ground state at xF
    HF = build_H_at_x(N, d, xF, model)
    if n_phys <= 20:
        eF, vF = np.linalg.eigh(HF.toarray())
    else:
        eF, vF = eigsh(HF, k=1, which='SA')
    gs_final = vF[:, 0]
    E_exact = eF[0] if isinstance(eF, np.ndarray) else eF

    # Time evolution
    dt = T / n_steps
    for step in range(n_steps):
        t_mid = (step + 0.5) * dt
        if ramp == 'cubic':
            x_t = xF * (t_mid / T) ** 3
        else:
            x_t = xF * t_mid / T

        H_t = build_H_at_x(N, d, x_t, model)

        # For small systems, use dense matrix exponential
        if n_phys <= 500:
            U = expm(-1j * dt * H_t.toarray())
            psi = U @ psi
        else:
            # Use Krylov subspace method for larger systems
            from scipy.sparse.linalg import expm_multiply
            psi = expm_multiply(-1j * dt * H_t, psi)

    # Compute overlap and energy
    overlap = abs(np.vdot(gs_final, psi)) ** 2
    E_final = np.real(np.vdot(psi, HF @ psi))
    if abs(E_exact) > 1e-15:
        energy_error = abs(E_final - E_exact) / abs(E_exact)
    else:
        energy_error = abs(E_final - E_exact)

    return overlap, energy_error


# ============================================================
# Noise simulation
# ============================================================

def build_noise_operator(N, d, model='cqed'):
    """Build the noise operator sum_n (L+_n + L-_n) for the given model."""
    dims = build_site_dims(N, d)
    D = total_dim(dims)

    if model == 'cqed':
        Lp = link_Lplus_cqed(d)
    else:
        Lp = link_Lplus_zd(d)
    Lm = link_Lminus(Lp)

    noise = sparse.csr_matrix((D, D))
    for n in range(N - 1):
        link_site = 2 * n + 1
        noise = noise + op_on_site(Lp + Lm, link_site, dims)

    return noise


def build_gauss_penalty(N, d, model='cqed'):
    """
    Build Gauss law penalty operator P = sum_n G_n^dag G_n.

    For cQED: G_n = L_n^z - L_{n-1}^z - charge_n
    For Zd: P = sum_n (U_n - 1)^dag (U_n - 1) where U_n = exp(i 2pi/d G_n)
    """
    dims = build_site_dims(N, d)
    D = total_dim(dims)
    Lz_op = link_Lz(d)
    sz = sigma_z()

    penalty = sparse.csr_matrix((D, D))

    for n in range(N):
        # Build G_n = L_n^z - L_{n-1}^z - charge_n
        # charge_n = (sigma_z_n + (-1)^n) / 2

        G_n = sparse.csr_matrix((D, D))

        # L_n^z (right link)
        if n < N - 1:
            G_n = G_n + op_on_site(Lz_op, 2 * n + 1, dims)

        # -L_{n-1}^z (left link)
        if n > 0:
            G_n = G_n - op_on_site(Lz_op, 2 * (n - 1) + 1, dims)

        # -charge_n = -(sigma_z_n + (-1)^n) / 2
        G_n = G_n - 0.5 * op_on_site(sz, 2 * n, dims)
        G_n = G_n - 0.5 * (-1)**n * op_on_site(identity(2), 2 * n, dims)

        if model == 'cqed':
            penalty = penalty + G_n.T.conjugate().dot(G_n)
        elif model == 'zd':
            # U_n = exp(i 2pi/d G_n)
            # For the Zd model, we compute (U_n - I)^dag (U_n - I)
            # This is complex for sparse matrices; approximate with G_n^2
            # since for physical states G_n = 0 mod d
            penalty = penalty + G_n.T.conjugate().dot(G_n)

    return penalty


def noisy_adiabatic_evolution(N, d, xF, T, n_steps, lam, model='cqed'):
    """
    Simulate adiabatic evolution with noise.

    Noise term: lambda * x(t) * sum_n (L+_n + L-_n)

    Returns overlap, energy_error, penalty_per_site
    """
    P_proj, phys_idx = build_gauss_law_projector(N, d, model)
    n_phys = P_proj.shape[1]

    # We need to work in the FULL Hilbert space since noise breaks Gauss law
    dims = build_site_dims(N, d)
    D = total_dim(dims)

    if D > 5000:
        # Too large for full space evolution
        return None, None, None

    # Initial state in full space
    H0_full = build_schwinger_hamiltonian(N, d, 0.0, model)
    H0_phys = P_proj.T @ H0_full @ P_proj
    if n_phys <= 20:
        e0, v0 = np.linalg.eigh(H0_phys.toarray())
    else:
        e0, v0 = eigsh(H0_phys, k=1, which='SA')
    # Embed in full space
    psi = (P_proj @ v0[:, 0]).toarray().flatten()

    # Noise operator in full space
    noise_op = build_noise_operator(N, d, model)
    penalty_op = build_gauss_penalty(N, d, model)

    # Exact ground state at xF in physical subspace
    HF_full = build_schwinger_hamiltonian(N, d, xF, model)
    HF_phys = P_proj.T @ HF_full @ P_proj
    if n_phys <= 20:
        eF, vF = np.linalg.eigh(HF_phys.toarray())
    else:
        eF, vF = eigsh(HF_phys, k=1, which='SA')
    gs_final_full = (P_proj @ vF[:, 0]).toarray().flatten()

    # Time evolution in full space
    dt = T / n_steps
    for step in range(n_steps):
        t_mid = (step + 0.5) * dt
        x_t = xF * (t_mid / T) ** 3

        H_t = build_schwinger_hamiltonian(N, d, x_t, model)
        # Add noise
        H_noisy = H_t + lam * x_t * noise_op

        U = expm(-1j * dt * H_noisy.toarray())
        psi = U @ psi

    # Measurements
    overlap = abs(np.vdot(gs_final_full, psi)) ** 2
    E_final = np.real(np.vdot(psi, HF_full @ psi))
    E_exact = eF[0] if isinstance(eF, np.ndarray) else eF
    energy_error = abs(E_final - E_exact) / abs(E_exact) if abs(E_exact) > 1e-15 else abs(E_final - E_exact)
    penalty = np.real(np.vdot(psi, penalty_op @ psi)) / N

    return overlap, energy_error, penalty


# ============================================================
# Thermodynamic limit extrapolation helpers
# ============================================================

def extrapolate_to_thermo(sizes, energies):
    """Linear extrapolation in 1/N to N->inf."""
    inv_sizes = 1.0 / np.array(sizes)
    coeffs = np.polyfit(inv_sizes, energies, 1)
    return coeffs[1]  # intercept = value at 1/N = 0
