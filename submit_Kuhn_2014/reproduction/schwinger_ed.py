"""
Exact diagonalization for the Schwinger model with finite-dimensional link variables.

Implements two models from Kuhn, Cirac, Banuls (2014):
  1. Truncated compact QED (cQED) - U(1) gauge symmetry preserved
  2. Zd model - discrete Zd gauge symmetry

Uses direct construction in the Gauss law constrained subspace for efficiency.
The key insight: once the fermion configuration is fixed, the link electric
fields are completely determined by the Gauss law (for cQED) or determined
modulo d (for Zd).
"""

import numpy as np
from scipy.sparse import csr_matrix, eye as speye, kron as skron
from scipy.sparse.linalg import eigsh
from scipy.linalg import expm


def get_physical_states(N, d, model='cqed'):
    """
    Enumerate physical states in the zero-charge sector satisfying the Gauss law.

    For each fermion configuration {s_0, ..., s_{N-1}} (s_n in {0,1} = spin up/down),
    the Gauss law determines the link electric fields.

    Returns list of (fermion_config, link_config) tuples.
    fermion_config: array of 0/1 (spin up = 0, spin down = 1)
    link_config: array of L^z values on each link
    """
    J = (d - 1) / 2.0
    states = []

    for bits in range(2**N):
        # Decode fermion configuration
        spins = np.zeros(N, dtype=int)
        for n in range(N):
            spins[n] = (bits >> (N - 1 - n)) & 1  # 0 = up (sz=+1), 1 = down (sz=-1)

        # sigma_z values: +1 for spin up, -1 for spin down
        sz = 1 - 2 * spins  # [+1 or -1]

        # Check total charge = 0
        total_charge = sum((sz[n] + (-1)**n) / 2.0 for n in range(N))
        if abs(total_charge) > 1e-10:
            continue

        # Compute required link electric fields from Gauss law
        # L_n^z = L_{n-1}^z + charge_n
        # charge_n = (sz_n + (-1)^n) / 2
        # L_{-1}^z = 0 (open BC)
        links = np.zeros(N - 1)
        prev_Lz = 0.0
        valid = True

        for n in range(N):
            charge_n = (sz[n] + (-1)**n) / 2.0
            curr_Lz = prev_Lz + charge_n

            if n < N - 1:
                if model == 'cqed':
                    # Link value must be in [-J, J]
                    if curr_Lz < -J - 1e-10 or curr_Lz > J + 1e-10:
                        valid = False
                        break
                    links[n] = curr_Lz
                elif model == 'zd':
                    # Link value is modulo d, mapped to [-J, J]
                    # The Zd model: Gauss law is mod d
                    # For the unique physical state, we take curr_Lz mod d
                    # mapped to the range [-J, J]
                    lz_mod = curr_Lz % d
                    if lz_mod > J + 0.5:
                        lz_mod -= d
                    links[n] = lz_mod
                prev_Lz = links[n] if n < N - 1 else curr_Lz

        if valid:
            states.append((spins.copy(), links.copy()))

    return states


def state_to_index(spins, links, N, d):
    """Convert (spins, links) to an index in the full Hilbert space."""
    J = (d - 1) / 2.0
    # Layout: spin_0, link_0, spin_1, link_1, ..., spin_{N-1}
    idx = 0
    stride = 1
    # Build from right to left
    for n in range(N - 1, -1, -1):
        # spin site
        idx += spins[n] * stride
        stride *= 2
        # link site (if exists)
        if n > 0:
            link_val = int(round(links[n - 1] + J))  # convert from L^z value to index
            idx += link_val * stride
            stride *= d
    return idx


def build_hamiltonian_physical(N, d, x, model='cqed', mu=0.0):
    """
    Build the Schwinger Hamiltonian directly in the physical (Gauss law) subspace.

    H_W = sum_n (L_n^z)^2 + mu * sum_n (-1)^n (1+sz_n)/2
          + x * sum_n (sigma+_n L+_n sigma-_{n+1} + h.c.)

    In the physical subspace, the first two terms are diagonal.
    The hopping term connects different physical states.

    Parameters:
    -----------
    N : int - number of sites (must be even)
    d : int - link Hilbert space dimension (odd)
    x : float - dimensionless parameter x = 1/(ag)^2
    model : str - 'cqed' or 'zd'
    mu : float - dimensionless mass parameter (0 for massless)
    """
    J = (d - 1) / 2.0
    states = get_physical_states(N, d, model)
    n_phys = len(states)

    if n_phys == 0:
        return csr_matrix((0, 0)), states

    # Build state lookup: map fermion config to index
    state_map = {}
    for i, (spins, links) in enumerate(states):
        key = tuple(spins)
        if key not in state_map:
            state_map[key] = []
        state_map[key].append(i)

    rows, cols, vals = [], [], []

    for i, (spins, links) in enumerate(states):
        # Diagonal: electric field energy
        E_field = np.sum(links**2)

        # Diagonal: mass term
        mass_term = 0.0
        if mu != 0:
            for n in range(N):
                sz = 1 - 2 * spins[n]
                mass_term += mu * (-1)**n * (1 + sz) / 2.0

        rows.append(i)
        cols.append(i)
        vals.append(E_field + mass_term)

        # Off-diagonal: hopping terms
        # sigma+_n L+_n sigma-_{n+1}: flip spin n from down to up, spin n+1 from up to down
        # sigma+ |down> = |up>, sigma- |up> = |down>
        # In our convention: spin=0 is up (sz=+1), spin=1 is down (sz=-1)
        # sigma+ |1> = |0>, sigma- |0> = |1>
        for n in range(N - 1):
            # Forward hop: sigma+_n L+_n sigma-_{n+1}
            # Need: spin n is down (1), spin n+1 is up (0)
            if spins[n] == 1 and spins[n + 1] == 0:
                new_spins = spins.copy()
                new_spins[n] = 0
                new_spins[n + 1] = 1

                # The new state must also satisfy Gauss law
                # Compute new links
                new_links = np.zeros(N - 1)
                prev_Lz = 0.0
                new_valid = True
                for m in range(N - 1):
                    new_sz = 1 - 2 * new_spins[m]
                    charge_m = (new_sz + (-1)**m) / 2.0
                    new_Lz = prev_Lz + charge_m

                    if model == 'cqed':
                        if new_Lz < -J - 1e-10 or new_Lz > J + 1e-10:
                            new_valid = False
                            break
                        new_links[m] = new_Lz
                    elif model == 'zd':
                        lz_mod = new_Lz % d
                        if lz_mod > J + 0.5:
                            lz_mod -= d
                        new_links[m] = lz_mod
                    prev_Lz = new_links[m]

                if new_valid:
                    new_key = tuple(new_spins)
                    if new_key in state_map:
                        for j in state_map[new_key]:
                            if np.allclose(states[j][1], new_links):
                                # Matrix element for L+ between old and new link states
                                # For cQED: <l+1|L+|l> = 1 (for truncated, 0 if l=J)
                                # For Zd: <(l+1)%d|L+|l> = 1
                                old_l = links[n]
                                if model == 'cqed':
                                    if old_l < J - 1e-10:
                                        mel = 1.0
                                    else:
                                        mel = 0.0
                                elif model == 'zd':
                                    mel = 1.0  # Always valid for cyclic

                                if abs(mel) > 1e-10:
                                    rows.append(j)
                                    cols.append(i)
                                    vals.append(x * mel)
                                break

            # Backward hop (hermitian conjugate): sigma-_n L-_n sigma+_{n+1}
            # Need: spin n is up (0), spin n+1 is down (1)
            if spins[n] == 0 and spins[n + 1] == 1:
                new_spins = spins.copy()
                new_spins[n] = 1
                new_spins[n + 1] = 0

                new_links = np.zeros(N - 1)
                prev_Lz = 0.0
                new_valid = True
                for m in range(N - 1):
                    new_sz = 1 - 2 * new_spins[m]
                    charge_m = (new_sz + (-1)**m) / 2.0
                    new_Lz = prev_Lz + charge_m

                    if model == 'cqed':
                        if new_Lz < -J - 1e-10 or new_Lz > J + 1e-10:
                            new_valid = False
                            break
                        new_links[m] = new_Lz
                    elif model == 'zd':
                        lz_mod = new_Lz % d
                        if lz_mod > J + 0.5:
                            lz_mod -= d
                        new_links[m] = lz_mod
                    prev_Lz = new_links[m]

                if new_valid:
                    new_key = tuple(new_spins)
                    if new_key in state_map:
                        for j in state_map[new_key]:
                            if np.allclose(states[j][1], new_links):
                                old_l = links[n]
                                if model == 'cqed':
                                    if old_l > -J + 1e-10:
                                        mel = 1.0
                                    else:
                                        mel = 0.0
                                elif model == 'zd':
                                    mel = 1.0

                                if abs(mel) > 1e-10:
                                    rows.append(j)
                                    cols.append(i)
                                    vals.append(x * mel)
                                break

    H = csr_matrix((vals, (rows, cols)), shape=(n_phys, n_phys))
    # Ensure Hermitian
    H = (H + H.T.conjugate()) / 2.0
    return H, states


def compute_ground_state(N, d, x, model='cqed', n_states=1):
    """Compute ground state(s) in the Gauss law subspace."""
    H, states = build_hamiltonian_physical(N, d, x, model)
    n_phys = H.shape[0]

    if n_phys == 0:
        return np.array([]), np.array([])

    if n_phys <= 2 * n_states + 1:
        H_dense = H.toarray()
        evals, evecs = np.linalg.eigh(H_dense)
        return evals[:n_states], evecs[:, :n_states]

    evals, evecs = eigsh(H, k=min(n_states, n_phys - 1), which='SA')
    idx = np.argsort(evals)
    return evals[idx], evecs[:, idx]


def energy_density(N, d, x, model='cqed'):
    """Compute energy density omega = E0/N."""
    evals, _ = compute_ground_state(N, d, x, model, n_states=1)
    if len(evals) == 0:
        return np.nan
    return evals[0] / N


def energy_gap(N, d, x, model='cqed'):
    """Compute gap between ground and first excited state."""
    evals, _ = compute_ground_state(N, d, x, model, n_states=2)
    if len(evals) < 2:
        return 0.0
    return evals[1] - evals[0]


def adiabatic_evolution(N, d, xF, T, n_steps, model='cqed'):
    """
    Simulate adiabatic evolution from x=0 to x=xF over time T using cubic ramp.

    Returns (overlap, energy_error).
    """
    # Initial state: ground state at x=0
    H0, states0 = build_hamiltonian_physical(N, d, 0.0, model)
    n_phys = H0.shape[0]
    if n_phys <= 1:
        return 1.0, 0.0

    if n_phys <= 20:
        e0, v0 = np.linalg.eigh(H0.toarray())
    else:
        e0, v0 = eigsh(H0, k=1, which='SA')
    psi = v0[:, 0].astype(complex).copy()

    # Target: ground state at xF
    HF, _ = build_hamiltonian_physical(N, d, xF, model)
    if n_phys <= 20:
        eF, vF = np.linalg.eigh(HF.toarray())
    else:
        eF, vF = eigsh(HF, k=1, which='SA')
    gs_final = vF[:, 0]
    E_exact = eF[0]

    # Time evolution with Suzuki-Trotter
    dt = T / n_steps
    for step in range(n_steps):
        t_mid = (step + 0.5) * dt
        x_t = xF * (t_mid / T) ** 3
        H_t, _ = build_hamiltonian_physical(N, d, x_t, model)

        if n_phys <= 300:
            U = expm(-1j * dt * H_t.toarray())
            psi = U @ psi
        else:
            from scipy.sparse.linalg import expm_multiply
            psi = expm_multiply(-1j * dt * H_t, psi)

    overlap = abs(np.vdot(gs_final, psi)) ** 2
    E_final = np.real(np.vdot(psi, HF @ psi))
    energy_error = abs(E_final - E_exact) / abs(E_exact) if abs(E_exact) > 1e-15 else abs(E_final - E_exact)

    return overlap, energy_error


def build_noise_hamiltonian(N, d, x_val, lam, model='cqed'):
    """
    Build the noisy Hamiltonian in the FULL Hilbert space (not just physical).
    Needed because noise breaks Gauss law.

    Returns H_noisy as dense matrix (only feasible for small systems).
    """
    J = (d - 1) / 2.0
    # Build site dims
    dims = []
    for n in range(N):
        dims.append(2)
        if n < N - 1:
            dims.append(d)

    D = 1
    for dd in dims:
        D *= dd

    if D > 5000:
        return None, None, D

    # Build H as dense matrix in full space
    H = np.zeros((D, D), dtype=complex)

    # We need to map each basis state to its local quantum numbers
    # and build the Hamiltonian matrix element by element

    n_sites = len(dims)

    def decode_state(idx):
        local = []
        remaining = idx
        for i in range(n_sites - 1, -1, -1):
            local.append(remaining % dims[i])
            remaining //= dims[i]
        local.reverse()
        return local

    def encode_state(local):
        idx = 0
        stride = 1
        for i in range(n_sites - 1, -1, -1):
            idx += local[i] * stride
            stride *= dims[i]
        return idx

    # Electric field energy (diagonal)
    for idx in range(D):
        local = decode_state(idx)
        E = 0.0
        for n in range(N - 1):
            link_idx = 2 * n + 1
            lz_val = local[link_idx] - J
            E += lz_val ** 2
        H[idx, idx] += E

    # Hopping terms
    for n in range(N - 1):
        spin_n_idx = 2 * n
        link_n_idx = 2 * n + 1
        spin_np1_idx = 2 * (n + 1)

        for idx in range(D):
            local = decode_state(idx)

            # sigma+_n L+_n sigma-_{n+1}: spin_n: 1->0, link_n: l->l+1, spin_{n+1}: 0->1
            if local[spin_n_idx] == 1 and local[spin_np1_idx] == 0:
                new_local = local.copy()
                new_local[spin_n_idx] = 0
                new_local[spin_np1_idx] = 1
                lz_old = local[link_n_idx]

                if model == 'cqed':
                    if lz_old < d - 1:
                        new_local[link_n_idx] = lz_old + 1
                        jdx = encode_state(new_local)
                        H[jdx, idx] += x_val
                elif model == 'zd':
                    new_local[link_n_idx] = (lz_old + 1) % d
                    jdx = encode_state(new_local)
                    H[jdx, idx] += x_val

            # sigma-_n L-_n sigma+_{n+1}: spin_n: 0->1, link_n: l->l-1, spin_{n+1}: 1->0
            if local[spin_n_idx] == 0 and local[spin_np1_idx] == 1:
                new_local = local.copy()
                new_local[spin_n_idx] = 1
                new_local[spin_np1_idx] = 0
                lz_old = local[link_n_idx]

                if model == 'cqed':
                    if lz_old > 0:
                        new_local[link_n_idx] = lz_old - 1
                        jdx = encode_state(new_local)
                        H[jdx, idx] += x_val
                elif model == 'zd':
                    new_local[link_n_idx] = (lz_old - 1) % d
                    jdx = encode_state(new_local)
                    H[jdx, idx] += x_val

    # Noise term: lam * x * sum_n (L+_n + L-_n)
    if lam > 0:
        for n in range(N - 1):
            link_n_idx = 2 * n + 1
            for idx in range(D):
                local = decode_state(idx)
                lz_old = local[link_n_idx]

                # L+ term
                if model == 'cqed':
                    if lz_old < d - 1:
                        new_local = local.copy()
                        new_local[link_n_idx] = lz_old + 1
                        jdx = encode_state(new_local)
                        H[jdx, idx] += lam * x_val
                elif model == 'zd':
                    new_local = local.copy()
                    new_local[link_n_idx] = (lz_old + 1) % d
                    jdx = encode_state(new_local)
                    H[jdx, idx] += lam * x_val

                # L- term
                if model == 'cqed':
                    if lz_old > 0:
                        new_local = local.copy()
                        new_local[link_n_idx] = lz_old - 1
                        jdx = encode_state(new_local)
                        H[jdx, idx] += lam * x_val
                elif model == 'zd':
                    new_local = local.copy()
                    new_local[link_n_idx] = (lz_old - 1) % d
                    jdx = encode_state(new_local)
                    H[jdx, idx] += lam * x_val

    # Ensure Hermitian
    H = (H + H.T.conj()) / 2.0
    return H, dims, D


def build_gauss_penalty_full(N, d, dims, D, model='cqed'):
    """Build Gauss law penalty operator in full Hilbert space."""
    J = (d - 1) / 2.0
    n_sites = len(dims)

    def decode_state(idx):
        local = []
        remaining = idx
        for i in range(n_sites - 1, -1, -1):
            local.append(remaining % dims[i])
            remaining //= dims[i]
        local.reverse()
        return local

    penalty = np.zeros(D)
    for idx in range(D):
        local = decode_state(idx)
        p = 0.0
        prev_Lz = 0.0
        for n in range(N):
            spin_idx = 2 * n
            sz = 1.0 if local[spin_idx] == 0 else -1.0
            charge_n = (sz + (-1)**n) / 2.0

            if n < N - 1:
                link_idx = 2 * n + 1
                actual_Lz = local[link_idx] - J
                Gn = actual_Lz - prev_Lz - charge_n
                p += Gn ** 2
                prev_Lz = actual_Lz
        penalty[idx] = p

    return np.diag(penalty)


def noisy_adiabatic_evolution(N, d, xF, T, n_steps, lam, model='cqed'):
    """
    Simulate adiabatic evolution with noise in the full Hilbert space.

    Returns (overlap, energy_error, penalty_per_site).
    """
    # Check feasibility
    test_H, test_dims, D = build_noise_hamiltonian(N, d, 0.0, 0.0, model)
    if test_H is None:
        return None, None, None

    # Initial state: ground state of H(x=0) in physical subspace, embedded in full space
    # At x=0, H = sum (Lz)^2, ground state has all Lz=0
    # Use physical subspace to find it
    H_phys, phys_states = build_hamiltonian_physical(N, d, 0.0, model)
    n_phys = H_phys.shape[0]
    if n_phys <= 20:
        e0, v0 = np.linalg.eigh(H_phys.toarray())
    else:
        e0, v0 = eigsh(H_phys, k=1, which='SA')

    # Embed physical ground state in full space
    psi = np.zeros(D, dtype=complex)
    for i, coeff in enumerate(v0[:, 0]):
        spins, links = phys_states[i]
        full_idx = state_to_full_index(spins, links, N, d, test_dims)
        psi[full_idx] = coeff

    # Exact ground state at xF in physical subspace
    HF_phys, _ = build_hamiltonian_physical(N, d, xF, model)
    if n_phys <= 20:
        eF, vF = np.linalg.eigh(HF_phys.toarray())
    else:
        eF, vF = eigsh(HF_phys, k=1, which='SA')

    gs_final_full = np.zeros(D, dtype=complex)
    for i, coeff in enumerate(vF[:, 0]):
        spins, links = phys_states[i]
        full_idx = state_to_full_index(spins, links, N, d, test_dims)
        gs_final_full[full_idx] = coeff

    E_exact = eF[0]

    # Time evolution
    dt = T / n_steps
    for step in range(n_steps):
        t_mid = (step + 0.5) * dt
        x_t = xF * (t_mid / T) ** 3
        H_t, _, _ = build_noise_hamiltonian(N, d, x_t, lam, model)
        U = expm(-1j * dt * H_t)
        psi = U @ psi

    # Measurements
    overlap = abs(np.vdot(gs_final_full, psi)) ** 2

    HF_full, _, _ = build_noise_hamiltonian(N, d, xF, 0.0, model)
    E_final = np.real(np.vdot(psi, HF_full @ psi))
    energy_error = abs(E_final - E_exact) / abs(E_exact) if abs(E_exact) > 1e-15 else abs(E_final - E_exact)

    P_op = build_gauss_penalty_full(N, d, test_dims, D, model)
    penalty = np.real(np.vdot(psi, P_op @ psi)) / N

    return overlap, energy_error, penalty


def state_to_full_index(spins, links, N, d, dims):
    """Convert (spins, links) to index in full Hilbert space."""
    J = (d - 1) / 2.0
    local = []
    for n in range(N):
        local.append(spins[n])
        if n < N - 1:
            local.append(int(round(links[n] + J)))

    idx = 0
    stride = 1
    for i in range(len(local) - 1, -1, -1):
        idx += local[i] * stride
        stride *= dims[i]
    return idx


def extrapolate_to_thermo(sizes, energies):
    """Linear extrapolation in 1/N to N->inf."""
    inv_sizes = 1.0 / np.array(sizes)
    coeffs = np.polyfit(inv_sizes, energies, 1)
    return coeffs[1]
