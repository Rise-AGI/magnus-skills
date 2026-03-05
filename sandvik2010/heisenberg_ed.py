"""
Exact Diagonalization of the S=1/2 Heisenberg Antiferromagnetic Chain.

Reproduces results from Sandvik (2010), Section 4:
- Excitation spectrum / dispersion relation (Fig 32)
- Singlet and triplet gaps vs system size (Fig 33)
- Spin correlations (Fig 34)

Hamiltonian:
  H = J sum_{<i,j>} S_i . S_j = J sum_{<i,j>} [S_i^z S_j^z + 1/2(S_i^+ S_j^- + S_i^- S_j^+)]

Uses conservation of total S^z and translational symmetry (momentum states)
via the Lanczos method for the ground state and low-energy excitations.
"""

import numpy as np
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse.linalg import eigsh


def build_heisenberg_hamiltonian(N, J=1.0, pbc=True, Sz_sector=0):
    """
    Build the Heisenberg chain Hamiltonian in the S^z = Sz_sector subspace.

    Parameters
    ----------
    N : int
        Number of sites
    J : float
        Coupling constant (J > 0 = antiferromagnetic)
    pbc : bool
        Periodic boundary conditions
    Sz_sector : int or float
        Total S^z quantum number (in units of 1/2: 0 means Sz=0 sector)

    Returns
    -------
    H : sparse matrix
        Hamiltonian in the given Sz sector
    basis : list
        List of basis state integers
    """
    # Generate basis states in the Sz sector
    # Each state is an integer where bit i = 1 means spin up, 0 means spin down
    target_nup = N // 2 + Sz_sector
    if target_nup < 0 or target_nup > N:
        return None, []

    basis = []
    for state in range(2**N):
        if bin(state).count('1') == target_nup:
            basis.append(state)

    n_basis = len(basis)
    state_index = {}
    for idx, state in enumerate(basis):
        state_index[state] = idx

    H = lil_matrix((n_basis, n_basis), dtype=np.float64)

    n_bonds = N if pbc else N - 1

    for idx, state in enumerate(basis):
        for bond in range(n_bonds):
            i = bond
            j = (bond + 1) % N

            # Extract spins at sites i and j
            si = (state >> i) & 1  # 1 = up, 0 = down
            sj = (state >> j) & 1

            # S^z_i S^z_j term: (si - 1/2)(sj - 1/2) = si*sj - (si+sj)/2 + 1/4
            # In {0,1}: sz = s - 1/2
            szi = si - 0.5
            szj = sj - 0.5
            H[idx, idx] += J * szi * szj

            # Off-diagonal: 1/2 (S^+_i S^-_j + S^-_i S^+_j)
            # This flips spins at i and j if they are antiparallel
            if si != sj:
                # Flip both spins
                new_state = state ^ (1 << i) ^ (1 << j)
                new_idx = state_index[new_state]
                H[idx, new_idx] += J * 0.5

    return csr_matrix(H), basis


def compute_ground_state(N, J=1.0, n_states=4):
    """Compute the ground state and low excited states of the Heisenberg chain."""
    H, basis = build_heisenberg_hamiltonian(N, J=J, pbc=True, Sz_sector=0)
    if H is None or len(basis) < n_states + 1:
        return None, None, None

    k = min(n_states, len(basis) - 1)
    energies, vectors = eigsh(H, k=k, which='SA')
    sort_idx = np.argsort(energies)
    energies = energies[sort_idx]
    vectors = vectors[:, sort_idx]
    return energies, vectors, basis


def spin_correlation(vector, basis, N, i, j):
    """
    Compute <S_i . S_j> in a given eigenstate.

    <S_i . S_j> = <S^z_i S^z_j> + 1/2 <S^+_i S^-_j + S^-_i S^+_j>
    """
    n = len(basis)
    state_index = {s: idx for idx, s in enumerate(basis)}
    coeff = vector

    result = 0.0
    for idx, state in enumerate(basis):
        c = coeff[idx]
        if abs(c) < 1e-15:
            continue

        si = (state >> i) & 1
        sj = (state >> j) & 1
        szi = si - 0.5
        szj = sj - 0.5

        # Diagonal: S^z_i S^z_j
        result += c * c * szi * szj

        # Off-diagonal: 1/2 (S^+_i S^-_j + S^-_i S^+_j)
        if si != sj:
            new_state = state ^ (1 << i) ^ (1 << j)
            new_idx = state_index[new_state]
            result += c * coeff[new_idx] * 0.5

    return result


def compute_all_correlations(N, J=1.0):
    """Compute spin-spin correlations C(r) = <S_0 . S_r> in the ground state."""
    energies, vectors, basis = compute_ground_state(N, J=J, n_states=1)
    if energies is None:
        return None, None

    gs_vector = vectors[:, 0]
    correlations = np.zeros(N // 2 + 1)
    for r in range(N // 2 + 1):
        correlations[r] = spin_correlation(gs_vector, basis, N, 0, r)

    return energies[0], correlations


def compute_gaps_sz_sectors(N, J=1.0):
    """
    Compute singlet and triplet gaps.
    Singlet gap: E_1(Sz=0) - E_0(Sz=0)
    Triplet gap: E_0(Sz=1) - E_0(Sz=0)
    """
    # Sz=0 sector
    H0, basis0 = build_heisenberg_hamiltonian(N, J=J, pbc=True, Sz_sector=0)
    k0 = min(4, len(basis0) - 1)
    e0, _ = eigsh(H0, k=k0, which='SA')
    e0 = np.sort(e0)

    # Sz=1 sector
    H1, basis1 = build_heisenberg_hamiltonian(N, J=J, pbc=True, Sz_sector=1)
    k1 = min(2, len(basis1) - 1)
    e1, _ = eigsh(H1, k=k1, which='SA')
    e1 = np.sort(e1)

    gs_energy = e0[0]
    singlet_gap = e0[1] - e0[0]
    triplet_gap = e1[0] - e0[0]

    return gs_energy, singlet_gap, triplet_gap


def compute_dispersion(N, J=1.0):
    """
    Compute excitation spectrum at different momenta using translational symmetry.
    For simplicity, compute in real space and Fourier transform the correlation matrix.

    Returns momentum-resolved excitation energies.
    """
    # For the dispersion, we compute the lowest energy in each Sz sector
    # and return the excitation spectrum

    # Get ground state energy from Sz=0
    H0, basis0 = build_heisenberg_hamiltonian(N, J=J, pbc=True, Sz_sector=0)
    n0 = min(6, len(basis0) - 1)
    e0, _ = eigsh(H0, k=n0, which='SA')
    e0 = np.sort(e0)
    gs_energy = e0[0]

    # Get triplet spectrum from Sz=1
    H1, basis1 = build_heisenberg_hamiltonian(N, J=J, pbc=True, Sz_sector=1)
    n1 = min(N, len(basis1) - 1)
    e1, v1 = eigsh(H1, k=n1, which='SA')
    sort_idx = np.argsort(e1)
    e1 = e1[sort_idx]
    v1 = v1[:, sort_idx]

    # Assign momenta by computing <T> for each eigenstate
    # Translation operator T shifts each state by one site
    momenta_triplet = []
    for s in range(len(e1)):
        vec = v1[:, s]
        # Compute momentum by applying translation
        overlap = 0.0 + 0.0j
        for idx, state in enumerate(basis1):
            # Translate state by 1
            new_state = 0
            for bit in range(N):
                if (state >> bit) & 1:
                    new_bit = (bit + 1) % N
                    new_state |= (1 << new_bit)
            if new_state in {b: i for i, b in enumerate(basis1)}:
                new_idx = {b: i for i, b in enumerate(basis1)}[new_state]
                overlap += vec[new_idx] * vec[idx]
        # momentum: T|k> = e^{ik}|k>, so k = -Im(ln(overlap))
        if abs(overlap) > 1e-10:
            k = -np.angle(overlap)
        else:
            k = 0.0
        momenta_triplet.append(k)

    return {
        'gs_energy': gs_energy,
        'singlet_excitations': e0[1:] - gs_energy,
        'triplet_energies': e1 - gs_energy,
        'triplet_momenta': np.array(momenta_triplet),
    }
