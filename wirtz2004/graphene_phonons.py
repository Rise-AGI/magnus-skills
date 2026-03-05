"""
Core module for graphene phonon dispersion using force-constant models.

Implements:
- 4th-nearest-neighbor force constant (4NNFC) model
- Valence-force-field (VFF) model of Aizawa et al.
- Dynamical matrix construction and diagonalization
- Phonon dispersion along high-symmetry paths
- Vibrational density of states (vDOS)

Reference: Wirtz & Rubio, cond-mat/0404637 (2004)
"""

import numpy as np


# ---------- Constants ----------

AMU = 1.66053906660e-27  # kg
M_C = 12.011 * AMU       # carbon mass

# Unit conversions
CM_INV_TO_THZ = 0.02998   # 1 cm^-1 = 0.02998 THz
THZ_TO_CM_INV = 1.0 / CM_INV_TO_THZ
DYN_CM_TO_N_M = 1e-5      # 1 dyn/cm = 1e-5 N/m = 1e-5 kg/s^2


# ---------- Lattice geometry ----------

def lattice_vectors(a):
    """Return lattice vectors a1, a2 for graphene (in Angstrom)."""
    a1 = a * np.array([1.0, 0.0])
    a2 = a * np.array([0.5, np.sqrt(3) / 2])
    return a1, a2


def reciprocal_vectors(a):
    """Return reciprocal lattice vectors b1, b2."""
    b1 = (2 * np.pi / a) * np.array([1.0, -1.0 / np.sqrt(3)])
    b2 = (2 * np.pi / a) * np.array([0.0, 2.0 / np.sqrt(3)])
    return b1, b2


def high_symmetry_path(a, n_points=200):
    """Generate q-points along Gamma-M-K-Gamma path.

    Returns:
        q_points: (N, 2) array of q vectors
        q_dist: (N,) array of cumulative distance along path
        tick_positions: list of distances at Gamma, M, K, Gamma
        tick_labels: list of labels
    """
    b1, b2 = reciprocal_vectors(a)

    G = np.array([0.0, 0.0])
    M = b1 / 2
    K = (2 * b1 + b2) / 3

    segments = [(G, M), (M, K), (K, G)]
    q_points = []
    q_dist = []
    tick_positions = [0.0]
    tick_labels = [r'$\Gamma$', 'M', 'K', r'$\Gamma$']

    cumulative = 0.0
    for i, (start, end) in enumerate(segments):
        n = n_points
        for j in range(n):
            if i > 0 and j == 0:
                continue  # avoid duplicate at segment boundaries
            t = j / (n - 1)
            q = start + t * (end - start)
            q_points.append(q)
            if len(q_dist) > 0:
                cumulative += np.linalg.norm(q - prev_q)
            q_dist.append(cumulative)
            prev_q = q
        tick_positions.append(cumulative)

    return np.array(q_points), np.array(q_dist), tick_positions, tick_labels


def neighbor_positions(a):
    """Compute neighbor positions for 4 nearest-neighbor shells.

    Returns dict with keys 1..4, each value is a list of:
        (delta, sublattice) where sublattice is 'AB' or 'AA'
    delta is the position vector from atom A to the neighbor.
    """
    a1, a2 = lattice_vectors(a)
    # B atom in unit cell (fractional 1/3, 1/3)
    tau = (a1 + a2) / 3  # = a*(1/2, sqrt(3)/6)

    neighbors = {}

    # 1st NN: A->B, 3 atoms at distance d = a/sqrt(3)
    nn1 = [
        (tau, 'AB'),
        (tau - a1, 'AB'),
        (tau - a2, 'AB'),
    ]
    neighbors[1] = nn1

    # 2nd NN: A->A, 6 atoms at distance a
    nn2 = [
        (a1, 'AA'), (-a1, 'AA'),
        (a2, 'AA'), (-a2, 'AA'),
        (a1 - a2, 'AA'), (a2 - a1, 'AA'),
    ]
    neighbors[2] = nn2

    # 3rd NN: A->B, 3 atoms at distance 2d
    nn3 = [
        (-tau + a1 + a2, 'AB'),   # = a(1, 2*sqrt(3)/3) - tau
        (tau + a1 - a2 - tau + tau, 'AB'),  # need to compute carefully
    ]
    # Let me recompute: 3rd NN B atoms from A at origin
    # B atoms at positions: tau + n*a1 + m*a2
    # We need distance = 2d = 2a/sqrt(3)
    # tau - a1 + a2 direction: tau - a1 + a2 = (a1+a2)/3 - a1 + a2 = (-2a1+4a2)/3
    # |(-2a1+4a2)/3| = ... Let me compute numerically
    d = a / np.sqrt(3)
    target_3 = 2 * d

    all_b_neighbors = []
    for n in range(-3, 4):
        for m in range(-3, 4):
            pos = tau + n * a1 + m * a2
            dist = np.linalg.norm(pos)
            if abs(dist - target_3) < 0.01:
                all_b_neighbors.append(pos)

    nn3 = [(pos, 'AB') for pos in all_b_neighbors]
    neighbors[3] = nn3

    # 4th NN: A->A, 6 atoms at distance a*sqrt(3)
    target_4 = a * np.sqrt(3)
    all_a_neighbors = []
    for n in range(-3, 4):
        for m in range(-3, 4):
            if n == 0 and m == 0:
                continue
            pos = n * a1 + m * a2
            dist = np.linalg.norm(pos)
            if abs(dist - target_4) < 0.01:
                all_a_neighbors.append(pos)

    nn4 = [(pos, 'AA') for pos in all_a_neighbors]
    neighbors[4] = nn4

    return neighbors


# ---------- Force constant rotation ----------

def rotate_fc_matrix(C_bond, theta):
    """Rotate a 3x3 force constant matrix from bond frame to global frame.

    C_bond is in the frame where x = bond direction.
    theta is the angle of the bond direction w.r.t. global x-axis.
    """
    c, s = np.cos(theta), np.sin(theta)
    R = np.array([
        [c, -s, 0],
        [s,  c, 0],
        [0,  0, 1]
    ])
    return R @ C_bond @ R.T


def fc_matrix_4nnfc(phi_l, phi_ti, phi_to, xi=0.0):
    """Build force constant matrix for 4NNFC model (Eq. 4/5).

    Args:
        phi_l: longitudinal force constant
        phi_ti: transverse in-plane force constant
        phi_to: transverse out-of-plane force constant
        xi: off-diagonal coupling (0 for standard 4NNFC)
    """
    return np.array([
        [phi_l,  xi,    0],
        [-xi,    phi_ti, 0],
        [0,      0,     phi_to]
    ])


# ---------- Dynamical matrix ----------

def build_dynamical_matrix(q, a, fc_params, model='4nnfc'):
    """Build the 6x6 dynamical matrix at wavevector q.

    Args:
        q: (2,) wavevector
        a: lattice constant
        fc_params: dict of force constant parameters
        model: '4nnfc' or 'vff'

    Returns:
        D: (6, 6) complex dynamical matrix (not yet divided by mass)
    """
    nbrs = neighbor_positions(a)

    if model == '4nnfc':
        fc_matrices = {}
        for n in range(1, 5):
            phi_l = fc_params[f'phi_l_{n}']
            phi_ti = fc_params[f'phi_ti_{n}']
            phi_to = fc_params[f'phi_to_{n}']
            xi = fc_params.get(f'xi_{n}', 0.0)
            fc_matrices[n] = fc_matrix_4nnfc(phi_l, phi_ti, phi_to, xi)
    elif model == 'vff':
        fc_matrices = build_vff_fc_matrices(a, fc_params)

    D = np.zeros((6, 6), dtype=complex)

    # The C_n in the paper (Eq. 5) are spring constant matrices (positive).
    # The interatomic force constant Phi_AB = -C_n (opposite atoms)
    # The on-site Phi_AA(0) = +sum(C_n) over all neighbors (acoustic sum rule)
    #
    # D_ss'(q) = sum_R Phi_ss'(R) exp(iq.R)
    # D_AB(q) = sum_R (-C_rotated) exp(iq.R)
    # D_AA(q) = sum_{R!=0} (-C_AA_rotated) exp(iq.R) + Phi_AA(0)
    #         = sum_all C_rotated - sum_AA C_rotated*exp(iq.R)

    D_AB = np.zeros((3, 3), dtype=complex)
    D_AA = np.zeros((3, 3), dtype=complex)
    sum_all_C = np.zeros((3, 3), dtype=float)  # for on-site term

    for n in range(1, 5):
        C_bond = fc_matrices[n]
        for delta, sublattice in nbrs[n]:
            theta = np.arctan2(delta[1], delta[0])
            C_global = rotate_fc_matrix(C_bond, theta)
            phase = np.exp(1j * np.dot(q, delta))

            sum_all_C += C_global.real  # accumulate for on-site

            if sublattice == 'AB':
                D_AB -= C_global * phase  # Phi_AB = -C
            else:  # AA
                D_AA -= C_global * phase  # Phi_AA(R!=0) = -C

    # Add on-site: Phi_AA(0) = +sum_all C
    D_AA += sum_all_C

    # Assemble 6x6 matrix
    D[0:3, 0:3] = D_AA
    D[0:3, 3:6] = D_AB
    D[3:6, 0:3] = D_AB.conj().T
    D[3:6, 3:6] = D_AA.copy()  # BB = AA by symmetry

    return D


def build_vff_fc_matrices(a, params):
    """Build force constant matrices for VFF model (Eqs. 7-10).

    Args:
        a: lattice constant
        params: dict with keys alpha1, alpha2, gamma1, gamma2, delta

    Returns:
        dict of 3x3 matrices for neighbors 1-4
    """
    d = a / np.sqrt(3)  # CC bond length in Angstrom

    alpha1 = params['alpha1']  # 10^4 dyn/cm
    alpha2 = params['alpha2']
    gamma1 = params['gamma1']  # 10^-13 erg
    gamma2 = params['gamma2']
    delta_p = params['delta']

    # Unit conversion: gamma (10^-13 erg) / length^2 (Angstrom^2)
    # = gamma*1e-13 / (length*1e-8)^2 = gamma/length^2 * 1e3 dyn/cm
    # In units of 10^4 dyn/cm: multiply by 0.1
    conv = 0.1  # conversion factor

    # In-plane bending: gamma1/d^2
    g1 = gamma1 * conv / (d ** 2)

    # Out-of-plane: Eq. 8 explicitly uses a^2 for gamma2 term.
    # For consistency, all out-of-plane terms (gamma2, delta) use a^2.
    g2a = gamma2 * conv / (a ** 2)
    da = delta_p * conv / (a ** 2)

    # Eqs. 7-10 with corrected length scales:
    C1 = np.array([
        [alpha1, 0, 0],
        [0, 4.5 * g1, 0],
        [0, 0, 6 * g2a * 3]  # 18*gamma2/(a^2) = 6*g2a*3
    ])
    # C1[2,2] = 18*gamma2*conv/a^2 using lattice constant for z-modes
    C1[2, 2] = 18 * g2a

    sqrt3 = np.sqrt(3)
    C2 = np.array([
        [alpha2 + 0.75 * g1, 0.75 * sqrt3 * g1, 0],
        [-0.75 * sqrt3 * g1, -2.25 * g1, 0],
        [0, 0, -3 * g2a + da]
    ])

    C3 = np.array([
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 2 * da]
    ])

    C4 = np.array([
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, -da]
    ])

    return {1: C1, 2: C2, 3: C3, 4: C4}


# ---------- Phonon dispersion ----------

def phonon_frequencies(q, a, fc_params, model='4nnfc'):
    """Compute phonon frequencies at wavevector q.

    Args:
        q: wavevector (2D)
        a: lattice constant (Angstrom)
        fc_params: force constant parameters
        model: '4nnfc' or 'vff'

    Returns:
        freqs: 6 frequencies in cm^-1, sorted ascending
    """
    D = build_dynamical_matrix(q, a, fc_params, model)

    # Convert to SI units:
    # D is in units of 10^4 dyn/cm
    # 1 dyn/cm = 10^-5 N / 10^-2 m = 10^-3 N/m
    # So 10^4 dyn/cm = 10^4 * 10^-3 N/m = 10 N/m
    D_SI = D * 10.0  # N/m = kg/s^2
    D_over_M = D_SI / M_C  # s^-2

    eigenvalues = np.linalg.eigvalsh(D_over_M.real)  # D should be Hermitian

    # Handle numerical noise (small negative eigenvalues)
    freqs_hz = np.sqrt(np.abs(eigenvalues)) * np.sign(eigenvalues)

    # Convert rad/s to cm^-1
    c_cgs = 2.99792458e10  # cm/s
    freqs_cm = freqs_hz / (2 * np.pi * c_cgs)

    return np.sort(freqs_cm)


def compute_dispersion(a, fc_params, model='4nnfc', n_points=200):
    """Compute phonon dispersion along Gamma-M-K-Gamma.

    Returns:
        q_dist: distances along path
        frequencies: (n_q, 6) array of frequencies in cm^-1
        tick_positions, tick_labels: for plotting
    """
    q_points, q_dist, ticks, labels = high_symmetry_path(a, n_points)

    frequencies = np.zeros((len(q_points), 6))
    for i, q in enumerate(q_points):
        frequencies[i] = phonon_frequencies(q, a, fc_params, model)

    return q_dist, frequencies, ticks, labels


# ---------- Vibrational density of states ----------

def compute_vdos(a, fc_params, model='4nnfc', n_grid=50, n_bins=500,
                 freq_max=1700):
    """Compute vibrational density of states on a uniform q-grid.

    Args:
        n_grid: number of q-points along each reciprocal axis
        n_bins: number of frequency bins
        freq_max: maximum frequency in cm^-1

    Returns:
        freq_bins: bin centers in cm^-1
        dos: density of states (normalized)
    """
    b1, b2 = reciprocal_vectors(a)

    freq_bins = np.linspace(0, freq_max, n_bins)
    dos = np.zeros(n_bins)
    bin_width = freq_bins[1] - freq_bins[0]

    for i in range(n_grid):
        for j in range(n_grid):
            q = (i / n_grid) * b1 + (j / n_grid) * b2
            freqs = phonon_frequencies(q, a, fc_params, model)
            for f in freqs:
                if 0 < f < freq_max:
                    idx = int(f / bin_width)
                    if 0 <= idx < n_bins:
                        dos[idx] += 1

    # Normalize
    total = np.sum(dos) * bin_width
    if total > 0:
        dos /= total

    # Gaussian smoothing
    sigma = 10.0  # cm^-1
    n_sigma = int(3 * sigma / bin_width)
    kernel = np.exp(-0.5 * (np.arange(-n_sigma, n_sigma + 1) * bin_width / sigma) ** 2)
    kernel /= kernel.sum()
    dos_smooth = np.convolve(dos, kernel, mode='same')

    return freq_bins, dos_smooth


# ---------- Parameter sets from Table 3 ----------

def params_4nnfc_jishi1993():
    """4NNFC parameters from Jishi et al. (Ref. [5])."""
    return {
        'phi_l_1': 36.50, 'phi_l_2': 8.80, 'phi_l_3': 3.00, 'phi_l_4': -1.92,
        'phi_ti_1': 24.50, 'phi_ti_2': -3.23, 'phi_ti_3': -5.25, 'phi_ti_4': 2.29,
        'phi_to_1': 9.82, 'phi_to_2': -0.40, 'phi_to_3': 0.15, 'phi_to_4': -0.58,
    }


def params_4nnfc_gruneis2002():
    """4NNFC parameters from Gruneis et al. (Ref. [3])."""
    return {
        'phi_l_1': 40.37, 'phi_l_2': 2.76, 'phi_l_3': 0.05, 'phi_l_4': 1.31,
        'phi_ti_1': 25.18, 'phi_ti_2': 2.22, 'phi_ti_3': -8.99, 'phi_ti_4': 0.22,
        'phi_to_1': 9.40, 'phi_to_2': -0.08, 'phi_to_3': -0.06, 'phi_to_4': -0.63,
    }


def params_4nnfc_gga_diag():
    """4NNFC diagonal fit to GGA (this work, Table 3 col 3)."""
    return {
        'phi_l_1': 39.87, 'phi_l_2': 7.29, 'phi_l_3': -2.64, 'phi_l_4': 0.10,
        'phi_ti_1': 17.28, 'phi_ti_2': -4.61, 'phi_ti_3': 3.31, 'phi_ti_4': 0.79,
        'phi_to_1': 9.89, 'phi_to_2': -0.82, 'phi_to_3': 0.58, 'phi_to_4': -0.52,
    }


def params_4nnfc_gga_offdiag():
    """4NNFC + off-diagonal fit to GGA (this work, Table 3 col 4)."""
    return {
        'phi_l_1': 40.98, 'phi_l_2': 7.42, 'phi_l_3': -3.32, 'phi_l_4': 0.65,
        'phi_ti_1': 14.50, 'phi_ti_2': -4.08, 'phi_ti_3': 5.01, 'phi_ti_4': 0.55,
        'phi_to_1': 9.89, 'phi_to_2': -0.82, 'phi_to_3': 0.58, 'phi_to_4': -0.52,
        'xi_2': -0.91,
    }


def params_vff_aizawa1990():
    """VFF parameters from Aizawa et al. (Ref. [11])."""
    return {
        'alpha1': 36.4, 'alpha2': 6.2,
        'gamma1': 83.0, 'gamma2': 33.8, 'delta': 31.7,
    }


def params_vff_siebentritt1997():
    """VFF parameters from Siebentritt et al. (Ref. [12])."""
    return {
        'alpha1': 34.4, 'alpha2': 6.2,
        'gamma1': 93.0, 'gamma2': 30.8, 'delta': 41.7,
    }


def params_vff_gga():
    """VFF fit to GGA (this work, Table 3)."""
    return {
        'alpha1': 39.9, 'alpha2': 5.7,
        'gamma1': 60.8, 'gamma2': 32.8, 'delta': 34.6,
    }


# Lattice constants
A_LDA = 2.449  # Angstrom
A_GGA = 2.457  # Angstrom
