"""
Core module for reproducing figures from:
  "Detailed Analysis of the Three Quark Potential in SU(3) Lattice QCD"
  T.T. Takahashi, H. Suganuma, Y. Nemoto, H. Matsufuru
  Phys. Rev. D 65, 114509 (2002)

Implements:
  - L_min (minimal flux-tube length) via Fermat point geometry
  - Lattice Coulomb potential V^LC(n)
  - Y-ansatz and Delta-ansatz fit functions
  - Data tables extracted from the paper
"""
import numpy as np
from scipy import optimize, integrate


# =============================================================================
# Physical constants and lattice parameters (Table 1)
# =============================================================================
LATTICE_PARAMS = {
    5.7: {"a_fm": 0.19, "lattice": (12, 12, 12, 24), "N3Q": 16, "Nconf": 210,
           "alpha_smear": 2.3, "Nsmr": 12},
    5.8: {"a_fm": 0.14, "lattice": (16, 16, 16, 32), "N3Q": 139, "Nconf": 200,
           "alpha_smear": 2.3, "Nsmr": 22},
    6.0: {"a_fm": 0.10, "lattice": (16, 16, 16, 32), "N3Q": 155, "Nconf": 150,
           "alpha_smear": 2.3, "Nsmr": 42},
}

# Table 2: Best-fit parameters in lattice units
FIT_PARAMS = {
    5.7: {
        "3QY":    {"sigma": 0.1524, "A": 0.1331, "C": 0.9182},
        "3QY_LC": {"sigma": 0.1556, "A": 0.1185, "C": 0.8876},
        "QQ":     {"sigma": 0.1629, "A": 0.2793, "C": 0.6203},
        "QQ_LC":  {"sigma": 0.1603, "A": 0.2627, "C": 0.6271},
    },
    5.8: {
        "3QY":    {"sigma": 0.1027, "A": 0.1230, "C": 0.9085},
        "3QY_LC": {"sigma": 0.1031, "A": 0.1141, "C": 0.8999},
        "QQ":     {"sigma": 0.1079, "A": 0.2607, "C": 0.6115},
        "QQ_LC":  {"sigma": 0.1018, "A": 0.2795, "C": 0.6596},
    },
    6.0: {
        "3QY":    {"sigma": 0.0460, "A": 0.1366, "C": 0.9599},
        "3QY_LC": {"sigma": 0.0467, "A": 0.1256, "C": 0.9467},
        "QQ":     {"sigma": 0.0506, "A": 0.2768, "C": 0.6374},
        "QQ_LC":  {"sigma": 0.0500, "A": 0.2557, "C": 0.6373},
    },
}

# Delta-ansatz fit parameters (Table 12)
DELTA_FIT_PARAMS = {
    5.7: {"sigma": 0.0876, "A": 0.1431, "C": 0.9378, "chi2": 10.1},
    5.8: {"sigma": 0.0586, "A": 0.1260, "C": 0.9127, "chi2": 13.7},
    6.0: {"sigma": 0.0256, "A": 0.1356, "C": 0.9421, "chi2": 7.85},
}

# Generalized Y-ansatz fit parameters (Table 13, best R_core)
GEN_Y_FIT_PARAMS = {
    5.8: {"R_core_fm": 0.08, "sigma": 0.1021, "A": 0.1220, "C": 0.9074, "chi2": 4.07},
    6.0: {"R_core_fm": 0.08, "sigma": 0.0458, "A": 0.1358, "C": 0.9588, "chi2": 2.17},
}


# =============================================================================
# Table 3: Lattice QCD data at beta=5.7 (Type I: quarks at (i,0,0),(0,j,0),(0,0,k))
# =============================================================================
DATA_BETA_57 = [
    # (i, j, k), V_3Q_latt, error, C_bar, V_3Q_fit_Y
    ((0, 1, 1), 0.8457, 0.0038, 0.9338, 0.8524),
    ((0, 1, 2), 1.0973, 0.0043, 0.9295, 1.1025),
    ((0, 1, 3), 1.2929, 0.0041, 0.8987, 1.2929),
    ((0, 2, 2), 1.3158, 0.0044, 0.9151, 1.3270),
    ((0, 2, 3), 1.5040, 0.0063, 0.9041, 1.5076),
    ((0, 3, 3), 1.6756, 0.0043, 0.8718, 1.6815),
    ((1, 1, 1), 1.0238, 0.0040, 0.9345, 1.0092),
    ((1, 1, 2), 1.2185, 0.0062, 0.9067, 1.2151),
    ((1, 1, 3), 1.4161, 0.0049, 0.9297, 1.3964),
    ((1, 2, 2), 1.3866, 0.0048, 0.9012, 1.3895),
    ((1, 2, 3), 1.5594, 0.0063, 0.8880, 1.5588),
    ((1, 3, 3), 1.7145, 0.0043, 0.8553, 1.7202),
    ((2, 2, 2), 1.5234, 0.0037, 0.8925, 1.5238),
    ((2, 2, 3), 1.6750, 0.0118, 0.8627, 1.6763),
    ((2, 3, 3), 1.8239, 0.0056, 0.8443, 1.8175),
    ((3, 3, 3), 1.9607, 0.0093, 0.8197, 1.9442),
]

# =============================================================================
# Table 4-7: Lattice QCD data at beta=5.8 (Type I, selected subset)
# =============================================================================
DATA_BETA_58 = [
    ((0, 1, 1), 0.7697, 0.0012, 0.9554),
    ((0, 1, 2), 0.9639, 0.0028, 0.9269),
    ((0, 1, 3), 1.1112, 0.0060, 0.9274),
    ((0, 1, 4), 1.2337, 0.0119, 0.9106),
    ((0, 2, 2), 1.1370, 0.0016, 0.9342),
    ((0, 2, 3), 1.2659, 0.0021, 0.9145),
    ((0, 3, 3), 1.3925, 0.0091, 0.9168),
    ((0, 3, 4), 1.5005, 0.0041, 0.8862),
    ((0, 3, 5), 1.6130, 0.0025, 0.8810),
    ((0, 3, 6), 1.7171, 0.0032, 0.8581),
    ((0, 4, 4), 1.6077, 0.0072, 0.8655),
    ((0, 4, 5), 1.7163, 0.0032, 0.8581),
    ((0, 4, 6), 1.8262, 0.0040, 0.8482),
    ((0, 5, 5), 1.8193, 0.0044, 0.8372),
    ((0, 5, 6), 1.9282, 0.0047, 0.8265),
    ((0, 6, 6), 2.0322, 0.0062, 0.8083),
    ((0, 7, 7), 2.2461, 0.0101, 0.7813),
    ((0, 8, 8), 2.4191, 0.0177, 0.6949),
    ((1, 1, 1), 0.9140, 0.0032, 0.9424),
    ((1, 1, 2), 1.0647, 0.0042, 0.9290),
    ((1, 1, 3), 1.1914, 0.0086, 0.8917),
    ((1, 2, 2), 1.1865, 0.0033, 0.9186),
    ((1, 2, 3), 1.3126, 0.0124, 0.9411),
    ((1, 3, 3), 1.4175, 0.0033, 0.9020),
    ((2, 2, 2), 1.2771, 0.0073, 0.9000),
    ((2, 2, 3), 1.3783, 0.0080, 0.8755),
    ((2, 3, 3), 1.4739, 0.0120, 0.8636),
    ((2, 3, 4), 1.5831, 0.0022, 0.8718),
    ((3, 3, 3), 1.5566, 0.0072, 0.8434),
    ((3, 3, 4), 1.6474, 0.0066, 0.8125),
    ((3, 3, 5), 1.7641, 0.0037, 0.8353),
    ((3, 3, 6), 1.8685, 0.0044, 0.8196),
    ((3, 3, 7), 1.9696, 0.0052, 0.7965),
    ((3, 3, 8), 2.0753, 0.0065, 0.7812),
    ((4, 4, 4), 1.8377, 0.0049, 0.8044),
    ((4, 4, 5), 1.9371, 0.0055, 0.7900),
    ((4, 4, 6), 2.0367, 0.0061, 0.7703),
    ((4, 5, 5), 2.0278, 0.0068, 0.7638),
    ((4, 5, 6), 2.1301, 0.0069, 0.7503),
    ((5, 5, 5), 2.1192, 0.0087, 0.7417),
    ((5, 5, 6), 2.2247, 0.0090, 0.7334),
    ((6, 6, 6), 2.4166, 0.0223, 0.6943),
]


# =============================================================================
# Geometry: Fermat point and L_min
# =============================================================================

def fermat_point(p1, p2, p3):
    """
    Find the Fermat point of triangle with vertices p1, p2, p3.
    The Fermat point minimizes the sum of distances to the three vertices.
    If any angle >= 120 degrees, the Fermat point is at that vertex.
    """
    p1, p2, p3 = np.array(p1, float), np.array(p2, float), np.array(p3, float)
    vertices = [p1, p2, p3]

    # Check for degenerate cases (two or more coincident points)
    sides = [np.linalg.norm(p2 - p3), np.linalg.norm(p1 - p3), np.linalg.norm(p1 - p2)]
    for i in range(3):
        if sides[i] < 1e-12:
            return vertices[i]

    # Check if any angle >= 120 degrees
    for i in range(3):
        j, k = (i + 1) % 3, (i + 2) % 3
        v1 = vertices[j] - vertices[i]
        v2 = vertices[k] - vertices[i]
        cos_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
        cos_angle = np.clip(cos_angle, -1, 1)
        if cos_angle <= np.cos(2 * np.pi / 3):  # angle >= 120 deg
            return vertices[i].copy()

    # Use optimization to find the Fermat point
    def total_dist(xy):
        p = np.array([xy[0], xy[1], 0.0]) if len(p1) == 3 else np.array(xy)
        return sum(np.linalg.norm(p - v) for v in vertices)

    if len(p1) == 3:
        x0 = (p1 + p2 + p3) / 3
        res = optimize.minimize(lambda xy: sum(np.linalg.norm(np.array([xy[0], xy[1], xy[2]]) - v)
                                               for v in vertices), x0, method='Nelder-Mead')
        return res.x
    else:
        x0 = (p1 + p2 + p3) / 3
        res = optimize.minimize(total_dist, x0[:2], method='Nelder-Mead')
        return np.array([res.x[0], res.x[1]])


def compute_lmin(i, j, k):
    """
    Compute L_min for a Type-I 3Q configuration:
    quarks at (i,0,0), (0,j,0), (0,0,k) in lattice units.
    L_min = sum of distances from Fermat point to each quark.
    """
    p1 = np.array([float(i), 0.0, 0.0])
    p2 = np.array([0.0, float(j), 0.0])
    p3 = np.array([0.0, 0.0, float(k)])

    fp = fermat_point(p1, p2, p3)
    lmin = (np.linalg.norm(fp - p1) + np.linalg.norm(fp - p2) +
            np.linalg.norm(fp - p3))
    return lmin


def compute_lmin_formula(a, b, c):
    """
    Compute L_min using the closed-form formula (Eq. 2 in paper).
    a, b, c are side lengths of the 3Q triangle.
    """
    # Check if any angle > 120 degrees
    sides = sorted([a, b, c])
    # If longest side forms angle > 120 deg at opposite vertex
    # Using law of cosines: c^2 = a^2 + b^2 - 2ab cos(C)
    # For the angle at the vertex opposite the longest side:
    max_side = sides[2]
    other1, other2 = sides[0], sides[1]
    cos_angle = (other1**2 + other2**2 - max_side**2) / (2 * other1 * other2 + 1e-30)
    if cos_angle <= np.cos(2 * np.pi / 3):
        # Angle >= 120 deg: L_min = sum of two shorter sides
        return other1 + other2

    # Formula from paper (Eq. 2)
    s = a + b + c
    area_sq = s * (-a + b + c) * (a - b + c) * (a + b - c)
    if area_sq < 0:
        area_sq = 0
    lmin = np.sqrt(0.5 * (a**2 + b**2 + c**2) + np.sqrt(3) / 2 * np.sqrt(area_sq))
    return lmin


def quark_distances_type1(i, j, k):
    """Return the three inter-quark distances for Type-I config."""
    p1 = np.array([float(i), 0.0, 0.0])
    p2 = np.array([0.0, float(j), 0.0])
    p3 = np.array([0.0, 0.0, float(k)])
    return (np.linalg.norm(p1 - p2), np.linalg.norm(p1 - p3), np.linalg.norm(p2 - p3))


# =============================================================================
# Lattice Coulomb Potential
# =============================================================================

def lattice_coulomb_potential(n_vec, N_grid=128):
    """
    Compute the lattice Coulomb potential V^LC(n) via numerical integration.
    V^LC(n) = pi * integral d^3q/(2pi)^3 exp(-i p.n) / sum_i sin^2(p_i/2)
    where the integral is over [-pi, pi]^3.
    Note: lattice spacing a=1 (lattice units).
    """
    n1, n2, n3 = n_vec

    # Use numerical integration via midpoint rule on a grid
    dp = 2 * np.pi / N_grid
    p = np.linspace(-np.pi + dp / 2, np.pi - dp / 2, N_grid)
    px, py, pz = np.meshgrid(p, p, p, indexing='ij')

    denom = np.sin(px / 2)**2 + np.sin(py / 2)**2 + np.sin(pz / 2)**2
    denom = np.maximum(denom, 1e-30)  # avoid division by zero

    # Exponential factor
    phase = px * n1 + py * n2 + pz * n3
    integrand = np.cos(phase) / denom  # imaginary part vanishes by symmetry

    # The integral: pi * (dp/(2*pi))^3 * sum
    vlc = np.pi * (dp / (2 * np.pi))**3 * np.sum(integrand)
    return vlc


def lattice_coulomb_table(max_n=10, N_grid=64):
    """
    Compute V^LC for various lattice vectors and return as dict.
    """
    vlc = {}
    # Compute for unique |n| values
    computed = set()
    for n1 in range(max_n + 1):
        for n2 in range(n1 + 1):
            for n3 in range(n2 + 1):
                key = (n1, n2, n3)
                if key not in computed:
                    vlc[key] = lattice_coulomb_potential(key, N_grid)
                    computed.add(key)
    return vlc


def get_vlc(n_vec, vlc_table):
    """Look up V^LC for a given vector, using symmetry."""
    key = tuple(sorted(np.abs(n_vec), reverse=True))
    return vlc_table.get(key, 1.0 / max(np.linalg.norm(n_vec), 1e-10))


# =============================================================================
# Potential ansatze
# =============================================================================

def v3q_y_ansatz(config, A, sigma, C):
    """Y-ansatz: V_3Q = -A * sum_{i<j} 1/|r_i - r_j| + sigma * L_min + C"""
    i, j, k = config
    d12, d13, d23 = quark_distances_type1(i, j, k)
    coulomb = 0
    for d in [d12, d13, d23]:
        if d > 1e-10:
            coulomb += 1.0 / d
    lmin = compute_lmin(i, j, k)
    return -A * coulomb + sigma * lmin + C


def v3q_delta_ansatz(config, A, sigma_delta, C):
    """Delta-ansatz: V_3Q = -A * sum 1/|r_i-r_j| + sigma_delta * sum |r_i-r_j| + C"""
    i, j, k = config
    d12, d13, d23 = quark_distances_type1(i, j, k)
    coulomb = 0
    perimeter = 0
    for d in [d12, d13, d23]:
        if d > 1e-10:
            coulomb += 1.0 / d
        perimeter += d
    return -A * coulomb + sigma_delta * perimeter + C


def vqq_potential(r, A, sigma, C):
    """Q-Qbar potential: V = -A/r + sigma*r + C"""
    return -A / np.maximum(r, 1e-10) + sigma * r + C


# =============================================================================
# Fit routines
# =============================================================================

def fit_y_ansatz(data):
    """
    Fit 3Q potential data with Y-ansatz.
    data: list of ((i,j,k), V_latt, error, ...)
    Returns: (A, sigma, C), cov
    """
    configs = [d[0] for d in data]
    v_latt = np.array([d[1] for d in data])
    errors = np.array([d[2] for d in data])

    # Precompute Lmin and Coulomb sums
    lmins = []
    coulombs = []
    for cfg in configs:
        i, j, k = cfg
        lmins.append(compute_lmin(i, j, k))
        d12, d13, d23 = quark_distances_type1(i, j, k)
        coul = sum(1.0 / d for d in [d12, d13, d23] if d > 1e-10)
        coulombs.append(coul)
    lmins = np.array(lmins)
    coulombs = np.array(coulombs)

    def model(x, A, sigma, C):
        return -A * coulombs + sigma * lmins + C

    popt, pcov = optimize.curve_fit(model, None, v_latt, p0=[0.13, 0.15, 0.9],
                                     sigma=errors, absolute_sigma=True)
    return popt, pcov


def fit_delta_ansatz(data):
    """Fit with Delta-ansatz."""
    configs = [d[0] for d in data]
    v_latt = np.array([d[1] for d in data])
    errors = np.array([d[2] for d in data])

    perimeters = []
    coulombs = []
    for cfg in configs:
        i, j, k = cfg
        d12, d13, d23 = quark_distances_type1(i, j, k)
        coul = sum(1.0 / d for d in [d12, d13, d23] if d > 1e-10)
        peri = sum([d12, d13, d23])
        coulombs.append(coul)
        perimeters.append(peri)
    coulombs = np.array(coulombs)
    perimeters = np.array(perimeters)

    def model(x, A, sigma_d, C):
        return -A * coulombs + sigma_d * perimeters + C

    popt, pcov = optimize.curve_fit(model, None, v_latt, p0=[0.13, 0.06, 0.9],
                                     sigma=errors, absolute_sigma=True)
    return popt, pcov


def chi2_per_dof(v_latt, v_fit, errors, n_params=3):
    """Compute chi^2/NDF."""
    residuals = (v_latt - v_fit) / errors
    return np.sum(residuals**2) / max(len(v_latt) - n_params, 1)
