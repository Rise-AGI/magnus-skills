"""
Core tight-binding module for graphene band structure and density of states.

Implements the tight-binding model with up to 3rd nearest neighbor
hopping and overlap integrals, following:
  R. Kundu, "Tight binding parameters for graphene" (2024)
"""

import numpy as np

# Lattice constants
A0 = 1.42  # nearest neighbour C-C distance in Angstroms
A = np.sqrt(3) * A0  # lattice constant ~ 2.46 Angstroms

# High symmetry points in reciprocal space (units of 1/Angstrom)
GAMMA = np.array([0.0, 0.0])
K_POINT = np.array([2 * np.pi / (np.sqrt(3) * A), 2 * np.pi / (3 * A)])
M_POINT = np.array([2 * np.pi / (np.sqrt(3) * A), 0.0])


# ---------- Structure factors ----------

def g_func(kx, ky):
    """g(k) = |f(k)|^2 = 1 + 4 cos^2(ky*a/2) + 4 cos(sqrt3*kx*a/2) cos(ky*a/2)"""
    return (1.0
            + 4.0 * np.cos(ky * A / 2.0)**2
            + 4.0 * np.cos(np.sqrt(3) * kx * A / 2.0) * np.cos(ky * A / 2.0))


def u_func(kx, ky):
    """u(k) = 2 cos(ky*a) + 4 cos(kx*a*sqrt3/2) cos(ky*a/2)"""
    return (2.0 * np.cos(ky * A)
            + 4.0 * np.cos(kx * A * np.sqrt(3) / 2.0) * np.cos(ky * A / 2.0))


def g2_func(kx, ky):
    """g(2k) — g evaluated at doubled arguments."""
    return (1.0
            + 4.0 * np.cos(ky * A)**2
            + 4.0 * np.cos(np.sqrt(3) * kx * A) * np.cos(ky * A))


def t_func(kx, ky):
    """t(k) cross-term for 3rd NN."""
    c1 = np.cos(kx * A * np.sqrt(3))
    c2 = np.cos(ky * A)
    c3 = np.cos(ky * A / 2.0)
    c4 = np.cos(kx * A * np.sqrt(3) / 2.0)
    return 2.0 * c1 + 4.0 * c2 + 4.0 * c3 * c4 + 8.0 * c2 * c3 * c4


# ---------- Dispersion relations ----------

def band_nearest(kx, ky, E2p=0.0, gamma0=-2.74, s0=0.065):
    """Nearest neighbour dispersion (Eq. 5)."""
    gk = g_func(kx, ky)
    sqrtg = np.sqrt(np.maximum(gk, 0.0))
    num_plus = (E2p - s0 * gamma0 * gk) + (gamma0 - s0 * E2p) * sqrtg
    num_minus = (E2p - s0 * gamma0 * gk) - (gamma0 - s0 * E2p) * sqrtg
    denom = 1.0 - s0**2 * gk
    # Avoid division by zero
    denom = np.where(np.abs(denom) < 1e-15, 1e-15, denom)
    return num_plus / denom, num_minus / denom


def band_nearest_no_overlap(kx, ky, E2p=0.0, gamma0=-2.74):
    """Nearest neighbour without overlap (s0=0)."""
    gk = g_func(kx, ky)
    sqrtg = np.sqrt(np.maximum(gk, 0.0))
    return E2p + gamma0 * sqrtg, E2p - gamma0 * sqrtg


def band_2nn(kx, ky, E2p=-0.21, gamma0=-2.74, gamma1=-0.07, s0=0.065, s1=0.002):
    """Second nearest neighbour dispersion (Eq. 6)."""
    gk = g_func(kx, ky)
    sqrtg = np.sqrt(np.maximum(gk, 0.0))
    uk = u_func(kx, ky)
    # Note the mp sign: valence band uses -, conduction uses +
    num_val = E2p + gamma1 * uk - gamma0 * sqrtg
    den_val = 1.0 + s1 * uk - s0 * sqrtg
    num_cond = E2p + gamma1 * uk + gamma0 * sqrtg
    den_cond = 1.0 + s1 * uk + s0 * sqrtg
    den_val = np.where(np.abs(den_val) < 1e-15, 1e-15, den_val)
    den_cond = np.where(np.abs(den_cond) < 1e-15, 1e-15, den_cond)
    return num_cond / den_cond, num_val / den_val


def band_3nn(kx, ky, E2p=-0.21, gamma0=-2.74, gamma1=-0.07, gamma2=-0.015,
             s0=0.065, s1=0.002, s2=0.001):
    """Third nearest neighbour dispersion (general formula Eq. 4 with Eq. 7)."""
    gk = g_func(kx, ky)
    uk = u_func(kx, ky)
    g2k = g2_func(kx, ky)
    tk = t_func(kx, ky)

    HAA = E2p + gamma1 * uk
    SAA = 1.0 + s1 * uk

    E0 = SAA * HAA
    E1 = (2.0 * s0 * gamma0 * gk
           + (s0 * gamma2 + gamma0 * s2) * tk
           + 2.0 * s2 * gamma2 * g2k)
    E2 = (HAA**2
           - (gamma0**2 * gk + gamma0 * gamma2 * tk + gamma2**2 * g2k))
    E3 = (SAA**2
           - (s0**2 * gk + s0 * s2 * tk + s2**2 * g2k))

    # Quadratic: E3*E^2 - (2E0 - E1)*E + E2 = 0
    discriminant = (2.0 * E0 - E1)**2 - 4.0 * E3 * E2
    # Clamp to avoid sqrt of negative due to numerics
    discriminant = np.maximum(discriminant, 0.0)
    sqrt_disc = np.sqrt(discriminant)

    denom = 2.0 * E3
    denom = np.where(np.abs(denom) < 1e-15, 1e-15, denom)
    Eplus = ((2.0 * E0 - E1) + sqrt_disc) / denom
    Eminus = ((2.0 * E0 - E1) - sqrt_disc) / denom
    return Eplus, Eminus


# ---------- k-path generation ----------

def kpath_GKMGamma(npts=300):
    """Generate k-path along Gamma-K-M-Gamma and return (kx, ky, distances)."""
    segments = [
        (GAMMA, K_POINT),
        (K_POINT, M_POINT),
        (M_POINT, GAMMA),
    ]
    kx_all, ky_all, dist_all = [], [], []
    cumulative = 0.0
    for i, (start, end) in enumerate(segments):
        n = npts
        t_arr = np.linspace(0, 1, n, endpoint=(i == len(segments) - 1))
        kx_seg = start[0] + t_arr * (end[0] - start[0])
        ky_seg = start[1] + t_arr * (end[1] - start[1])
        # Compute distances along the path
        dk = np.sqrt((kx_seg[1:] - kx_seg[:-1])**2 + (ky_seg[1:] - ky_seg[:-1])**2)
        d_seg = np.zeros(len(t_arr))
        d_seg[1:] = np.cumsum(dk)
        d_seg += cumulative
        cumulative = d_seg[-1]

        kx_all.append(kx_seg)
        ky_all.append(ky_seg)
        dist_all.append(d_seg)

    return np.concatenate(kx_all), np.concatenate(ky_all), np.concatenate(dist_all)


def get_tick_positions(npts=300):
    """Return distances for Gamma, K, M, Gamma tick marks."""
    kx, ky, dist = kpath_GKMGamma(npts)
    # Find approximate positions
    dGK = np.linalg.norm(K_POINT - GAMMA)
    dKM = np.linalg.norm(M_POINT - K_POINT)
    dMG = np.linalg.norm(GAMMA - M_POINT)
    return [0.0, dGK, dGK + dKM, dGK + dKM + dMG]


# ---------- Density of states ----------

def compute_dos(band_func, params, n_kpoints=500, n_bins=500, energy_range=(-10, 15)):
    """Compute DOS by sampling over the full Brillouin zone.

    Uses a uniform grid in kx, ky over the first Brillouin zone
    and histograms the energies.
    """
    # Sample over rectangular region covering the BZ
    kx_max = 2 * np.pi / (np.sqrt(3) * A) * 1.2
    ky_max = 2 * np.pi / A * 0.7
    kx = np.linspace(-kx_max, kx_max, n_kpoints)
    ky = np.linspace(-ky_max, ky_max, n_kpoints)
    KX, KY = np.meshgrid(kx, ky)
    KX_flat = KX.ravel()
    KY_flat = KY.ravel()

    Eplus, Eminus = band_func(KX_flat, KY_flat, **params)

    energies = np.concatenate([Eplus, Eminus])
    energies = energies[np.isfinite(energies)]

    bins = np.linspace(energy_range[0], energy_range[1], n_bins + 1)
    hist, bin_edges = np.histogram(energies, bins=bins)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    return bin_centers, hist.astype(float)
