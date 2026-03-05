"""
Core module for reproducing Ferrenberg, Xu & Landau,
Phys. Rev. E 97, 043301 (2018).

Implements:
- 3D Ising model Wolff cluster algorithm
- Histogram reweighting
- Thermodynamic observables (energy, magnetization, cumulants)
- Finite-size scaling fitting utilities

All figures import from this module to avoid code duplication.
"""

import numpy as np
from scipy.optimize import curve_fit


# ---------------------------------------------------------------------------
# Published data from Tables II-IX of the paper
# ---------------------------------------------------------------------------

# Table II: nu vs L_min with 1, 2, 3 correction exponents
TABLE_II = {
    "L_min": [16, 24, 32, 48, 64, 80, 96, 112, 128, 144, 160],
    "nu_1corr":   [0.631814, 0.631046, 0.630722, 0.630350, 0.630319, 0.630285, 0.630250, 0.630140, 0.630040, 0.629850, 0.629950],
    "nu_1corr_err": [0.000018, 0.000026, 0.000033, 0.000048, 0.000062, 0.000078, 0.000100, 0.000130, 0.000150, 0.000180, 0.000220],
    "nu_2corr":   [0.630806, 0.630513, 0.630241, 0.630278, 0.630210, 0.630100, 0.629930, 0.630010, 0.630040, 0.629850, 0.629950],
    "nu_2corr_err": [0.000030, 0.000040, 0.000055, 0.000078, 0.000110, 0.000150, 0.000180, 0.000170, 0.000150, 0.000180, 0.000220],
    "nu_3corr":   [0.630072, 0.630049, 0.629980, 0.629990, 0.630060, 0.629930, 0.629900, 0.629930, 0.629840, 0.629960, None],
    "nu_3corr_err": [0.000045, 0.000057, 0.000077, 0.000110, 0.000150, 0.000210, 0.000290, 0.000180, 0.000220, 0.000260, None],
}

# Table IV: Kc vs L_min with one correction term (no fixed omega)
TABLE_IV = {
    "L_min": [16, 24, 32, 48, 64, 80, 96, 112, 128, 144, 160],
    "Kc": [0.2216546218, 0.2216546239, 0.2216546249, 0.2216546234, 0.2216546253,
           0.2216546261, 0.2216546300, 0.2216546302, 0.2216546302, 0.2216546280, 0.2216546303],
    "Kc_err": [1.3e-9, 1.6e-9, 1.9e-9, 2.7e-9, 4.5e-9, 6.2e-9, 7.8e-9, 6.9e-9, 6.3e-9, 1.3e-8, 8.5e-9],
}

# Table V: Kc vs L_min with 1, 2, 3 fixed correction exponents
TABLE_V = {
    "L_min": [16, 24, 32, 48, 64, 80, 96, 112, 128, 144, 160],
    "Kc_1corr": [0.2216546562, 0.2216546388, 0.2216546343, 0.2216546307, 0.2216546284,
                 0.2216546275, 0.2216546260, 0.2216546259, 0.2216546258, 0.2216546270, 0.2216546263],
    "Kc_1corr_err": [1.0e-9, 1.1e-9, 1.1e-9, 1.2e-9, 1.3e-9, 1.4e-9, 1.7e-9, 1.8e-9, 2.1e-9, 2.5e-9, 2.3e-9],
    "Kc_2corr": [0.2216546393, 0.2216546308, 0.2216546307, 0.2216546305, 0.2216546284,
                 0.2216546275, 0.2216546260, 0.2216546260, 0.2216546258, 0.2216546270, 0.2216546264],
    "Kc_2corr_err": [1.1e-9, 1.2e-9, 1.2e-9, 1.2e-9, 1.3e-9, 1.5e-9, 1.6e-9, 1.8e-9, 2.1e-9, 2.5e-9, 2.3e-9],
    "Kc_3corr": [0.2216546257, 0.2216546257, 0.2216546253, 0.2216546232, 0.2216546234,
                 0.2216546250, 0.2216546279, 0.2216546250, 0.2216546263, 0.2216546271, None],
    "Kc_3corr_err": [2.1e-9, 2.4e-9, 3.2e-9, 3.0e-9, 6.0e-9, 7.5e-9, 9.7e-9, 4.9e-9, 4.8e-9, 3.4e-9, None],
}

# Table VII: Kc from cumulant crossing with one correction term
TABLE_VII = {
    "L_min": [16, 24, 32, 48, 64, 80, 96, 112, 128, 144, 160, 192],
    "Kc": [0.22165462872, 0.22165462685, 0.22165462617, 0.22165462463, 0.22165462544,
           0.22165462630, 0.22165462770, 0.22165462800, 0.22165462840, 0.22165462780,
           0.22165462930, 0.22165462950],
    "Kc_err": [4.1e-10, 5.0e-10, 5.8e-10, 7.5e-10, 9.1e-10, 1.1e-9, 1.2e-9, 1.5e-9,
               1.7e-9, 2.2e-9, 2.4e-9, 3.3e-9],
    "chi2_per_dof": [1.64, 1.10, 1.06, 0.88, 0.89, 0.90, 0.84, 0.93, 1.04, 1.18, 1.29, 1.70],
}

# Table VIII: Kc from cumulant crossing with two correction terms
TABLE_VIII = {
    "L_min": [16, 24, 32, 48, 64],
    "Kc": [0.22165462483, 0.22165462490, 0.22165462450, 0.22165462463, 0.22165462540],
    "Kc_err": [9.5e-10, 1.0e-9, 8.0e-10, 8.5e-10, 1.0e-9],
}


# ---------------------------------------------------------------------------
# Wolff cluster algorithm for 3D Ising model
# ---------------------------------------------------------------------------

class Ising3D:
    """3D Ising model on a simple cubic lattice with periodic BCs."""

    def __init__(self, L, K, rng=None):
        self.L = L
        self.K = K
        self.N = L ** 3
        self.spins = np.ones(self.N, dtype=np.int8)
        self.rng = rng or np.random.default_rng()

        # Pre-compute neighbor indices for periodic boundary conditions
        self._neighbors = self._build_neighbors()

    def _build_neighbors(self):
        L = self.L
        N = self.N
        nbrs = np.zeros((N, 6), dtype=np.int32)
        for idx in range(N):
            z, rem = divmod(idx, L * L)
            y, x = divmod(rem, L)
            nbrs[idx, 0] = z * L * L + y * L + (x + 1) % L          # +x
            nbrs[idx, 1] = z * L * L + y * L + (x - 1) % L          # -x
            nbrs[idx, 2] = z * L * L + ((y + 1) % L) * L + x        # +y
            nbrs[idx, 3] = z * L * L + ((y - 1) % L) * L + x        # -y
            nbrs[idx, 4] = ((z + 1) % L) * L * L + y * L + x        # +z
            nbrs[idx, 5] = ((z - 1) % L) * L * L + y * L + x        # -z
        return nbrs

    def initialize_cold(self):
        self.spins[:] = 1

    def initialize_hot(self):
        self.spins = self.rng.choice(np.array([-1, 1], dtype=np.int8), size=self.N)

    def wolff_step(self):
        """Perform one Wolff cluster flip."""
        p_add = 1.0 - np.exp(-2.0 * self.K)
        seed = self.rng.integers(0, self.N)
        cluster_spin = self.spins[seed]

        stack = [seed]
        in_cluster = np.zeros(self.N, dtype=np.bool_)
        in_cluster[seed] = True
        cluster_size = 0

        while stack:
            site = stack.pop()
            cluster_size += 1
            for nbr in self._neighbors[site]:
                if not in_cluster[nbr] and self.spins[nbr] == cluster_spin:
                    if self.rng.random() < p_add:
                        in_cluster[nbr] = True
                        stack.append(nbr)

        # Flip the cluster
        self.spins[in_cluster] *= -1
        return cluster_size

    def energy(self):
        """Compute dimensionless energy E = -sum_{<i,j>} s_i s_j."""
        E = 0
        for dim_offsets in [(1, self.L), (self.L, self.L * self.L), (self.L * self.L, self.N)]:
            pass
        # More efficient: use neighbor list
        E = 0.0
        for i in range(self.N):
            for d in range(3):  # only positive direction to avoid double counting
                j = self._neighbors[i, 2 * d]
                E -= self.spins[i] * self.spins[j]
        return E

    def energy_fast(self):
        """Fast energy calculation using array operations."""
        L = self.L
        s = self.spins.reshape((L, L, L))
        E = -(np.sum(s * np.roll(s, 1, axis=0)) +
              np.sum(s * np.roll(s, 1, axis=1)) +
              np.sum(s * np.roll(s, 1, axis=2)))
        return int(E)

    def magnetization(self):
        """Compute total magnetization."""
        return int(np.sum(self.spins))

    def measure(self):
        """Return (energy_per_site, |magnetization|_per_site)."""
        E = self.energy_fast()
        M = abs(self.magnetization())
        return E / self.N, M / self.N


# ---------------------------------------------------------------------------
# Histogram reweighting
# ---------------------------------------------------------------------------

def histogram_reweight(energies, K0, K_target, L):
    """
    Reweight energy histogram from simulation at K0 to target K.
    energies: array of total dimensionless energy measurements.
    Returns reweighted expectation values of <E>, <E^2>, <|M|>, etc.
    """
    N = L ** 3
    dK = K0 - K_target
    weights = np.exp(dK * energies)
    Z = np.mean(weights)
    return weights / Z


def compute_u4(m2_samples, m4_samples, weights=None):
    """Compute 4th order Binder cumulant U_4 = 1 - <m^4>/(3<m^2>^2)."""
    if weights is not None:
        m2_avg = np.average(m2_samples, weights=weights)
        m4_avg = np.average(m4_samples, weights=weights)
    else:
        m2_avg = np.mean(m2_samples)
        m4_avg = np.mean(m4_samples)
    if m2_avg == 0:
        return 0.0
    return 1.0 - m4_avg / (3.0 * m2_avg ** 2)


# ---------------------------------------------------------------------------
# Finite-size scaling fitting functions
# ---------------------------------------------------------------------------

def fss_one_correction(L, X0, Kc, A0, A1, omega1, nu):
    """K_c(L) = K_c + A0 * L^(-1/nu) * (1 + A1 * L^(-omega1))"""
    return Kc + A0 * L ** (-1.0 / nu) * (1.0 + A1 * L ** (-omega1))


def fss_three_corrections(L, X0, nu, a1, a2, a3):
    """X_max = X0 * L^(1/nu) * (1 + a1*L^{-0.83} + a2*L^{-4} + a3*L^{-1.6})"""
    return X0 * L ** (1.0 / nu) * (1.0 + a1 * L ** (-0.83) + a2 * L ** (-4.0) + a3 * L ** (-1.6))


def fit_kc_vs_lmin(L_arr, Kc_arr, nu=0.629912, omega1=0.83):
    """Fit K_c(L) with one correction term, fixed nu and omega."""
    def model(L, Kc, A0, A1):
        return Kc + A0 * L ** (-1.0 / nu) * (1.0 + A1 * L ** (-omega1))
    try:
        popt, pcov = curve_fit(model, L_arr, Kc_arr, p0=[0.22165463, 1.0, 0.1])
        return popt, pcov
    except RuntimeError:
        return None, None


# ---------------------------------------------------------------------------
# Constants and parameters from the paper
# ---------------------------------------------------------------------------

KC_BEST = 0.221654626   # Best estimate for Kc
NU_BEST = 0.629912      # Best estimate for nu
GAMMA_BEST = 1.23708    # Best estimate for gamma
BETA_BEST = 0.32630     # Best estimate for beta

# Correction exponents
OMEGA_1 = 0.83          # Leading confluent correction
OMEGA_2 = 4.0           # Sub-leading confluent correction
OMEGA_NU = 1.6          # Non-linear scaling field correction

# Lattice sizes studied
LATTICE_SIZES = [16, 24, 32, 48, 64, 80, 96, 112, 128, 144, 160, 192, 256, 384, 512, 768, 1024]

# Simulation parameters
K0_SIM = 0.221654       # Simulation coupling (near Kc)
