"""
Core 2D Ising Model with Metropolis Algorithm (checkerboard decomposition).

Implements the model from:
  J. Kotze, "Introduction to Monte Carlo methods for an Ising Model of a
  Ferromagnet", arXiv:0803.0217 (2008).

Uses vectorized numpy operations via checkerboard (black/white) sublattice
updates for performance.  Only numpy, scipy, matplotlib allowed.
"""

import numpy as np


# ---------------------------------------------------------------------------
# Lattice helpers
# ---------------------------------------------------------------------------

def init_lattice(L, cold=False, rng=None):
    """Return an L x L lattice of +/-1 spins."""
    if rng is None:
        rng = np.random.default_rng()
    if cold:
        return np.ones((L, L), dtype=np.int8)
    return rng.choice(np.array([-1, 1], dtype=np.int8), size=(L, L))


def _neighbour_sum(lat):
    """Sum of four nearest neighbours for every site (periodic BC)."""
    return (np.roll(lat, 1, axis=0) + np.roll(lat, -1, axis=0) +
            np.roll(lat, 1, axis=1) + np.roll(lat, -1, axis=1))


# ---------------------------------------------------------------------------
# Single checkerboard sweep
# ---------------------------------------------------------------------------

def metropolis_sweep(lat, beta, rng):
    """One full Metropolis sweep using checkerboard decomposition.

    Black sites (i+j even) are updated first, then white sites (i+j odd).
    Each sublattice update is fully vectorized.
    """
    L = lat.shape[0]
    rows, cols = np.meshgrid(np.arange(L), np.arange(L), indexing='ij')

    for parity in (0, 1):
        mask = (rows + cols) % 2 == parity
        nn = _neighbour_sum(lat)
        dE = 2 * lat * nn          # energy change if flipped
        rand = rng.random((L, L))
        accept = (dE <= 0) | (rand < np.exp(-beta * dE))
        flip = mask & accept
        lat[flip] *= -1


# ---------------------------------------------------------------------------
# Full simulation
# ---------------------------------------------------------------------------

def run_simulation(L, temperatures, n_eq=2000, n_mc=10000, seed=42):
    """Run Ising MC for a range of temperatures.

    Parameters
    ----------
    L : int
        Lattice linear size.
    temperatures : array-like
        Temperatures to simulate (in units of J/k_B).
    n_eq : int
        Number of equilibration sweeps per temperature.
    n_mc : int
        Number of measurement sweeps per temperature.
    seed : int
        RNG seed.

    Returns
    -------
    dict with keys:
        T, E_avg, E2_avg, M_avg, Mabs_avg, M2_avg, M4_avg, C, chi_prime, U_L
    """
    rng = np.random.default_rng(seed)
    N = L * L
    temps = np.asarray(temperatures, dtype=float)

    results = {k: np.zeros(len(temps)) for k in
               ['T', 'E_avg', 'E2_avg', 'M_avg', 'Mabs_avg',
                'M2_avg', 'M4_avg', 'C', 'chi_prime', 'chi', 'U_L']}
    results['T'] = temps.copy()

    # Start from high-T random config; sweep T from high to low
    idx_order = np.argsort(-temps)  # descending
    lat = init_lattice(L, cold=False, rng=rng)

    for count, idx in enumerate(idx_order):
        T = temps[idx]
        beta = 1.0 / T

        # Equilibrate
        for _ in range(n_eq):
            metropolis_sweep(lat, beta, rng)

        # Measure
        e_sum = 0.0
        e2_sum = 0.0
        m_sum = 0.0
        mabs_sum = 0.0
        m2_sum = 0.0
        m4_sum = 0.0

        for _ in range(n_mc):
            metropolis_sweep(lat, beta, rng)
            nn = _neighbour_sum(lat)
            E = -0.5 * np.sum(lat * nn)  # factor 1/2 avoids double-counting
            M = np.sum(lat)
            e = E / N
            m = M / N
            e_sum += e
            e2_sum += e * e
            m_sum += m
            mabs_sum += abs(m)
            m2_sum += m * m
            m4_sum += m * m * m * m

        inv = 1.0 / n_mc
        e_avg = e_sum * inv
        e2_avg = e2_sum * inv
        m_avg = m_sum * inv
        mabs_avg = mabs_sum * inv
        m2_avg = m2_sum * inv
        m4_avg = m4_sum * inv

        results['E_avg'][idx] = e_avg
        results['E2_avg'][idx] = e2_avg
        results['M_avg'][idx] = m_avg
        results['Mabs_avg'][idx] = mabs_avg
        results['M2_avg'][idx] = m2_avg
        results['M4_avg'][idx] = m4_avg
        results['C'][idx] = (e2_avg - e_avg**2) * N / (T**2)
        results['chi_prime'][idx] = (m2_avg - mabs_avg**2) * N / T
        results['chi'][idx] = (m2_avg - m_avg**2) * N / T
        denom = 3.0 * m2_avg**2
        results['U_L'][idx] = 1.0 - m4_avg / denom if denom > 0 else 0.0

        if (count + 1) % 10 == 0:
            print(f"  L={L}: {count+1}/{len(temps)} temperatures done (T={T:.2f})")

    return results


# ---------------------------------------------------------------------------
# Exact 2x2 solutions
# ---------------------------------------------------------------------------

def exact_2x2(temperatures):
    """Exact analytic results for 2x2 Ising model.

    Returns dict with T, E_avg (per spin), C (per spin),
    Mabs_avg (per spin), chi_prime (per spin).
    """
    temps = np.asarray(temperatures, dtype=float)
    beta = 1.0 / temps
    N = 4

    Z = 2 * np.exp(8 * beta) + 12 + 2 * np.exp(-8 * beta)
    E_total = -(2 * 8 * np.exp(8 * beta) + 2 * (-8) * np.exp(-8 * beta)) / Z
    E2_total = (2 * 64 * np.exp(8 * beta) + 2 * 64 * np.exp(-8 * beta)) / Z
    Mabs_total = (2 * 4 * np.exp(8 * beta) + 8 * 2) / Z
    M2_total = (2 * 16 * np.exp(8 * beta) + 8 * 4) / Z

    E_per_spin = E_total / N
    C_per_spin = (E2_total - E_total**2) / (N * temps**2)
    Mabs_per_spin = Mabs_total / N
    M2_per_spin = M2_total / N**2
    Mabs_sq_per_spin = (Mabs_total / N)**2
    chi_prime_per_spin = (M2_per_spin - Mabs_sq_per_spin) * N / temps

    return {
        'T': temps,
        'E_avg': E_per_spin,
        'C': C_per_spin,
        'Mabs_avg': Mabs_per_spin,
        'chi_prime': chi_prime_per_spin,
    }
