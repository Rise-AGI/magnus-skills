"""
2D Ising Model Monte Carlo Simulation using Wolff Cluster Algorithm.

Reproduces results from Sandvik (2010), Section 3:
- Average squared magnetization vs temperature (Fig 13)
- Susceptibility vs temperature (Fig 14)
- Binder cumulant vs temperature (Fig 16)

The 2D Ising model Hamiltonian:
  H = -J sum_{<i,j>} s_i s_j
where s_i = +/- 1 and the sum is over nearest neighbors on a square lattice
with periodic boundary conditions.

Exact critical temperature: T_c / J = 2 / ln(1 + sqrt(2)) ~ 2.26919
Exact critical exponents: nu = 1, gamma = 7/4, beta = 1/8, eta = 1/4
"""

import numpy as np


class IsingModel:
    """2D Ising model on L x L square lattice with periodic BCs."""

    def __init__(self, L, T, J=1.0, seed=None):
        self.L = L
        self.N = L * L
        self.T = T
        self.J = J
        self.beta = 1.0 / T
        self.rng = np.random.default_rng(seed)
        # Initialize all spins up (cold start)
        self.spins = np.ones((L, L), dtype=np.int8)
        # Precompute neighbor indices
        self._build_neighbors()
        # Wolff cluster probability
        self.p_add = 1.0 - np.exp(-2.0 * self.J * self.beta)

    def _build_neighbors(self):
        L = self.L
        self.neighbors = np.zeros((L, L, 4, 2), dtype=np.int32)
        for x in range(L):
            for y in range(L):
                self.neighbors[x, y, 0] = [(x + 1) % L, y]
                self.neighbors[x, y, 1] = [(x - 1) % L, y]
                self.neighbors[x, y, 2] = [x, (y + 1) % L]
                self.neighbors[x, y, 3] = [x, (y - 1) % L]

    def set_temperature(self, T):
        self.T = T
        self.beta = 1.0 / T
        self.p_add = 1.0 - np.exp(-2.0 * self.J * self.beta)

    def cold_start(self):
        self.spins[:] = 1

    def hot_start(self):
        self.spins = self.rng.choice(np.array([-1, 1], dtype=np.int8),
                                      size=(self.L, self.L))

    def wolff_step(self):
        """One Wolff cluster flip."""
        L = self.L
        rng = self.rng
        # Pick random seed site
        sx, sy = rng.integers(0, L, size=2)
        seed_spin = self.spins[sx, sy]
        cluster = [(sx, sy)]
        in_cluster = np.zeros((L, L), dtype=bool)
        in_cluster[sx, sy] = True
        idx = 0
        while idx < len(cluster):
            cx, cy = cluster[idx]
            for d in range(4):
                nx, ny = self.neighbors[cx, cy, d]
                if not in_cluster[nx, ny] and self.spins[nx, ny] == seed_spin:
                    if rng.random() < self.p_add:
                        cluster.append((nx, ny))
                        in_cluster[nx, ny] = True
            idx += 1
        # Flip the cluster
        for cx, cy in cluster:
            self.spins[cx, cy] = -self.spins[cx, cy]
        return len(cluster)

    def sweep(self, n_clusters=None):
        """Perform enough cluster flips to flip ~N spins on average."""
        if n_clusters is None:
            total_flipped = 0
            while total_flipped < self.N:
                total_flipped += self.wolff_step()
        else:
            for _ in range(n_clusters):
                self.wolff_step()

    def magnetization(self):
        return np.sum(self.spins) / self.N

    def magnetization_abs(self):
        return np.abs(np.sum(self.spins)) / self.N

    def energy_per_site(self):
        s = self.spins
        L = self.L
        # Sum over bonds (each counted once)
        e = 0.0
        e -= np.sum(s * np.roll(s, -1, axis=0))  # x-direction
        e -= np.sum(s * np.roll(s, -1, axis=1))  # y-direction
        return self.J * e / self.N

    def measure(self):
        """Return magnetization and energy per site."""
        m = self.magnetization()
        e = self.energy_per_site()
        return m, e


def run_ising_simulation(L, T, n_therm=200, n_meas=2000, seed=None):
    """Run Ising MC simulation and return measurements."""
    model = IsingModel(L, T, seed=seed)
    model.hot_start()

    # Thermalize
    for _ in range(n_therm):
        model.sweep()

    # Measure
    m_list = []
    e_list = []
    for _ in range(n_meas):
        model.sweep()
        m, e = model.measure()
        m_list.append(m)
        e_list.append(e)

    m_arr = np.array(m_list)
    e_arr = np.array(e_list)

    m2 = np.mean(m_arr**2)
    m4 = np.mean(m_arr**4)
    m_abs = np.mean(np.abs(m_arr))
    chi = model.N * (np.mean(m_arr**2) - np.mean(np.abs(m_arr))**2) / T
    # Binder cumulant: U_2 = 1 - <m^4> / (3 <m^2>^2)
    binder = 1.0 - m4 / (3.0 * m2**2) if m2 > 0 else 0.0

    return {
        'm2': m2,
        'm4': m4,
        'm_abs': m_abs,
        'chi': chi,
        'binder': binder,
        'energy': np.mean(e_arr),
    }
