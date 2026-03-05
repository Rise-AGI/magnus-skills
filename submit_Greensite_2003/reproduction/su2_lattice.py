"""
SU(2) Lattice Gauge Theory Monte Carlo Simulation
Based on: Greensite, J. (2003). The confinement problem in lattice gauge theory.
          Prog. Part. Nucl. Phys. 51, 1-83.

Core module implementing the heat bath algorithm for SU(2) lattice gauge theory
on a 4D hypercubic lattice with periodic boundary conditions, plus
analytical strong/weak coupling expansions.
"""

import numpy as np
from typing import Tuple, Optional


# ============================================================
# SU(2) quaternion algebra
# ============================================================

def su2_identity():
    return np.array([1.0, 0.0, 0.0, 0.0])

def su2_random(rng):
    a = rng.standard_normal(4)
    return a / np.linalg.norm(a)

def su2_multiply(u1, u2):
    a0, a1, a2, a3 = u1
    b0, b1, b2, b3 = u2
    return np.array([
        a0*b0 - a1*b1 - a2*b2 - a3*b3,
        a0*b1 + a1*b0 + a2*b3 - a3*b2,
        a0*b2 + a2*b0 + a3*b1 - a1*b3,
        a0*b3 + a3*b0 + a1*b2 - a2*b1
    ])

def su2_inverse(u):
    return np.array([u[0], -u[1], -u[2], -u[3]])

def su2_trace(u):
    return 2.0 * u[0]


# ============================================================
# Lattice class
# ============================================================

class SU2Lattice:
    def __init__(self, L, beta, seed=None):
        self.L = L
        self.beta = beta
        self.rng = np.random.default_rng(seed)
        self.links = np.zeros((L, L, L, L, 4, 4))

    def initialize_cold(self):
        for idx in np.ndindex(self.L, self.L, self.L, self.L):
            for mu in range(4):
                self.links[idx[0], idx[1], idx[2], idx[3], mu] = su2_identity()

    def initialize_hot(self):
        for idx in np.ndindex(self.L, self.L, self.L, self.L):
            for mu in range(4):
                self.links[idx[0], idx[1], idx[2], idx[3], mu] = su2_random(self.rng)

    def _shift(self, site, mu, d=1):
        s = list(site)
        s[mu] = (s[mu] + d) % self.L
        return tuple(s)

    def get_link(self, site, mu):
        return self.links[site[0], site[1], site[2], site[3], mu]

    def set_link(self, site, mu, u):
        self.links[site[0], site[1], site[2], site[3], mu] = u

    def get_staple_sum(self, site, mu):
        staple_sum = np.zeros(4)
        for nu in range(4):
            if nu == mu:
                continue
            site_mu = self._shift(site, mu)
            site_nu = self._shift(site, nu)
            u1 = self.get_link(site_mu, nu)
            u2 = su2_inverse(self.get_link(site_nu, mu))
            u3 = su2_inverse(self.get_link(site, nu))
            staple_sum += su2_multiply(su2_multiply(u1, u2), u3)

            site_mu_mnu = self._shift(site_mu, nu, -1)
            site_mnu = self._shift(site, nu, -1)
            u1 = su2_inverse(self.get_link(site_mu_mnu, nu))
            u2 = su2_inverse(self.get_link(site_mnu, mu))
            u3 = self.get_link(site_mnu, nu)
            staple_sum += su2_multiply(su2_multiply(u1, u2), u3)
        return staple_sum

    def heat_bath_update(self, site, mu):
        staple_sum = self.get_staple_sum(site, mu)
        k = np.sqrt(np.sum(staple_sum**2))
        if k < 1e-10:
            self.set_link(site, mu, su2_random(self.rng))
            return
        u_bar = staple_sum / k
        beta_k = self.beta * k
        while True:
            x_min = np.exp(-2.0 * beta_k)
            x = x_min + (1.0 - x_min) * self.rng.random()
            a0 = 1.0 + np.log(x) / beta_k
            if self.rng.random() < np.sqrt(max(0.0, 1.0 - a0**2)):
                break
        r = np.sqrt(max(0.0, 1.0 - a0**2))
        phi = 2.0 * np.pi * self.rng.random()
        cos_theta = 2.0 * self.rng.random() - 1.0
        sin_theta = np.sqrt(max(0.0, 1.0 - cos_theta**2))
        u = np.array([a0, r*sin_theta*np.cos(phi), r*sin_theta*np.sin(phi), r*cos_theta])
        u_new = su2_multiply(u, su2_inverse(u_bar))
        u_new = u_new / np.linalg.norm(u_new)
        self.set_link(site, mu, u_new)

    def sweep(self):
        for idx in np.ndindex(self.L, self.L, self.L, self.L):
            for mu in range(4):
                self.heat_bath_update(idx, mu)

    def plaquette_value(self, site, mu, nu):
        site_mu = self._shift(site, mu)
        site_nu = self._shift(site, nu)
        u1 = self.get_link(site, mu)
        u2 = self.get_link(site_mu, nu)
        u3 = su2_inverse(self.get_link(site_nu, mu))
        u4 = su2_inverse(self.get_link(site, nu))
        plaq = su2_multiply(su2_multiply(su2_multiply(u1, u2), u3), u4)
        return 0.5 * su2_trace(plaq)

    def average_plaquette(self):
        total = 0.0
        count = 0
        for idx in np.ndindex(self.L, self.L, self.L, self.L):
            for mu in range(4):
                for nu in range(mu+1, 4):
                    total += self.plaquette_value(idx, mu, nu)
                    count += 1
        return total / count

    def wilson_loop(self, R, T):
        total = 0.0
        count = 0
        for start in np.ndindex(self.L, self.L, self.L, self.L):
            for mu in range(4):
                for nu in range(mu+1, 4):
                    product = su2_identity()
                    site = start
                    for _ in range(R):
                        product = su2_multiply(product, self.get_link(site, mu))
                        site = self._shift(site, mu)
                    for _ in range(T):
                        product = su2_multiply(product, self.get_link(site, nu))
                        site = self._shift(site, nu)
                    for _ in range(R):
                        site = self._shift(site, mu, -1)
                        product = su2_multiply(product, su2_inverse(self.get_link(site, mu)))
                    for _ in range(T):
                        site = self._shift(site, nu, -1)
                        product = su2_multiply(product, su2_inverse(self.get_link(site, nu)))
                    total += 0.5 * su2_trace(product)
                    count += 1
        return total / count

    def square_wilson_loop(self, S):
        return self.wilson_loop(S, S)


# ============================================================
# Analytical predictions
# ============================================================

def strong_coupling_plaquette(beta):
    """Strong coupling expansion for average plaquette (1/2)Tr(U_plaq).
    To leading orders: <P> = beta/4 + (beta/4)^2/2 + ...
    where P = (1/2)Tr(U_plaq) and action plaquette = 1 - P.
    """
    x = beta / 4.0
    return x * (1.0 + x/2.0 + x**2/3.0)

def weak_coupling_plaquette(beta):
    """Weak coupling (perturbative) prediction for average plaquette.
    <P> = 1 - 3/(4*beta) + ...
    """
    return 1.0 - 3.0 / (4.0 * beta)

def strong_coupling_string_tension(beta):
    """Strong coupling expansion for string tension.
    a^2 * K = -ln(beta/4) + O(beta^2)
    """
    return -np.log(beta / 4.0)

def asymptotic_freedom_string_tension(beta):
    """Asymptotic freedom (weak coupling) prediction for string tension.
    a^2 * K ~ C * exp(-6*pi^2/11 * (beta - 2))
    Normalization chosen so curves meet near beta ~ 2.2
    """
    return 0.5 * np.exp(-6.0 * np.pi**2 / 11.0 * (beta - 2.0))


def run_simulation(L, beta, n_therm, n_meas, start='cold', seed=None):
    """Run MC simulation and return measurements."""
    lattice = SU2Lattice(L, beta, seed)
    if start == 'cold':
        lattice.initialize_cold()
    else:
        lattice.initialize_hot()

    plaquettes_all = []
    for i in range(n_therm):
        lattice.sweep()
        plaquettes_all.append(lattice.average_plaquette())

    plaquettes = []
    for i in range(n_meas):
        lattice.sweep()
        plaquettes.append(lattice.average_plaquette())

    return {
        'plaquettes_therm': np.array(plaquettes_all),
        'plaquettes': np.array(plaquettes),
        'lattice': lattice,
        'beta': beta,
        'L': L
    }
