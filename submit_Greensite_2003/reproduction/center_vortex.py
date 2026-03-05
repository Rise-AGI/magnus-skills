"""
Random Center Vortex Model for SU(2) Lattice Gauge Theory

Implements the Z_2 random center vortex model described in Greensite (2003),
Section 5. In this model, thin center vortices pierce plaquettes randomly,
and a Wilson loop receives a factor of -1 for each vortex piercing.

The key prediction is that Wilson loops follow an area law:
    W(C) = exp(-sigma * A)
where sigma = -ln(1 - 2*f) and f is the probability of a plaquette being
pierced by a vortex.
"""

import numpy as np


def vortex_wilson_loop_analytic(area, f):
    """
    Analytical prediction for Wilson loop in the random vortex model.

    W(C) = (1 - 2f)^A

    where f is the vortex piercing probability per plaquette and
    A is the minimal area of the loop.

    Parameters
    ----------
    area : float or array
        Minimal area enclosed by the Wilson loop (in lattice units).
    f : float
        Probability that a vortex pierces a given plaquette.

    Returns
    -------
    float or array
        Wilson loop value.
    """
    return (1.0 - 2.0 * f) ** area


def vortex_string_tension(f):
    """
    String tension from the random vortex model.

    sigma = -ln(1 - 2f)

    Parameters
    ----------
    f : float
        Vortex piercing probability.

    Returns
    -------
    float
        String tension in lattice units.
    """
    return -np.log(1.0 - 2.0 * f)


def simulate_vortex_model_2d(Lx, Ly, f, n_configs, rng):
    """
    Monte Carlo simulation of the Z_2 random vortex model on a 2D lattice.

    Each plaquette is independently pierced by a vortex with probability f.
    Wilson loops are computed by counting the number of vortex piercings.

    Parameters
    ----------
    Lx, Ly : int
        Lattice dimensions.
    f : float
        Vortex piercing probability.
    n_configs : int
        Number of independent configurations to average over.
    rng : numpy.random.Generator
        Random number generator.

    Returns
    -------
    dict
        Wilson loop values for various loop sizes.
    """
    max_R = min(Lx, Ly) // 2
    wilson_sums = {}
    wilson_counts = {}

    for S in range(1, max_R + 1):
        wilson_sums[S] = 0.0
        wilson_counts[S] = 0

    for _ in range(n_configs):
        # Generate random vortex configuration
        # vortex[x, y] = 1 if pierced, 0 otherwise
        vortex = (rng.random((Lx, Ly)) < f).astype(int)

        for S in range(1, max_R + 1):
            for x0 in range(Lx):
                for y0 in range(Ly):
                    # Count vortices piercing the S x S loop starting at (x0, y0)
                    n_piercings = 0
                    for dx in range(S):
                        for dy in range(S):
                            n_piercings += vortex[(x0 + dx) % Lx, (y0 + dy) % Ly]
                    # Wilson loop = (-1)^n_piercings
                    wilson_sums[S] += (-1.0) ** n_piercings
                    wilson_counts[S] += 1

    results = {}
    for S in wilson_sums:
        results[S] = wilson_sums[S] / wilson_counts[S]
    return results


def simulate_vortex_model_4d(L, f, n_configs, rng, max_loop=4):
    """
    Monte Carlo simulation of Z_2 random vortex model on a 4D lattice.

    In 4D, vortices are surfaces on the dual lattice. For each plane (mu,nu),
    each plaquette is independently pierced with probability f.
    Wilson loops in a given plane are determined by the piercings in that plane.

    Parameters
    ----------
    L : int
        Linear lattice size (L^4).
    f : float
        Vortex piercing probability per plaquette per plane.
    n_configs : int
        Number of configurations.
    rng : numpy.random.Generator
        Random number generator.
    max_loop : int
        Maximum loop size to compute.

    Returns
    -------
    dict
        Wilson loop values keyed by loop size.
    """
    wilson_sums = {S: 0.0 for S in range(1, max_loop + 1)}
    wilson_counts = {S: 0 for S in range(1, max_loop + 1)}

    n_planes = 6  # 4 choose 2

    for _ in range(n_configs):
        # For each plane orientation, generate independent piercings
        for plane_idx in range(n_planes):
            # 2D slice of piercings
            piercings = (rng.random((L, L)) < f).astype(int)

            for S in range(1, max_loop + 1):
                for x0 in range(L):
                    for y0 in range(L):
                        n_p = 0
                        for dx in range(S):
                            for dy in range(S):
                                n_p += piercings[(x0+dx)%L, (y0+dy)%L]
                        wilson_sums[S] += (-1.0)**n_p
                        wilson_counts[S] += 1

    return {S: wilson_sums[S] / wilson_counts[S] for S in wilson_sums}
