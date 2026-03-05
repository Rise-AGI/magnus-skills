"""
Core module for the Bolder et al. (2000) lattice QCD paper reproduction.
Implements:
  - All-R Approach (ARA) vector enumeration (Table I)
  - Bresenham algorithm for lattice path construction (Section II)
  - Cornell potential parametrization (Section III)
"""

import numpy as np


# ===========================================================================
# R-vector enumeration (Table I)
# ===========================================================================

def enumerate_r_vectors(R_min, R_max, C_max_limit):
    """
    Enumerate all integer 3-vectors R = (C1, C2, C3) satisfying:
      R_min^2 <= C1^2 + C2^2 + C3^2 <= R_max^2
      |Ci| <= C_max_limit
    Returns dict mapping R^2 -> list of (Cmin, Cmid, Cmax) tuples (sorted).
    """
    vectors_by_R2 = {}
    for c1 in range(-C_max_limit, C_max_limit + 1):
        for c2 in range(-C_max_limit, C_max_limit + 1):
            for c3 in range(-C_max_limit, C_max_limit + 1):
                R2 = c1*c1 + c2*c2 + c3*c3
                if R_min**2 <= R2 <= R_max**2:
                    if R2 not in vectors_by_R2:
                        vectors_by_R2[R2] = []
                    vectors_by_R2[R2].append((c1, c2, c3))
    return vectors_by_R2


def classify_vectors(vectors_by_R2):
    """
    Classify vectors into 'standard' (multiples of (1,0,0),(1,1,0),(1,1,1),
    (2,1,0),(2,1,1),(2,2,1)) and 'ARA' (all vectors).
    Returns (standard_stats, ara_stats) dicts with R-value counts.
    """
    standard_directions = [
        (1, 0, 0), (1, 1, 0), (1, 1, 1),
        (2, 1, 0), (2, 1, 1), (2, 2, 1)
    ]

    def is_standard(c1, c2, c3):
        """Check if (c1,c2,c3) is a multiple of a standard direction."""
        absc = sorted([abs(c1), abs(c2), abs(c3)])
        for d in standard_directions:
            ds = sorted(d)
            if ds[0] == 0 and ds[1] == 0:
                # (n,0,0) type
                if absc[0] == 0 and absc[1] == 0:
                    return True
            elif ds[0] == 0:
                # (n,n,0) or (2n,n,0) type
                if absc[0] == 0 and absc[1] > 0 and absc[2] > 0:
                    r1 = absc[1] / ds[1]
                    r2 = absc[2] / ds[2]
                    if r1 == r2 and r1 == int(r1):
                        return True
            else:
                if absc[0] > 0:
                    r0 = absc[0] / ds[0]
                    r1 = absc[1] / ds[1]
                    r2 = absc[2] / ds[2]
                    if r0 == r1 == r2 and r0 == int(r0):
                        return True
        return False

    standard_by_R2 = {}
    for R2, vecs in vectors_by_R2.items():
        std_vecs = [v for v in vecs if is_standard(*v)]
        if std_vecs:
            standard_by_R2[R2] = std_vecs

    n_R_ara = len(vectors_by_R2)
    n_vec_ara = sum(len(v) for v in vectors_by_R2.values())
    n_R_std = len(standard_by_R2)
    n_vec_std = sum(len(v) for v in standard_by_R2.values())

    return {
        'standard': {'n_R': n_R_std, 'n_vectors': n_vec_std,
                     'avg_per_R': n_vec_std / max(n_R_std, 1)},
        'ARA': {'n_R': n_R_ara, 'n_vectors': n_vec_ara,
                'avg_per_R': n_vec_ara / max(n_R_ara, 1)}
    }


# ===========================================================================
# Bresenham algorithm (Section II, Figure 1)
# ===========================================================================

def bresenham_2d(cmax, cmin):
    """
    2D Bresenham path from (0,0) to (cmax, cmin).
    Both cmax, cmin >= 0 and cmax >= cmin.
    Returns list of (x, y) positions along the path.
    """
    path = [(0, 0)]
    x, y = 0, 0
    chi = 2 * cmin - cmax
    for _ in range(cmax):
        x += 1
        if chi >= 0:
            chi -= 2 * cmax
            y += 1
        chi += 2 * cmin
        path.append((x, y))
    return path


def bresenham_3d(c1, c2, c3):
    """
    3D Bresenham path from origin to (c1, c2, c3).
    All components >= 0, c1 >= c2 >= c3 (max, mid, min).
    Returns list of (x, y, z) positions.
    """
    path = [(0, 0, 0)]
    x, y, z = 0, 0, 0
    chi_mid = 2 * c2 - c1
    chi_min = 2 * c3 - c1
    for _ in range(c1):
        x += 1
        if chi_mid >= 0:
            chi_mid -= 2 * c1
            y += 1
        if chi_min >= 0:
            chi_min -= 2 * c1
            z += 1
        chi_mid += 2 * c2
        chi_min += 2 * c3
        path.append((x, y, z))
    return path


def bresenham_general(C):
    """
    General Bresenham for arbitrary integer vector C = (c1, c2, c3).
    Handles signs and ordering internally.
    Returns list of (x, y, z) positions.
    """
    c = list(C)
    signs = [1 if ci >= 0 else -1 for ci in c]
    absc = [abs(ci) for ci in c]

    # Sort by magnitude: find ordering
    order = sorted(range(3), key=lambda i: absc[i], reverse=True)
    sorted_c = [absc[order[i]] for i in range(3)]

    raw_path = bresenham_3d(sorted_c[0], sorted_c[1], sorted_c[2])

    # Unsort and apply signs
    path = []
    for pt in raw_path:
        out = [0, 0, 0]
        for i in range(3):
            out[order[i]] = pt[i] * signs[order[i]]
        path.append(tuple(out))
    return path


# ===========================================================================
# Cornell potential (Section III)
# ===========================================================================

def cornell_potential(R, V0, K, e_param, R0=5.89):
    """
    Cornell potential: V(R) = V0 + K * R - e / R
    R in lattice units, R0 = Sommer radius in lattice units.
    Returns V(R) in lattice units.
    """
    return V0 + K * R - e_param / R


def string_breaking_threshold(m_PS_a):
    """
    String breaking threshold energy = 2 * m_PS * a.
    m_PS_a: pseudoscalar meson mass in lattice units.
    """
    return 2.0 * m_PS_a
