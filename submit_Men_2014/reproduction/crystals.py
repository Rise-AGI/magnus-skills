"""
Crystal structure definitions for 3D photonic crystals.

Men et al., Optics Express 22(19), 22632-22657 (2014)

Defines parametric dielectric functions for:
- SC5: Simple cubic hollow spheres + cylinders (Fig. 5)
- Diamond2: Diamond lattice cylinder bonds (Fig. 6)
- FCC8: FCC hollow spheres + rods to 12 neighbors (Fig. 7)

All epsilon functions take fractional coordinates (s1, s2, s3) in [0,1)
of the PRIMITIVE cell and return the dielectric constant.
"""

import numpy as np
from pwe3d import reciprocal_vectors


# ==================== Lattice definitions ====================

def sc_lattice():
    """Simple cubic lattice vectors (a=1)."""
    return (np.array([1, 0, 0], dtype=float),
            np.array([0, 1, 0], dtype=float),
            np.array([0, 0, 1], dtype=float))


def fcc_lattice():
    """FCC primitive lattice vectors (conventional a=1)."""
    return (0.5 * np.array([0, 1, 1], dtype=float),
            0.5 * np.array([1, 0, 1], dtype=float),
            0.5 * np.array([1, 1, 0], dtype=float))


# ==================== Helper: minimum image distance ====================

def _min_image_disp(pos_cart, ref_cart, A, A_inv):
    """
    Compute minimum image displacement vector from ref to pos.

    Parameters
    ----------
    pos_cart : (..., 3) positions in Cartesian
    ref_cart : (3,) reference point in Cartesian
    A : (3, 3) matrix [a1, a2, a3] column vectors
    A_inv : (3, 3) inverse

    Returns
    -------
    dp : (..., 3) minimum image displacement in Cartesian
    """
    dp = pos_cart - ref_cart
    # Convert to fractional
    ds = np.einsum('ij,...j->...i', A_inv, dp)
    # Wrap to [-0.5, 0.5)
    ds = ds - np.round(ds)
    # Convert back to Cartesian
    dp = np.einsum('ij,...j->...i', A, ds)
    return dp


def _point_to_cylinder_dist(pos_cart, start_cart, bond_cart, A, A_inv):
    """
    Compute distance from points to a cylinder (line segment).

    Parameters
    ----------
    pos_cart : (..., 3) positions
    start_cart : (3,) cylinder start
    bond_cart : (3,) bond vector from start to end
    A, A_inv : lattice matrices for periodic wrapping

    Returns
    -------
    perp_dist : (...) perpendicular distance to cylinder axis
    proj : (...) projection along bond (0=start, bond_length=end)
    """
    bond_len = np.linalg.norm(bond_cart)
    bond_hat = bond_cart / bond_len

    # Minimum image displacement from start to pos
    dp = _min_image_disp(pos_cart, start_cart, A, A_inv)

    # Project onto bond axis
    proj = np.einsum('...i,i->...', dp, bond_hat)

    # Perpendicular component
    perp_vec = dp - proj[..., np.newaxis] * bond_hat
    perp_dist = np.linalg.norm(perp_vec, axis=-1)

    return perp_dist, proj, bond_len


# ==================== Epsilon functions ====================

def make_sc5_eps(r1=0.14, r2=0.36, r3=0.105, eps_hi=12.96, eps_lo=1.0):
    """
    SC5 structure epsilon function.

    Hollow dielectric spheres at SC lattice sites connected by dielectric cylinders
    along the three cubic axes.

    Parameters: r1=inner air radius, r2=outer sphere radius, r3=cylinder radius.
    Optimal: (0.14, 0.36, 0.105) giving ~17% gap (Fig. 5b).
    """
    def eps_func(s1, s2, s3):
        # SC: fractional = Cartesian (a=1)
        x = s1 - np.round(s1)
        y = s2 - np.round(s2)
        z = s3 - np.round(s3)

        eps = np.full_like(x, eps_lo, dtype=float)

        r = np.sqrt(x**2 + y**2 + z**2)
        eps = np.where(r < r2, eps_hi, eps)

        # Cylinders along x, y, z axes
        eps = np.where(np.sqrt(y**2 + z**2) < r3, eps_hi, eps)
        eps = np.where(np.sqrt(x**2 + z**2) < r3, eps_hi, eps)
        eps = np.where(np.sqrt(x**2 + y**2) < r3, eps_hi, eps)

        # Air hole
        eps = np.where(r < r1, eps_lo, eps)

        return eps

    return eps_func


def make_diamond2_eps(r=0.1, eps_hi=12.96, eps_lo=1.0):
    """
    Diamond2 structure: diamond lattice of dielectric cylinder bonds.

    Two atoms per FCC primitive cell at:
      Atom 1: (0,0,0) in Cartesian
      Atom 2: (1/4, 1/4, 1/4) in Cartesian

    Four tetrahedral bonds from Atom 1 to Atom 2 images:
      (1/4, 1/4, 1/4), (1/4, -1/4, -1/4), (-1/4, 1/4, -1/4), (-1/4, -1/4, 1/4)

    Cylinder radius r in units of conventional cube edge a.
    Optimal: r=0.1 giving ~31.56% gap (Fig. 6b).
    """
    a1 = 0.5 * np.array([0, 1, 1], dtype=float)
    a2 = 0.5 * np.array([1, 0, 1], dtype=float)
    a3 = 0.5 * np.array([1, 1, 0], dtype=float)
    A = np.column_stack([a1, a2, a3])
    A_inv = np.linalg.inv(A)

    # Atom positions in Cartesian
    atom1 = np.array([0.0, 0.0, 0.0])
    atom2 = np.array([0.25, 0.25, 0.25])

    # Bond vectors from atom 1 to 4 neighbors (atom 2 images)
    bonds_from_1 = np.array([
        [1, 1, 1],
        [1, -1, -1],
        [-1, 1, -1],
        [-1, -1, 1]
    ], dtype=float) / 4.0

    # Bond vectors from atom 2 to 4 neighbors (atom 1 images)
    bonds_from_2 = np.array([
        [-1, -1, -1],
        [-1, 1, 1],
        [1, -1, 1],
        [1, 1, -1]
    ], dtype=float) / 4.0

    def eps_func(s1, s2, s3):
        shape = s1.shape
        eps = np.full(shape, eps_lo, dtype=float)

        # Convert fractional to Cartesian
        x = s1 * a1[0] + s2 * a2[0] + s3 * a3[0]
        y = s1 * a1[1] + s2 * a2[1] + s3 * a3[1]
        z = s1 * a1[2] + s2 * a2[2] + s3 * a3[2]
        pos = np.stack([x, y, z], axis=-1)  # (..., 3)

        # Bonds from atom 1
        for bond in bonds_from_1:
            perp_dist, proj, blen = _point_to_cylinder_dist(
                pos, atom1, bond, A, A_inv)
            in_cyl = (proj >= 0) & (proj <= blen) & (perp_dist < r)
            eps = np.where(in_cyl, eps_hi, eps)

        # Bonds from atom 2
        for bond in bonds_from_2:
            perp_dist, proj, blen = _point_to_cylinder_dist(
                pos, atom2, bond, A, A_inv)
            in_cyl = (proj >= 0) & (proj <= blen) & (perp_dist < r)
            eps = np.where(in_cyl, eps_hi, eps)

        return eps

    return eps_func


def make_fcc8_eps(r1=0.12, r2=0.19, r3=0.08, eps_hi=12.96, eps_lo=1.0):
    """
    FCC8 structure: FCC lattice of hollow dielectric spheres connected by rods
    to 12 nearest neighbors.

    One atom per FCC primitive cell at origin.
    12 FCC nearest neighbors: displacements (±1,±1,0)/2, (±1,0,±1)/2, (0,±1,±1)/2.

    Parameters: r1=air sphere radius, r2=dielectric sphere radius, r3=rod radius.
    Optimal: (0.12, 0.19, 0.08) giving ~18.3% gap (Fig. 7b).
    """
    a1 = 0.5 * np.array([0, 1, 1], dtype=float)
    a2 = 0.5 * np.array([1, 0, 1], dtype=float)
    a3 = 0.5 * np.array([1, 1, 0], dtype=float)
    A = np.column_stack([a1, a2, a3])
    A_inv = np.linalg.inv(A)

    atom = np.array([0.0, 0.0, 0.0])

    # 12 nearest neighbor half-bond vectors (we place half-bonds from atom)
    nn_half = []
    for i in range(3):
        for s1 in [1, -1]:
            for s2 in [1, -1]:
                v = np.zeros(3)
                v[(i + 1) % 3] = s1 * 0.5
                v[(i + 2) % 3] = s2 * 0.5
                nn_half.append(v)
    nn_half = np.array(nn_half)  # (12, 3)

    def eps_func(s1, s2, s3):
        shape = s1.shape
        eps = np.full(shape, eps_lo, dtype=float)

        # Convert fractional to Cartesian
        x = s1 * a1[0] + s2 * a2[0] + s3 * a3[0]
        y = s1 * a1[1] + s2 * a2[1] + s3 * a3[1]
        z = s1 * a1[2] + s2 * a2[2] + s3 * a3[2]
        pos = np.stack([x, y, z], axis=-1)

        # Minimum image distance to atom at origin
        dp = _min_image_disp(pos, atom, A, A_inv)
        dist = np.linalg.norm(dp, axis=-1)

        # Dielectric sphere
        eps = np.where(dist < r2, eps_hi, eps)

        # Connecting half-rods
        for nn in nn_half:
            perp_dist, proj, blen = _point_to_cylinder_dist(
                pos, atom, nn, A, A_inv)
            # Half bond: 0 to bond_length/2
            in_rod = (proj >= 0) & (proj <= blen) & (perp_dist < r3)
            eps = np.where(in_rod, eps_hi, eps)

        # Air hole
        eps = np.where(dist < r1, eps_lo, eps)

        return eps

    return eps_func


# ==================== Brillouin zone k-paths ====================

def sc_kpath(b1, b2, b3):
    """SC BZ high-symmetry path: Gamma-X-M-Gamma-R-X"""
    G = np.zeros(3)
    X = 0.5 * b1
    M = 0.5 * (b1 + b2)
    R = 0.5 * (b1 + b2 + b3)

    path = [G, X, M, G, R, X]
    labels = [r'$\Gamma$', 'X', 'M', r'$\Gamma$', 'R', 'X']
    return path, labels


def fcc_kpath(b1, b2, b3):
    """
    FCC BZ high-symmetry path: Gamma-X-U|K-Gamma-L-W-X

    Points expressed in Cartesian coordinates (same basis as G vectors).
    For FCC with conventional a=1:
      b1 = 2pi(-1,1,1), b2 = 2pi(1,-1,1), b3 = 2pi(1,1,-1)
    """
    twopi = 2 * np.pi
    G = np.zeros(3)
    # X = center of square face of BZ
    X = twopi * np.array([0, 1, 0])  # (b1+b2)/2
    # L = center of hexagonal face
    L = twopi * np.array([0.5, 0.5, 0.5])  # (b1+b2+b3)/2... wait

    # Standard FCC BZ high-symmetry points in Cartesian (units of 2pi/a):
    # Gamma: (0,0,0)
    # X: (0,0,1) → in our coords, this is (b1+b2)/2
    # L: (1/2,1/2,1/2) → (b1+b2+b3)/4... hmm

    # Let me use the standard definition.
    # For FCC with b1=2pi(-1,1,1), b2=2pi(1,-1,1), b3=2pi(1,1,-1):
    # X = (b3+b1)/2 = 2pi(0,2,0)/2 = 2pi(0,1,0)?
    # Wait: (b3+b1)/2 = (2pi(1,1,-1) + 2pi(-1,1,1))/2 = 2pi(0,2,0)/2 = 2pi(0,1,0)
    # Hmm that's correct. And |X| = 2pi, which is correct.

    # L = (b1+b2+b3)/2 = 2pi(-1+1+1, 1-1+1, 1+1-1)/2 = 2pi(1/2, 1/2, 1/2)
    # |L| = 2pi*sqrt(3)/2 ≈ 5.44. But L should be at half the body diagonal.

    # Actually for FCC, L is at the center of a hexagonal face:
    # L = b1/2 = 2pi(-1/2, 1/2, 1/2) in Cartesian
    # |L| = 2pi*sqrt(3)/2 ≈ 5.44.

    # Let me use the correct standard points:
    # X = (b1+b3)/2 = pi*(0, 2, 0) = (0, 2pi, 0), or equivalently (2pi, 0, 0)
    # by symmetry, all permutations are equivalent. Let me use X=(0,0,2pi).

    # Standard convention: express in units of 2pi/a
    X = twopi * np.array([0, 1, 0])
    L = 0.5 * b1  # = pi*(-1,1,1)
    W = twopi * np.array([0.5, 1, 0])
    K = twopi * np.array([0.75, 0.75, 0])
    U = twopi * np.array([0.25, 1, 0.25])

    path = [G, X, U, K, G, L, W, X]
    labels = [r'$\Gamma$', 'X', 'U|K', 'K', r'$\Gamma$', 'L', 'W', 'X']
    return path, labels


def sc_ibz_kpoints(b1, b2, b3, n_per_edge=8):
    """Sample k-points on edges of SC irreducible Brillouin zone."""
    G = np.zeros(3)
    X = 0.5 * b1
    M = 0.5 * (b1 + b2)
    R = 0.5 * (b1 + b2 + b3)

    edges = [(G, X), (X, M), (M, G), (G, R), (R, X), (M, R)]
    kpts = []
    for start, end in edges:
        for i in range(n_per_edge):
            t = i / n_per_edge
            kpts.append((1 - t) * start + t * end)
    kpts.append(R)
    return np.array(kpts)


def fcc_ibz_kpoints(b1, b2, b3, n_per_edge=8):
    """Sample k-points on edges of FCC irreducible Brillouin zone."""
    twopi = 2 * np.pi
    G = np.zeros(3)
    X = twopi * np.array([0, 1, 0])
    L = 0.5 * b1
    W = twopi * np.array([0.5, 1, 0])
    K = twopi * np.array([0.75, 0.75, 0])
    U = twopi * np.array([0.25, 1, 0.25])

    edges = [(G, X), (X, W), (W, L), (L, G), (G, K), (K, X), (U, L), (W, K)]
    kpts = []
    for start, end in edges:
        for i in range(n_per_edge):
            t = i / n_per_edge
            kpts.append((1 - t) * start + t * end)
    kpts.append(K)
    return np.array(kpts)
