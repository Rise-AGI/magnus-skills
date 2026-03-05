"""
Efficient 3D Plane-Wave Expansion (PWE) solver for photonic crystal band structures.

Men et al., Optics Express 22(19), 22632-22657 (2014)
"Robust topology optimization of three-dimensional photonic-crystal band-gap structures"

Solves: curl(1/eps(r) * curl H) = (omega/c)^2 H
in a plane-wave basis for periodic dielectric structures.
"""

import numpy as np
from scipy.linalg import eigh


def reciprocal_vectors(a1, a2, a3):
    """Compute reciprocal lattice vectors b_i such that a_i . b_j = 2*pi*delta_ij."""
    vol = np.dot(a1, np.cross(a2, a3))
    b1 = 2 * np.pi * np.cross(a2, a3) / vol
    b2 = 2 * np.pi * np.cross(a3, a1) / vol
    b3 = 2 * np.pi * np.cross(a1, a2) / vol
    return b1, b2, b3


class PWESolver3D:
    """
    3D Plane-Wave Expansion solver for photonic band structures.

    Parameters
    ----------
    a1, a2, a3 : lattice vectors (in units of a)
    eps_func : callable(x, y, z) -> epsilon, where x,y,z are fractional coordinates
    n_max : maximum Miller index for plane waves (total: (2*n_max+1)^3)
    grid_res : real-space grid resolution for FFT of 1/epsilon
    """

    def __init__(self, a1, a2, a3, eps_func, n_max=2, grid_res=24,
                 use_inverse_method=True):
        self.a1 = np.array(a1, dtype=float)
        self.a2 = np.array(a2, dtype=float)
        self.a3 = np.array(a3, dtype=float)
        self.b1, self.b2, self.b3 = reciprocal_vectors(a1, a2, a3)
        self.n_max = n_max
        self.grid_res = grid_res

        # Generate G vectors (Miller indices)
        ns = np.arange(-n_max, n_max + 1)
        n1, n2, n3 = np.meshgrid(ns, ns, ns, indexing='ij')
        self.miller = np.column_stack([n1.ravel(), n2.ravel(), n3.ravel()])
        self.n_g = len(self.miller)

        # G vectors in Cartesian coordinates
        self.G = (self.miller[:, 0:1] * self.b1 +
                  self.miller[:, 1:2] * self.b2 +
                  self.miller[:, 2:3] * self.b3)

        # Compute eta matrix using inverse method (Ho-Chan-Soukoulis)
        # or direct method
        self._compute_eta(eps_func, use_inverse_method)

    def _compute_eta(self, eps_func, use_inverse_method):
        """
        Compute eta matrix = Fourier representation of 1/epsilon.

        If use_inverse_method=True (recommended for high index contrast):
          eta = inverse of the matrix eps_{GG'} = FT(eps)(G-G')
        Otherwise (direct method):
          eta_{GG'} = FT(1/eps)(G-G')
        """
        N = self.grid_res
        s = np.linspace(0, 1, N, endpoint=False)
        s1, s2, s3 = np.meshgrid(s, s, s, indexing='ij')
        eps_grid = eps_func(s1, s2, s3)

        if use_inverse_method:
            # Inverse method: eta = [eps_matrix]^{-1}
            eps_fft = np.fft.fftn(eps_grid) / (N * N * N)
            dm = self.miller[:, np.newaxis, :] - self.miller[np.newaxis, :, :]
            eps_mat = eps_fft[dm[:, :, 0] % N, dm[:, :, 1] % N, dm[:, :, 2] % N]
            # Invert the epsilon matrix
            self.eta_mat = np.linalg.inv(eps_mat)
        else:
            # Direct method: eta = FT(1/eps)
            inv_eps_fft = np.fft.fftn(1.0 / eps_grid) / (N * N * N)
            dm = self.miller[:, np.newaxis, :] - self.miller[np.newaxis, :, :]
            self.eta_mat = inv_eps_fft[dm[:, :, 0] % N, dm[:, :, 1] % N, dm[:, :, 2] % N]

    def _build_eta_matrix(self):
        """Return pre-computed eta matrix."""
        return self.eta_mat

    def solve_k(self, k_vec, n_bands=10):
        """
        Solve for eigenfrequencies at a single k-point.

        Parameters
        ----------
        k_vec : Bloch wavevector in Cartesian coordinates (units of 2*pi/a are implicit)
        n_bands : number of bands to return

        Returns
        -------
        freqs : normalized frequencies omega*a/(2*pi*c), sorted ascending
        """
        k = np.array(k_vec, dtype=float)

        # k+G vectors
        kpG = k + self.G  # shape (n_g, 3)
        kpG_norm = np.linalg.norm(kpG, axis=1)

        # Compute polarization vectors for each k+G
        e1, e2 = self._compute_polarizations(kpG, kpG_norm)

        # Build Maxwell matrix using vectorized operations
        eta = self._build_eta_matrix()

        # Cross products: (k+Gi) x e_ia for a=1,2
        # cross_1[i] = kpG[i] x e1[i], shape (n_g, 3)
        cross_1 = np.cross(kpG, e1)
        cross_2 = np.cross(kpG, e2)

        # Matrix elements:
        # M[i*2+a, j*2+b] = eta[i,j] * cross_a[i] . cross_b[j]
        # We compute four blocks: (1,1), (1,2), (2,1), (2,2)

        # dot products: cross_a[i] . cross_b[j] = sum_k cross_a[i,k] * cross_b[j,k]
        # Shape: (n_g, n_g) for each block
        d11 = np.einsum('ik,jk->ij', cross_1, cross_1)
        d12 = np.einsum('ik,jk->ij', cross_1, cross_2)
        d21 = np.einsum('ik,jk->ij', cross_2, cross_1)
        d22 = np.einsum('ik,jk->ij', cross_2, cross_2)

        # Full matrix
        ng = self.n_g
        M = np.zeros((2 * ng, 2 * ng), dtype=complex)
        M[:ng, :ng] = eta * d11
        M[:ng, ng:] = eta * d12
        M[ng:, :ng] = eta * d21
        M[ng:, ng:] = eta * d22

        # Symmetrize
        M = 0.5 * (M + M.conj().T)

        # Solve eigenvalue problem (must use full complex Hermitian matrix)
        n_eig = min(2 * ng, n_bands + 2)
        try:
            eigenvalues = eigh(M, eigvals_only=True,
                               subset_by_index=[0, n_eig - 1])
        except Exception:
            eigenvalues = eigh(M, eigvals_only=True)

        # Convert to normalized frequency omega*a/(2*pi*c)
        freqs = np.sqrt(np.maximum(eigenvalues, 0)) / (2 * np.pi)
        return np.sort(freqs)[:n_bands]

    def _compute_polarizations(self, kpG, kpG_norm):
        """Compute two orthogonal unit vectors perpendicular to each k+G."""
        ng = len(kpG)
        e1 = np.zeros((ng, 3))
        e2 = np.zeros((ng, 3))

        for i in range(ng):
            if kpG_norm[i] < 1e-12:
                e1[i] = [1, 0, 0]
                e2[i] = [0, 1, 0]
            else:
                v = kpG[i] / kpG_norm[i]
                t = np.array([1, 0, 0]) if abs(v[0]) < 0.9 else np.array([0, 1, 0])
                e1[i] = np.cross(v, t)
                e1[i] /= np.linalg.norm(e1[i])
                e2[i] = np.cross(v, e1[i])
                e2[i] /= np.linalg.norm(e2[i])

        return e1, e2

    def compute_bands(self, k_path_cart, n_bands=10, verbose=False):
        """
        Compute band structure along a k-path.

        Parameters
        ----------
        k_path_cart : list of k-points in Cartesian coords
        n_bands : number of bands
        verbose : print progress

        Returns
        -------
        all_freqs : array shape (n_kpoints, n_bands)
        """
        n_k = len(k_path_cart)
        all_freqs = np.zeros((n_k, n_bands))

        for i, k in enumerate(k_path_cart):
            if verbose and i % 10 == 0:
                print(f"  k-point {i+1}/{n_k}")
            freqs = self.solve_k(k, n_bands)
            n_got = min(len(freqs), n_bands)
            all_freqs[i, :n_got] = freqs[:n_got]
            if n_got < n_bands:
                all_freqs[i, n_got:] = freqs[-1] if len(freqs) > 0 else 0

        return all_freqs


def make_k_path(k_points, n_per_segment=20):
    """
    Generate k-path with uniform sampling between high-symmetry points.

    Parameters
    ----------
    k_points : list of high-symmetry k-points (Cartesian)
    n_per_segment : points per segment

    Returns
    -------
    k_list : list of k-vectors
    k_dist : cumulative distances
    tick_pos : positions of high-symmetry points
    """
    k_list = []
    k_dist = []
    tick_pos = [0.0]
    dist = 0.0

    for seg in range(len(k_points) - 1):
        start = np.array(k_points[seg])
        end = np.array(k_points[seg + 1])
        n_pts = n_per_segment

        for j in range(n_pts):
            t = j / n_pts
            k = (1 - t) * start + t * end
            if len(k_list) > 0:
                dist += np.linalg.norm(k - k_list[-1])
            k_list.append(k)
            k_dist.append(dist)

        if seg == len(k_points) - 2:
            k_list.append(end)
            dist += np.linalg.norm(end - k_list[-2])
            k_dist.append(dist)

        tick_pos.append(dist if seg < len(k_points) - 2 else k_dist[-1])

    return k_list, np.array(k_dist), tick_pos


def compute_gap(all_freqs, band_below, band_above):
    """
    Compute fractional band gap.

    Parameters
    ----------
    all_freqs : shape (n_k, n_bands), normalized frequencies
    band_below : 0-based index of lower band
    band_above : 0-based index of upper band

    Returns
    -------
    gap_pct : fractional gap in percent
    """
    max_below = np.max(all_freqs[:, band_below])
    min_above = np.min(all_freqs[:, band_above])
    mid = 0.5 * (max_below + min_above)
    if mid < 1e-10:
        return 0.0
    return (min_above - max_below) / mid * 100.0
