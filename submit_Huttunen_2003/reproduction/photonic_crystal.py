"""
Core module for photonic crystal slab computations.

Implements:
- Transfer matrix method for 1D photonic crystal band structure
- Transfer matrix reflectance for finite PhC stacks
- 2D FDTD for TE polarization in 1D photonic crystal slab
"""

import numpy as np


# ===========================================================================
# Transfer Matrix Method
# ===========================================================================

def layer_transfer_matrix(k0, eps, d):
    """
    Transfer matrix for a single homogeneous layer (normal incidence, TE).

    Parameters
    ----------
    k0 : complex
        Free-space wave number (omega/c)
    eps : complex
        Dielectric constant
    d : float
        Layer thickness

    Returns
    -------
    M : ndarray (2x2, complex)
    """
    n = np.sqrt(eps + 0j)
    kz = k0 * n
    if abs(kz) < 1e-15:
        return np.eye(2, dtype=complex)
    phi = kz * d
    cos_phi = np.cos(phi)
    sin_phi = np.sin(phi)
    eta = n  # impedance factor for TE
    M = np.array([
        [cos_phi, 1j * sin_phi / eta],
        [1j * eta * sin_phi, cos_phi]
    ], dtype=complex)
    return M


def unit_cell_matrix(k0, eps1, eps2, d1, d2):
    """Transfer matrix for one period of a 1D PhC."""
    M1 = layer_transfer_matrix(k0, eps1, d1)
    M2 = layer_transfer_matrix(k0, eps2, d2)
    return M2 @ M1


def band_structure_trace(omega_norm, eps1, eps2, d1, d2):
    """
    Compute Tr(M)/2 for normalized frequency omega_norm = omega*P/(2*pi*c).
    P = d1 + d2 = 1 in our units.
    """
    k0 = omega_norm * 2.0 * np.pi  # since c=1, P=1
    M = unit_cell_matrix(k0, eps1, eps2, d1, d2)
    return 0.5 * np.real(M[0, 0] + M[1, 1])


def is_in_band(omega_norm, eps1, eps2, d1, d2):
    """Check if normalized frequency is in a propagating band."""
    half_trace = band_structure_trace(omega_norm, eps1, eps2, d1, d2)
    return abs(half_trace) <= 1.0


def find_bandgaps(eps1, eps2, d1, d2, n_omega=2000, omega_max=0.8):
    """Find band gaps of 1D photonic crystal at normal incidence."""
    omegas = np.linspace(0.001, omega_max, n_omega)
    in_band = np.array([is_in_band(w, eps1, eps2, d1, d2) for w in omegas])

    gaps = []
    in_gap = False
    gap_start = 0
    for i in range(n_omega):
        if not in_band[i] and not in_gap:
            in_gap = True
            gap_start = omegas[i]
        elif in_band[i] and in_gap:
            in_gap = False
            gaps.append((gap_start, omegas[i]))
    return gaps


def finite_stack_reflectance(omega_norm, eps1, eps2, d1, d2, n_periods, eps_in=1.0, eps_out=1.0):
    """
    Reflectance of a finite 1D PhC stack embedded between two semi-infinite media.

    Parameters
    ----------
    omega_norm : float
        Normalized frequency omega*P/(2*pi*c)
    eps1, eps2 : float
        Dielectric constants of the two layers
    d1, d2 : float
        Layer thicknesses (fractions of period)
    n_periods : int
        Number of unit cells
    eps_in, eps_out : float
        Dielectric constants of input/output media

    Returns
    -------
    R : float
        Power reflectance |r|^2
    """
    k0 = omega_norm * 2.0 * np.pi
    M_cell = unit_cell_matrix(k0, eps1, eps2, d1, d2)

    # Total matrix = M_cell^n_periods
    M_total = np.eye(2, dtype=complex)
    for _ in range(n_periods):
        M_total = M_cell @ M_total

    # Compute reflection coefficient
    n_in = np.sqrt(eps_in + 0j)
    n_out = np.sqrt(eps_out + 0j)

    # r = (M11*n_out + M12*n_in*n_out - M21 - M22*n_in) /
    #     (M11*n_out + M12*n_in*n_out + M21 + M22*n_in)
    a = M_total[0, 0] * n_out + M_total[0, 1] * n_in * n_out
    b = M_total[1, 0] + M_total[1, 1] * n_in
    r = (a - b) / (a + b)
    return float(np.abs(r) ** 2)


def slab_guided_mode_condition(omega_norm, eps_slab_eff, eps_boundary, h):
    """
    Check if a guided mode exists in a slab waveguide.

    For a symmetric slab waveguide with effective index n_eff = sqrt(eps_slab_eff)
    and boundary index n_b = sqrt(eps_boundary), guided modes exist when
    n_b < n_eff (total internal reflection).

    Parameters
    ----------
    omega_norm : float
        Normalized frequency
    eps_slab_eff : float
        Effective dielectric constant of the slab
    eps_boundary : float
        Dielectric constant of boundary material
    h : float
        Slab height (normalized to period)

    Returns
    -------
    has_mode : bool
    n_guided : float
        Effective guided mode index (0 if no mode)
    """
    n_slab = np.sqrt(eps_slab_eff)
    n_bound = np.sqrt(eps_boundary)

    if n_slab <= n_bound:
        return False, 0.0

    # For the lowest TE mode: kz*h = 2*arctan(gamma/kz)
    # where kz = k0*sqrt(n_slab^2 - n_eff^2), gamma = k0*sqrt(n_eff^2 - n_bound^2)
    k0 = omega_norm * 2.0 * np.pi

    # Search for guided mode index between n_bound and n_slab
    from scipy.optimize import brentq

    def mode_eq(n_eff):
        if n_eff <= n_bound or n_eff >= n_slab:
            return 1.0
        kz = k0 * np.sqrt(n_slab ** 2 - n_eff ** 2)
        gamma = k0 * np.sqrt(n_eff ** 2 - n_bound ** 2)
        return np.tan(kz * h / 2.0) - gamma / kz

    try:
        n_eff = brentq(mode_eq, n_bound + 0.001, n_slab - 0.001)
        return True, n_eff
    except (ValueError, RuntimeError):
        return False, 0.0


# ===========================================================================
# 2D FDTD for 1D Photonic Crystal Slab (TE: Hx, Ey, Ez)
# ===========================================================================

class FDTD2D_TE:
    """
    2D FDTD for TE polarization (Hx, Ey, Ez) in the y-z plane.
    Uses CPML absorbing boundary conditions.
    """

    def __init__(self, ny, nz, dy, dz, dt, eps_grid, pml_thick=10):
        self.ny = ny
        self.nz = nz
        self.dy = dy
        self.dz = dz
        self.dt = dt
        self.eps = eps_grid.astype(np.float64)
        self.pml_thick = pml_thick

        # Field arrays
        self.Hx = np.zeros((ny, nz))
        self.Ey = np.zeros((ny, nz))
        self.Ez = np.zeros((ny, nz))

        # PML conductivity profiles
        self.sigma_y = np.zeros(ny)
        self.sigma_z = np.zeros(nz)
        sigma_max = 0.8 * (pml_thick + 1) / (dy * np.sqrt(1.0))  # empirical

        for i in range(pml_thick):
            sigma = sigma_max * ((pml_thick - i) / pml_thick) ** 3
            self.sigma_y[i] = sigma
            self.sigma_y[ny - 1 - i] = sigma
            self.sigma_z[i] = sigma
            self.sigma_z[nz - 1 - i] = sigma

        # PML coefficients
        self.cy1 = np.exp(-self.sigma_y * dt)
        self.cy2 = np.where(self.sigma_y > 0,
                            (1.0 - self.cy1) / (self.sigma_y * dy),
                            dt / dy * np.ones(ny))
        self.cz1 = np.exp(-self.sigma_z * dt)
        self.cz2 = np.where(self.sigma_z > 0,
                            (1.0 - self.cz1) / (self.sigma_z * dz),
                            dt / dz * np.ones(nz))

        # PML auxiliary fields
        self.psi_hx_y = np.zeros((ny, nz))
        self.psi_hx_z = np.zeros((ny, nz))
        self.psi_ey_z = np.zeros((ny, nz))
        self.psi_ez_y = np.zeros((ny, nz))

    def step(self):
        dt = self.dt
        dy = self.dy
        dz = self.dz

        # Update Hx
        # Hx += dt * (dEy/dz - dEz/dy)
        dEy_dz = np.zeros_like(self.Hx)
        dEz_dy = np.zeros_like(self.Hx)
        dEy_dz[:, :-1] = (self.Ey[:, 1:] - self.Ey[:, :-1]) / dz
        dEz_dy[:-1, :] = (self.Ez[1:, :] - self.Ez[:-1, :]) / dy

        # PML update for H
        self.psi_hx_z[:, :-1] = self.cz1[:-1] * self.psi_hx_z[:, :-1] + \
            self.cz2[:-1] * (self.Ey[:, 1:] - self.Ey[:, :-1])
        self.psi_hx_y[:-1, :] = self.cy1[:-1, np.newaxis] * self.psi_hx_y[:-1, :] + \
            (self.cy2[:-1, np.newaxis]) * (self.Ez[1:, :] - self.Ez[:-1, :])

        self.Hx += dt * (dEy_dz - dEz_dy)

        # Update Ey: dEy/dt = (1/eps) * dHx/dz
        dHx_dz = np.zeros_like(self.Ey)
        dHx_dz[:, 1:] = (self.Hx[:, 1:] - self.Hx[:, :-1]) / dz

        self.psi_ey_z[:, 1:] = self.cz1[1:] * self.psi_ey_z[:, 1:] + \
            self.cz2[1:] * (self.Hx[:, 1:] - self.Hx[:, :-1])

        self.Ey += dt / self.eps * dHx_dz

        # Update Ez: dEz/dt = -(1/eps) * dHx/dy
        dHx_dy = np.zeros_like(self.Ez)
        dHx_dy[1:, :] = (self.Hx[1:, :] - self.Hx[:-1, :]) / dy

        self.psi_ez_y[1:, :] = self.cy1[1:, np.newaxis] * self.psi_ez_y[1:, :] + \
            (self.cy2[1:, np.newaxis]) * (self.Hx[1:, :] - self.Hx[:-1, :])

        self.Ez -= dt / self.eps * dHx_dy


def run_fdtd_reflection(eps_boundary, eps_rod=13.0, l_frac=0.2, h_frac=0.5,
                         omega_norm=0.3, resolution=20, n_periods=30,
                         n_steps=2000, return_snapshots=False):
    """
    Run 2D FDTD simulation for 1D PhC slab with boundary material change.

    Returns reflected energy fraction.
    """
    dy = 1.0 / resolution
    dz = 1.0 / resolution

    ny = n_periods * resolution + 200
    nz = int(3.0 / dz)

    # Courant
    dt = 0.4 / np.sqrt(1.0 / dy ** 2 + 1.0 / dz ** 2)

    # Slab geometry
    slab_h = int(h_frac / dz)
    slab_z_lo = nz // 2 - slab_h // 2
    slab_z_hi = slab_z_lo + slab_h

    boundary_y = ny // 2

    # Build epsilon
    eps_grid = np.ones((ny, nz))
    eps_grid[boundary_y:, :slab_z_lo] = eps_boundary
    eps_grid[boundary_y:, slab_z_hi:] = eps_boundary

    period_cells = resolution
    rod_cells = max(int(l_frac * resolution), 1)
    y_start = 100

    for p in range(n_periods):
        y0 = y_start + p * period_cells
        if y0 + rod_cells <= ny:
            eps_grid[y0:y0 + rod_cells, slab_z_lo:slab_z_hi] = eps_rod

    fdtd = FDTD2D_TE(ny, nz, dy, dz, dt, eps_grid)

    source_y = boundary_y - 8 * resolution
    detect_y = source_y - 3 * resolution

    omega = omega_norm * 2.0 * np.pi
    sigma_t = 12.0
    t_center = 35.0

    reflected = 0.0
    incident = 0.0
    snapshots = []
    pulse_done_step = int((t_center + 3 * sigma_t) / dt)

    for n in range(n_steps):
        t = n * dt

        # Gaussian pulse
        envelope = np.exp(-((t - t_center) ** 2) / (2 * sigma_t ** 2))
        if envelope > 1e-8:
            signal = 0.1 * np.sin(omega * t) * envelope
            fdtd.Ey[source_y, slab_z_lo:slab_z_hi] += signal

        fdtd.step()

        if n > pulse_done_step:
            reflected += np.sum(fdtd.Ey[detect_y, slab_z_lo:slab_z_hi] ** 2) * dt
        incident += np.sum(fdtd.Ey[source_y + 2, slab_z_lo:slab_z_hi] ** 2) * dt

        if return_snapshots and n in [int(n_steps * f) for f in [0.25, 0.5, 0.75, 0.95]]:
            z_mid = (slab_z_lo + slab_z_hi) // 2
            snapshots.append((t, fdtd.Ey[:, z_mid].copy(), boundary_y * dy))

    frac = reflected / max(incident, 1e-30)
    if return_snapshots:
        return frac, snapshots
    return frac
