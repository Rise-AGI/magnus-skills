"""
Electrostatic Energy-Conserving Particle-in-Cell (EC-PIC) code.

Translated from the Matlab/Octave code in Appendices A and B of:
  Markidis & Lapenta, "The Energy Conserving Particle-in-Cell Method",
  J. Comput. Phys. 230 (2011) 7037-7052.

Uses midpoint implicit time integration solved by scipy's Newton-Krylov solver.
"""

import numpy as np
from scipy.optimize import newton_krylov
from scipy.sparse import csr_matrix


def interpolate_to_grid(x, dx, NG):
    """Cloud-in-Cell interpolation: particle positions -> grid weights."""
    g1 = np.floor(x / dx - 0.5).astype(int) + 1
    fraz1 = 1.0 - np.abs(x / dx - g1 + 0.5)
    g = np.stack([g1, g1 + 1])
    fraz = np.stack([fraz1, 1.0 - fraz1])
    g[g < 1] += NG
    g[g > NG] -= NG
    # Convert to 0-indexed
    g = g - 1
    return g, fraz


def build_interp_matrix(x, dx, NG, N):
    """Build sparse interpolation matrix (N x NG)."""
    g, fraz = interpolate_to_grid(x, dx, NG)
    rows = np.concatenate([np.arange(N), np.arange(N)])
    cols = np.concatenate([g[0], g[1]])
    vals = np.concatenate([fraz[0], fraz[1]])
    mat = csr_matrix((vals, (rows, cols)), shape=(N, NG))
    return mat


def compute_residual(xkrylov, x0, v0, E0, Q, QM, dx, DT, NG, N, rho_back):
    """Compute residual for the EC-PIC system (electrostatic).
    xkrylov = [v_average(N), E_new(NG)]
    """
    v_avg = xkrylov[:N]
    E_new = xkrylov[N:N + NG]

    # Average position
    x_avg = x0 + v_avg * DT / 2.0
    L = dx * NG
    x_avg = x_avg % L

    # Interpolation matrix at average position
    mat = build_interp_matrix(x_avg, dx, NG, N)

    res = np.zeros(N + NG)

    # Velocity residual: v_avg - v0 - 0.25*QM*DT * W*(E0 + E_new) = 0
    E_avg = 0.5 * (E0 + E_new)
    E_at_particle = mat.dot(E_avg)
    res[:N] = v_avg - v0 - 0.5 * QM * DT * E_at_particle

    # Current density
    g, fraz = interpolate_to_grid(x_avg, dx, NG)
    rows = np.concatenate([np.arange(N), np.arange(N)])
    cols = np.concatenate([g[0], g[1]])
    vals_j = np.concatenate([fraz[0] * v_avg, fraz[1] * v_avg])
    J = np.zeros(NG)
    np.add.at(J, cols, vals_j * Q / dx)

    # Field residual: E_new - E0 + J*DT = 0
    res[N:N + NG] = E_new - E0 + J * DT

    return res


def solve_poisson(x0, Q, dx, NG, N, rho_back):
    """Solve Poisson equation to get initial E field satisfying Gauss' law."""
    mat = build_interp_matrix(x0, dx, NG, N)
    rho = (Q / dx) * np.array(mat.sum(axis=0)).flatten() + rho_back

    # Poisson matrix (periodic BC via spectral method)
    from numpy.fft import fft, ifft
    k = 2.0 * np.pi * np.fft.fftfreq(NG, d=dx)
    rho_hat = fft(rho)
    phi_hat = np.zeros(NG, dtype=complex)
    phi_hat[1:] = rho_hat[1:] / (k[1:] ** 2)  # Gauss law: -k^2 phi = -rho (CGS: 4pi*rho)
    # Actually in CGS: nabla^2 phi = -4*pi*rho, so phi_hat = 4*pi*rho_hat / k^2
    phi_hat[1:] = 4.0 * np.pi * rho_hat[1:] / (k[1:] ** 2)
    phi = np.real(ifft(phi_hat))

    # E = -grad(phi) using centered differences (periodic)
    E = np.zeros(NG)
    E = -(np.roll(phi, -1) - np.roll(phi, 1)) / (2.0 * dx)
    return E


def run_electrostatic_pic(L, DT, NT, NG, N, WP, QM, x0, v0, tol=1e-7,
                           save_phase_space_at=None, save_field_at=None):
    """Run electrostatic energy-conserving PIC simulation.

    Parameters
    ----------
    L : float - simulation box length
    DT : float - time step
    NT : int - number of time steps
    NG : int - number of grid cells
    N : int - number of particles
    WP : float - plasma frequency
    QM : float - charge-to-mass ratio
    x0 : array - initial particle positions
    v0 : array - initial particle velocities
    tol : float - solver tolerance
    save_phase_space_at : list of int - time steps to save phase space
    save_field_at : list of int - time steps to save field data

    Returns
    -------
    dict with energy history, saved phase spaces, etc.
    """
    dx = L / NG
    Q = WP ** 2 / (QM * N / L)
    rho_back = -Q * N / L  # background ion charge density (neutralizing)

    # Initial electric field from Poisson
    E0 = solve_poisson(x0, Q, dx, NG, N, rho_back)

    energy_history = []
    phase_spaces = {}
    field_data = {}

    for it in range(NT):
        # Total energy
        KE = 0.5 * np.abs(Q) * np.sum(v0 ** 2)
        FE = 0.5 * np.sum(E0 ** 2) * dx
        Etot = KE + FE
        energy_history.append(Etot)

        if save_phase_space_at and it in save_phase_space_at:
            phase_spaces[it] = (x0.copy(), v0.copy())
        if save_field_at and it in save_field_at:
            field_data[it] = E0.copy()

        # Initial guess for solver
        xkrylov0 = np.concatenate([v0, E0])

        # Define residual for this time step
        def res_func(xk):
            return compute_residual(xk, x0, v0, E0, Q, QM, dx, DT, NG, N, rho_back)

        try:
            sol = newton_krylov(res_func, xkrylov0,
                                f_tol=tol * np.sqrt(N + NG),
                                maxiter=40, verbose=False)
        except Exception:
            # If Newton-Krylov fails, try with relaxed tolerance
            try:
                sol = newton_krylov(res_func, xkrylov0,
                                    f_tol=tol * 10 * np.sqrt(N + NG),
                                    maxiter=60, verbose=False)
            except Exception:
                print(f"Warning: solver failed at step {it}, using last guess")
                sol = xkrylov0

        v_avg = sol[:N]
        E_new = sol[N:N + NG]

        # Update particles
        v0 = 2.0 * v_avg - v0
        x0 = x0 + v_avg * DT
        x0 = x0 % L

        # Update field
        E0 = E_new

    # Final energy
    KE = 0.5 * np.abs(Q) * np.sum(v0 ** 2)
    FE = 0.5 * np.sum(E0 ** 2) * dx
    energy_history.append(KE + FE)

    return {
        'energy': np.array(energy_history),
        'phase_spaces': phase_spaces,
        'field_data': field_data,
        'x_final': x0,
        'v_final': v0,
        'E_final': E0,
    }


def run_explicit_pic(L, DT, NT, NG, N, WP, QM, x0, v0):
    """Run explicit momentum-conserving leapfrog PIC (for comparison)."""
    dx = L / NG
    Q = WP ** 2 / (QM * N / L)
    rho_back = -Q * N / L

    E0 = solve_poisson(x0, Q, dx, NG, N, rho_back)

    energy_history = []

    # Half-step velocity push for leapfrog initialization
    mat = build_interp_matrix(x0, dx, NG, N)
    E_at_p = mat.dot(E0)
    v_half = v0 + 0.5 * QM * DT * E_at_p

    for it in range(NT):
        # Energy (at integer time)
        v_int = v_half - 0.5 * QM * DT * E_at_p
        KE = 0.5 * np.abs(Q) * np.sum(v_int ** 2)
        FE = 0.5 * np.sum(E0 ** 2) * dx
        energy_history.append(KE + FE)

        # Push positions
        x0 = (x0 + v_half * DT) % L

        # Deposit charge
        mat = build_interp_matrix(x0, dx, NG, N)
        rho = (Q / dx) * np.array(mat.sum(axis=0)).flatten() + rho_back

        # Solve Poisson for E
        from numpy.fft import fft, ifft
        k = 2.0 * np.pi * np.fft.fftfreq(NG, d=dx)
        rho_hat = fft(rho)
        phi_hat = np.zeros(NG, dtype=complex)
        phi_hat[1:] = 4.0 * np.pi * rho_hat[1:] / (k[1:] ** 2)
        phi = np.real(ifft(phi_hat))
        E0 = -(np.roll(phi, -1) - np.roll(phi, 1)) / (2.0 * dx)

        # Push velocity
        E_at_p = mat.dot(E0)
        v_half = v_half + QM * DT * E_at_p

    # Final energy
    v_int = v_half - 0.5 * QM * DT * E_at_p
    KE = 0.5 * np.abs(Q) * np.sum(v_int ** 2)
    FE = 0.5 * np.sum(E0 ** 2) * dx
    energy_history.append(KE + FE)

    return {
        'energy': np.array(energy_history),
        'x_final': x0,
        'v_final': v_half,
    }
