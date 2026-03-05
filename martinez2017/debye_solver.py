"""
Core solver for the Debye shielding comparison.

Solves the normalized Poisson equation (Eq. 3.3 of Martinez-Fuentes & Herrera-Velazquez 2017):
    d^2[rho * Psi(rho)] / drho^2 = -N_D * rho * [1 - exp(Psi(rho) / N_D)]

using 4th-order Runge-Kutta, integrating backward from large rho to small rho.
Compares with the approximate Yukawa solution: Psi_a(rho) = exp(-rho) / rho.
"""

import numpy as np


def psi_approx(rho):
    """Approximate (Yukawa) solution: Psi_a(rho) = exp(-rho) / rho."""
    return np.exp(-rho) / rho


def rhs(rho, u, du, N_D):
    """
    Right-hand side for the ODE system.

    Let w(rho) = rho * Psi(rho), then:
        w'' = -N_D * rho * [1 - exp(w / (rho * N_D))]

    Variables:
        u  = w = rho * Psi(rho)
        du = w' = d(rho * Psi) / drho
    """
    if rho < 1e-15:
        return 0.0
    psi = u / rho
    arg = psi / N_D
    if arg > 500:
        arg = 500.0
    return -N_D * rho * (1.0 - np.exp(arg))


def solve_exact(N_D, rho_start, rho_end, n_steps):
    """
    Solve the exact Poisson equation backward from rho_start to rho_end
    using 4th-order Runge-Kutta.

    Returns arrays (rho_arr, psi_arr) from rho_end to rho_start (reversed to ascending order).
    """
    h = (rho_end - rho_start) / n_steps  # h is negative (backward integration)

    rho_arr = np.zeros(n_steps + 1)
    psi_arr = np.zeros(n_steps + 1)

    rho = rho_start
    # Initial conditions from approximate solution at rho_start
    w = rho * psi_approx(rho)        # w = rho * Psi
    dw = psi_approx(rho) + rho * (-psi_approx(rho) / rho - np.exp(-rho) / rho)
    # dw/drho for Psi_a = exp(-rho)/rho:
    # d[rho * exp(-rho)/rho]/drho = d[exp(-rho)]/drho = -exp(-rho)
    dw = -np.exp(-rho)

    rho_arr[0] = rho
    psi_arr[0] = w / rho if rho > 1e-15 else 0.0

    for i in range(n_steps):
        # RK4 for the second-order ODE
        # System: w' = dw, dw' = f(rho, w, dw)
        k1_w = h * dw
        k1_dw = h * rhs(rho, w, dw, N_D)

        k2_w = h * (dw + 0.5 * k1_dw)
        k2_dw = h * rhs(rho + 0.5 * h, w + 0.5 * k1_w, dw + 0.5 * k1_dw, N_D)

        k3_w = h * (dw + 0.5 * k2_dw)
        k3_dw = h * rhs(rho + 0.5 * h, w + 0.5 * k2_w, dw + 0.5 * k2_dw, N_D)

        k4_w = h * (dw + k3_dw)
        k4_dw = h * rhs(rho + h, w + k3_w, dw + k3_dw, N_D)

        w += (k1_w + 2 * k2_w + 2 * k3_w + k4_w) / 6.0
        dw += (k1_dw + 2 * k2_dw + 2 * k3_dw + k4_dw) / 6.0
        rho += h

        rho_arr[i + 1] = rho
        psi_arr[i + 1] = w / rho if abs(rho) > 1e-15 else 0.0

    # Reverse to ascending order
    return rho_arr[::-1], psi_arr[::-1]


def find_breakdown_distance(N_D, tolerance, rho_start=None, n_steps=500000):
    """
    Find the distance rho_d at which |Psi_a - Psi| / |Psi| * 100 = tolerance (%).

    Parameters:
        N_D: plasma parameter
        tolerance: percentage (1, 5, or 10)
        rho_start: starting rho for backward integration (auto if None)
        n_steps: number of integration steps

    Returns:
        rho_d: the distance at which the tolerance is reached
    """
    # Choose starting rho based on N_D (from Table 1)
    if rho_start is None:
        if N_D >= 1e7:
            rho_start = 0.001
        elif N_D >= 1e6:
            rho_start = 0.01
        elif N_D >= 1e5:
            rho_start = 0.1
        elif N_D >= 1e4:
            rho_start = 0.5
        elif N_D >= 1e3:
            rho_start = 2.0
        elif N_D >= 100:
            rho_start = 5.0
        elif N_D >= 40:
            rho_start = 5.0
        elif N_D >= 5:
            rho_start = 7.0
        elif N_D >= 1:
            rho_start = 11.0
        else:
            rho_start = 13.0

    # Adjust rho_end and n_steps for very large N_D
    if N_D >= 1e7:
        rho_end = 1e-10
        n_steps = max(n_steps, 1000000)
    elif N_D >= 1e6:
        rho_end = 1e-9
        n_steps = max(n_steps, 800000)
    elif N_D >= 1e5:
        rho_end = 1e-8
    else:
        rho_end = 1e-7
    rho_arr, psi_exact = solve_exact(N_D, rho_start, rho_end, n_steps)
    psi_a = psi_approx(rho_arr)

    # Compute relative error
    mask = psi_exact > 0
    rel_error = np.abs(psi_a[mask] - psi_exact[mask]) / psi_exact[mask] * 100.0
    rho_masked = rho_arr[mask]

    # Find where rel_error crosses tolerance (going from large to small rho)
    # The error increases as rho decreases, so we look for the crossing
    idx_above = np.where(rel_error >= tolerance)[0]
    if len(idx_above) == 0:
        return rho_end  # tolerance never reached

    # The first index where error >= tolerance (from large rho side)
    # Since arrays are in ascending order, we want the last index where error >= tolerance
    # Actually, error increases as rho decreases. In ascending rho order,
    # error decreases. So first index from left where error < tolerance is the crossing.
    idx_first_below = np.where(rel_error < tolerance)[0]
    if len(idx_first_below) == 0:
        return rho_masked[-1]

    # Linear interpolation at the crossing
    idx = idx_first_below[0]
    if idx == 0:
        return rho_masked[0]

    # Interpolate between idx-1 (above) and idx (below)
    rho1, e1 = rho_masked[idx - 1], rel_error[idx - 1]
    rho2, e2 = rho_masked[idx], rel_error[idx]
    if abs(e1 - e2) < 1e-15:
        return rho1
    rho_d = rho1 + (tolerance - e1) * (rho2 - rho1) / (e2 - e1)
    return rho_d


def compute_table1():
    """
    Reproduce Table 1 from the paper.
    Returns dict with N_D values and corresponding rho_d for 1%, 5%, 10% tolerances.
    """
    ND_values = [0.1, 0.3, 1, 2, 5, 10, 40, 1e2, 5e2, 1e3, 1e4, 1e5, 1e6, 1e7]
    results = []

    for N_D in ND_values:
        rho_1 = find_breakdown_distance(N_D, 1.0)
        rho_5 = find_breakdown_distance(N_D, 5.0)
        rho_10 = find_breakdown_distance(N_D, 10.0)
        results.append({
            'N_D': N_D,
            'rho_1pct': rho_1,
            'rho_5pct': rho_5,
            'rho_10pct': rho_10,
        })

    return results
