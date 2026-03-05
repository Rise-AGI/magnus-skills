"""
Core module for beam-plasma instability growth rates.
Bret (2009), ApJ 699:990, arXiv:0903.2658.

Solves the full cold-fluid dielectric tensor for a relativistic electron beam
passing through a magnetized cold plasma with return current and mobile ions.
"""

import numpy as np


def plasma_params(alpha, gamma_b, R=1.0/1836.0):
    """Derive beta, gamma_p from alpha and gamma_b.

    Parameters
    ----------
    alpha : beam-to-plasma density ratio nb/np
    gamma_b : beam Lorentz factor
    R : electron-to-proton mass ratio

    Returns
    -------
    beta, gamma_p : beam velocity (v/c), return-current Lorentz factor
    """
    beta = np.sqrt(1.0 - 1.0 / gamma_b**2)
    beta_p = alpha * beta
    if beta_p >= 1.0:
        raise ValueError("Return current superluminal")
    gamma_p = 1.0 / np.sqrt(1.0 - beta_p**2)
    return beta, gamma_p


def newton_complex(f, x0, niter=300, tol=1e-13, dx=1e-8):
    """Newton's method for f(x) = 0 in the complex plane."""
    x = x0
    for _ in range(niter):
        fx = f(x)
        if abs(fx) < tol:
            return x
        if not np.isfinite(abs(fx)):
            return x0  # diverged to pole
        dfx = (f(x + dx) - f(x - dx)) / (2 * dx)
        if abs(dfx) < 1e-30 or not np.isfinite(abs(dfx)):
            return x
        step = fx / dfx
        # Limit step size to avoid jumping past poles
        if abs(step) > 1.0:
            step = step / abs(step) * 1.0
        x = x - step
    return x


# --- Parallel propagation (Zx = 0) ---

def Tzz_par(x, Zz, alpha, gamma_b, R):
    """Electrostatic dispersion for k parallel to beam (Eq. 9)."""
    beta, gamma_p = plasma_params(alpha, gamma_b, R)
    return (1.0
            - R * (1 + alpha) / x**2
            - alpha / ((x - Zz)**2 * gamma_b**3)
            - 1.0 / ((x + alpha * Zz)**2 * gamma_p**3))


def Txx_pm_iTxy_par(x, Zz, alpha, gamma_b, OmegaB, R, sign):
    """EM dispersion for k parallel: Txx + sign*i*Txy (Eq. 8)."""
    beta, gamma_p = plasma_params(alpha, gamma_b, R)
    xb = x - Zz
    xp = x + alpha * Zz

    # T2_xx (Eqs A3-A4)
    t2xx = (alpha * (-xb) / (xb * gamma_b + OmegaB)
            + alpha * xb / (-xb * gamma_b + OmegaB)
            - xp / (xp * gamma_p + OmegaB)
            + xp / (-xp * gamma_p + OmegaB)
            + R * x * (1 + alpha) / (-x + R * OmegaB)
            - R * x * (1 + alpha) / (x + R * OmegaB))

    # T2_xy (Eq A9) - purely imaginary
    t2xy = (2j * R**2 * x * OmegaB * (1 + alpha) / (-x**2 + R**2 * OmegaB**2)
            + 2j * OmegaB * (xp / (xp**2 * gamma_p**2 - OmegaB**2)
                             + alpha * (-xb) / (-xb**2 * gamma_b**2 + OmegaB**2)))

    # n^2 = Zz^2 / (x^2 * beta^2)
    n2 = Zz**2 / (x**2 * beta**2) if abs(x) > 1e-30 else 0.0

    Txx = 1.0 + t2xx - n2
    return Txx + sign * 1j * t2xy


def scan_parallel_modes(Zz_array, alpha, gamma_b, OmegaB, R,
                        n_re=60, n_im=30, x_re_max=None, x_im_max=None):
    """Scan growth rates of parallel modes (Zx=0) vs Zz.

    For each Zz, finds roots of:
      - Tzz = 0 (electrostatic: Two-Stream / Buneman)
      - Txx - iTxy = 0 (EM minus, sign=-1)
      - Txx + iTxy = 0 (EM plus, sign=+1)

    Uses continuation + grid scan + Newton refinement.
    """
    beta, gamma_p = plasma_params(alpha, gamma_b, R)
    n = len(Zz_array)
    delta_es = np.zeros(n)
    delta_em_m = np.zeros(n)
    delta_em_p = np.zeros(n)

    # Adaptive scan ranges
    if x_re_max is None:
        x_re_max = max(3.0, Zz_array.max() * 1.5)
    if x_im_max is None:
        x_im_max = 0.3

    prev_roots_es = []
    prev_roots_em_m = []
    prev_roots_em_p = []

    for i, Zz in enumerate(Zz_array):
        # --- Electrostatic ---
        def fes(x):
            return Tzz_par(x, Zz, alpha, gamma_b, R)

        guesses = list(prev_roots_es)
        # Beam resonance region
        for dr in [-0.15, -0.1, -0.05, 0.0, 0.05]:
            for di in [0.005, 0.01, 0.02, 0.05]:
                guesses.append(complex(Zz + dr, di))
        # Return current resonance region (Buneman)
        for dr in [-0.1, 0.0, 0.1]:
            for di in [0.01, 0.03, 0.056]:
                guesses.append(complex(-alpha * Zz + dr, di))
        # Proton resonance
        for di in [0.001, 0.005]:
            guesses.append(complex(0.0, di))

        best_d, best_x = _find_best_root(fes, guesses)
        delta_es[i] = best_d
        prev_roots_es = [best_x] if best_x is not None else []

        # --- EM minus (Txx - iTxy = 0) ---
        def fem_m(x):
            return Txx_pm_iTxy_par(x, Zz, alpha, gamma_b, OmegaB, R, -1)

        guesses = list(prev_roots_em_m)
        # Bell-like modes near Zz ~ OmegaB/gamma_b
        em_bound = R * OmegaB if OmegaB > 0 else 0.01
        for dr in [-0.1, 0.0, 0.1]:
            for di_frac in [0.1, 0.3, 0.5, 0.7, 0.95]:
                di = em_bound * di_frac
                if di > 1e-8:
                    guesses.append(complex(OmegaB / gamma_b + dr, di))
        for di in [1e-5, 1e-4, 5e-4]:
            guesses.append(complex(Zz * 0.5, di))

        # EM growth rates bounded by R*OmegaB (Eq. 13)
        em_max = R * OmegaB + 0.01 if OmegaB > 0 else 0.5
        best_d, best_x = _find_best_root(fem_m, guesses, max_growth=em_max)
        delta_em_m[i] = best_d
        prev_roots_em_m = [best_x] if best_x is not None else []

        # --- EM plus (Txx + iTxy = 0) ---
        def fem_p(x):
            return Txx_pm_iTxy_par(x, Zz, alpha, gamma_b, OmegaB, R, +1)

        guesses = list(prev_roots_em_p)
        if OmegaB > 0:
            for dr in [-0.1, 0.0, 0.1]:
                for di_frac in [0.1, 0.3, 0.5, 0.7]:
                    di = em_bound * di_frac
                    if di > 1e-8:
                        guesses.append(complex(0.5 * alpha / OmegaB + dr, di))
        for di in [1e-5, 1e-4, 5e-4]:
            guesses.append(complex(Zz * 0.3, di))

        best_d, best_x = _find_best_root(fem_p, guesses, max_growth=em_max)
        delta_em_p[i] = best_d
        prev_roots_em_p = [best_x] if best_x is not None else []

    return delta_es, delta_em_m, delta_em_p


def _find_best_root(f, guesses, tol=1e-8, max_growth=2.0):
    """Find the root with largest positive imaginary part.

    max_growth: physical upper bound on delta (rejects spurious poles).
    """
    best_d = 0.0
    best_x = None
    for g in guesses:
        try:
            root = newton_complex(f, g)
            fval = f(root)
            if (np.isfinite(abs(fval)) and abs(fval) < tol
                    and root.imag > best_d and root.imag < max_growth):
                best_d = root.imag
                best_x = root
        except (OverflowError, ZeroDivisionError, ValueError):
            pass
    return best_d, best_x


# --- Full 2D spectrum (general Zx, Zz) ---

def build_tensor(x, Zx, Zz, alpha, gamma_b, OmegaB, R):
    """Build the full 3x3 dispersion tensor."""
    beta, gamma_p = plasma_params(alpha, gamma_b, R)
    gb = gamma_b
    gp = gamma_p
    xb = x - Zz
    xp = x + alpha * Zz

    # Susceptibilities (Appendix A)
    t2xx = (alpha * (-xb) / (xb * gb + OmegaB)
            + alpha * xb / (-xb * gb + OmegaB)
            - xp / (xp * gp + OmegaB)
            + xp / (-xp * gp + OmegaB)
            + R * x * (1 + alpha) / (-x + R * OmegaB)
            - R * x * (1 + alpha) / (x + R * OmegaB))

    t2xy = (2j * R**2 * x * OmegaB * (1 + alpha) / (-x**2 + R**2 * OmegaB**2)
            + 2j * OmegaB * (xp / (xp**2 * gp**2 - OmegaB**2)
                             + alpha * (-xb) / (-xb**2 * gb**2 + OmegaB**2)))

    t2xz = -2 * alpha * Zx * (
        (-xb) * gb / (-xb**2 * gb**2 + OmegaB**2)
        + xp * gp / (-xp**2 * gp**2 + OmegaB**2))

    t2yz = 2j * alpha * Zx * OmegaB * (
        1.0 / (xp**2 * gp**2 - OmegaB**2)
        + 1.0 / (-xb**2 * gb**2 + OmegaB**2))

    # T2_zz (Eqs A7-A8)
    half_tzz = -R * (1 + alpha) + x**2 * (
        -alpha / (xb**2 * gb**3) - 1.0 / (xp**2 * gp**3))
    if abs(Zx) > 1e-15:
        num = gb * gp * (-alpha * xb**2 * gb - xp**2 * gp) + (gb + alpha * gp) * OmegaB**2
        den = (xb**2 * gb**2 - OmegaB**2) * (xp**2 * gp**2 - OmegaB**2)
        half_tzz += alpha * Zx**2 * num / den
    t2zz = 2 * half_tzz

    # Wave terms
    Z2 = Zx**2 + Zz**2
    n2 = Z2 / (x**2 * beta**2) if abs(x) > 1e-30 else 0.0

    T = np.zeros((3, 3), dtype=complex)
    if Z2 > 0:
        T[0, 0] = 1.0 + t2xx - n2 * (1 - Zx**2 / Z2)
        T[1, 1] = 1.0 + t2xx - n2
        T[2, 2] = 1.0 + t2zz - n2 * (1 - Zz**2 / Z2)
        T[0, 2] = t2xz + n2 * Zx * Zz / Z2
        T[2, 0] = T[0, 2]
    else:
        T[0, 0] = 1.0 + t2xx
        T[1, 1] = 1.0 + t2xx
        T[2, 2] = 1.0 + t2zz

    T[0, 1] = t2xy
    T[1, 0] = -t2xy
    T[1, 2] = t2yz
    T[2, 1] = -t2yz

    return T


def det_dispersion(x, Zx, Zz, alpha, gamma_b, OmegaB, R):
    """det(T) for the full dispersion equation."""
    T = build_tensor(x, Zx, Zz, alpha, gamma_b, OmegaB, R)
    return np.linalg.det(T)


def growth_rate_at_point(Zx, Zz, alpha, gamma_b, OmegaB, R,
                         x_re_pts=30, x_im_pts=15, x_im_max=0.3):
    """Find the maximum growth rate at a given (Zx, Zz) point."""
    beta, gamma_p = plasma_params(alpha, gamma_b, R)

    def f(x):
        return det_dispersion(x, Zx, Zz, alpha, gamma_b, OmegaB, R)

    # Build initial guess grid
    guesses = []
    x_re_range = np.linspace(-0.5, max(2.0, Zz * 1.5), x_re_pts)
    x_im_range = np.linspace(0.003, x_im_max, x_im_pts)

    for xr in x_re_range:
        for xi in x_im_range:
            guesses.append(complex(xr, xi))

    # Add physics-motivated guesses
    # Near beam resonance
    for di in [0.005, 0.01, 0.02, 0.05]:
        guesses.append(complex(Zz * 0.95, di))
    # Near return current resonance
    for di in [0.01, 0.05]:
        guesses.append(complex(-alpha * Zz, di))

    best_d = 0.0
    for g in guesses:
        try:
            root = newton_complex(f, g, niter=100)
            if abs(f(root)) < 1e-6 and root.imag > best_d:
                best_d = root.imag
        except (OverflowError, ZeroDivisionError, ValueError):
            pass

    return best_d


def compute_2d_spectrum(Zx_arr, Zz_arr, alpha, gamma_b, OmegaB, R, verbose=False):
    """Compute 2D growth rate map."""
    nZx = len(Zx_arr)
    nZz = len(Zz_arr)
    delta = np.zeros((nZx, nZz))

    for i, Zx in enumerate(Zx_arr):
        if verbose and i % 5 == 0:
            print(f"  Zx row {i}/{nZx}")
        for j, Zz in enumerate(Zz_arr):
            delta[i, j] = growth_rate_at_point(Zx, Zz, alpha, gamma_b, OmegaB, R)

    return delta
