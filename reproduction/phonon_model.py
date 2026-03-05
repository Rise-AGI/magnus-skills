"""
Core module for FCC-Al phonon calculations using Born-von Karman force constants.

Reproduces the computational methodology from:
Togo & Tanaka, Scripta Materialia 108, 1-5 (2015)
"First principles phonon calculations in materials science"

All figures import from this module.
"""

import numpy as np
from scipy.constants import k as kB, hbar, N_A, eV
from scipy.optimize import minimize, minimize_scalar

AMU = 1.66053906660e-27  # kg

# Al physical constants
AL_MASS = 26.9815385 * AMU
AL_A0 = 4.05e-10  # m, equilibrium lattice constant
AL_V0 = AL_A0**3 / 4  # m^3 per atom (FCC has 4 atoms per conventional cell)
AL_B0 = 79e9  # Pa, bulk modulus
AL_B0P = 4.4  # B' (pressure derivative of bulk modulus)


def generate_fcc_shells(a, n_shells=4):
    """Generate neighbor shells for FCC lattice."""
    max_n = 4
    sites = []
    for i in range(-max_n, max_n + 1):
        for j in range(-max_n, max_n + 1):
            for k in range(-max_n, max_n + 1):
                if i == 0 and j == 0 and k == 0:
                    continue
                if (i + j + k) % 2 == 0:
                    r = np.array([i, j, k], dtype=float) * a / 2
                    d = np.linalg.norm(r)
                    sites.append((d, r))

    sites.sort(key=lambda x: x[0])

    shells = []
    current_d = -1
    tol = a * 1e-6
    for d, r in sites:
        if abs(d - current_d) > tol:
            if len(shells) >= n_shells:
                break
            shells.append([])
            current_d = d
        if len(shells) <= n_shells:
            shells[-1].append(r)

    return shells[:n_shells]


def dynamical_matrix(q, shells, fc_params, mass):
    """Compute 3x3 dynamical matrix at wave vector q.

    D_ab(q) = (1/M) Sum_R Phi_ab(R) [1 - exp(iq.R)]
    """
    D = np.zeros((3, 3), dtype=complex)

    for shell, (alpha, beta) in zip(shells, fc_params):
        for R in shell:
            d = np.linalg.norm(R)
            dhat = R / d
            Phi = alpha * np.outer(dhat, dhat) + beta * (np.eye(3) - np.outer(dhat, dhat))
            phase = np.exp(1j * np.dot(q, R))
            D += Phi * (1.0 - phase)

    D /= mass
    return D


def phonon_frequencies_Hz(q, shells, fc_params, mass):
    """Compute phonon frequencies in Hz at wave vector q."""
    D = dynamical_matrix(q, shells, fc_params, mass)
    D = 0.5 * (D + D.conj().T)
    eigenvalues = np.linalg.eigvalsh(D)
    freqs = np.sign(eigenvalues) * np.sqrt(np.abs(eigenvalues)) / (2 * np.pi)
    return np.sort(freqs)


def fit_al_force_constants(n_shells=4):
    """Fit Born-von Karman force constants to match Al phonon dispersion."""
    a = AL_A0
    M = AL_MASS
    shells = generate_fcc_shells(a, n_shells=n_shells)

    # Target frequencies in THz at high-symmetry points
    # From experimental neutron scattering data
    targets = [
        ([1.0, 0.0, 0.0], [5.82, 5.82, 9.69]),   # X
        ([0.5, 0.5, 0.5], [4.18, 4.18, 9.64]),   # L
        ([1.0, 0.5, 0.0], [5.41, 7.82, 9.31]),   # W
        ([0.75, 0.75, 0.0], [4.85, 6.50, 8.15]), # K
        ([0.5, 0.0, 0.0], [2.80, 2.80, 5.80]),   # Delta midpoint
        ([0.25, 0.25, 0.25], [2.10, 2.10, 5.20]), # Lambda midpoint
        ([0.5, 0.25, 0.0], [3.50, 4.20, 7.40]),   # Sigma midpoint
    ]

    def objective(params):
        fc = [(params[i * 2], params[i * 2 + 1]) for i in range(n_shells)]
        error = 0
        for q_frac, target_freqs in targets:
            q = np.array(q_frac) * 2 * np.pi / a
            freqs = phonon_frequencies_Hz(q, shells, fc, M) * 1e-12
            for f_calc, f_target in zip(sorted(np.abs(freqs)), sorted(target_freqs)):
                error += (f_calc - f_target) ** 2
        # Regularization to keep small-shell parameters small
        for i in range(n_shells):
            error += 0.01 * (params[i * 2] ** 2 + params[i * 2 + 1] ** 2) / (i + 1) ** 3
        return error

    x0 = [23.0, -2.5, -2.0, 0.8, 0.5, -0.2, 0.1, -0.05]
    result = minimize(objective, x0, method='Nelder-Mead',
                      options={'maxiter': 100000, 'xatol': 1e-8, 'fatol': 1e-10})

    fc = [(result.x[i * 2], result.x[i * 2 + 1]) for i in range(n_shells)]
    return fc, result.fun


def get_al_model(a=None):
    """Get fitted Al phonon model. Returns (shells, fc_params, mass, a)."""
    if a is None:
        a = AL_A0
    shells = generate_fcc_shells(a, n_shells=4)
    fc, _ = fit_al_force_constants()
    return shells, fc, AL_MASS, a


def scale_force_constants(fc_base, a, a0, power=6.0):
    """Scale force constants for different lattice constant.

    Phi(a) = Phi(a0) * (a0/a)^power
    power ~ 6*gamma where gamma is the Gruneisen parameter.
    For Al, gamma ~ 2.2, so power ~ 13. But empirically power ~ 6-8 works
    better with the BvK model.
    """
    scale = (a0 / a) ** power
    return [(alpha * scale, beta * scale) for alpha, beta in fc_base]


# ---------------------------------------------------------------------------
# Phonon band structure
# ---------------------------------------------------------------------------

# FCC high-symmetry points in fractional reciprocal coordinates (units of 2pi/a)
FCC_POINTS = {
    'G': np.array([0.0, 0.0, 0.0]),
    'X': np.array([1.0, 0.0, 0.0]),
    'W': np.array([1.0, 0.5, 0.0]),
    'K': np.array([0.75, 0.75, 0.0]),
    'U': np.array([1.0, 0.25, 0.25]),
    'L': np.array([0.5, 0.5, 0.5]),
}


def compute_bandstructure(path_labels, n_per_segment, shells, fc, mass, a):
    """Compute phonon band structure along a path.

    path_labels: list of point labels, e.g. ['L', 'G', 'X', 'W', 'K', 'G']
    Returns: (distances, frequencies_THz, tick_positions, tick_labels)
    """
    points = [FCC_POINTS[label.replace('Gamma', 'G').replace(u'\u0393', 'G')]
              for label in path_labels]

    all_dist = []
    all_freq = []
    tick_pos = [0.0]
    total_dist = 0.0

    for seg in range(len(points) - 1):
        q_start = points[seg] * 2 * np.pi / a
        q_end = points[seg + 1] * 2 * np.pi / a

        for j in range(n_per_segment):
            if seg == len(points) - 2:
                t = j / max(n_per_segment - 1, 1)
            else:
                t = j / n_per_segment

            q = q_start + t * (q_end - q_start)
            freqs = phonon_frequencies_Hz(q, shells, fc, mass) * 1e-12  # THz

            if len(all_dist) > 0:
                dq = np.linalg.norm(q - prev_q)
                total_dist += dq

            all_dist.append(total_dist)
            all_freq.append(freqs)
            prev_q = q

        tick_pos.append(total_dist)

    return np.array(all_dist), np.array(all_freq), tick_pos, path_labels


# ---------------------------------------------------------------------------
# Phonon DOS
# ---------------------------------------------------------------------------

def compute_dos(shells, fc, mass, a, n_grid=30, n_bins=200, sigma_THz=0.15):
    """Compute phonon DOS using Gaussian smearing on a uniform q-grid.

    Returns: (freq_THz, dos) where dos is in states/THz/atom, normalized so
    integral = 3 (3 modes per atom).
    """
    all_freqs = []

    for i in range(n_grid):
        for j in range(n_grid):
            for k in range(n_grid):
                q = np.array([i, j, k], dtype=float) / n_grid * 2 * np.pi / a
                freqs = phonon_frequencies_Hz(q, shells, fc, mass) * 1e-12
                all_freqs.extend(freqs)

    all_freqs = np.array(all_freqs)

    freq_max = np.max(np.abs(all_freqs)) * 1.1
    freq_grid = np.linspace(0, freq_max, n_bins)
    dos = np.zeros(n_bins)

    n_total = n_grid ** 3
    for f in all_freqs:
        if f > 0:
            dos += np.exp(-0.5 * ((freq_grid - f) / sigma_THz) ** 2) / (
                sigma_THz * np.sqrt(2 * np.pi))

    dos /= n_total  # Normalize per q-point (3 modes per point)
    # Check: integral should be 3
    df = freq_grid[1] - freq_grid[0]
    integral = np.sum(dos) * df
    if integral > 0:
        dos *= 3.0 / integral

    return freq_grid, dos


# ---------------------------------------------------------------------------
# Thermal properties from phonon DOS
# ---------------------------------------------------------------------------

def thermal_properties(freq_THz, dos, T_array):
    """Compute thermal properties from phonon DOS.

    Returns dict with arrays: Cv (J/K/mol), S (J/K/mol), F (kJ/mol).
    """
    omega = freq_THz * 1e12 * 2 * np.pi  # rad/s
    df_THz = freq_THz[1] - freq_THz[0]

    Cv = np.zeros_like(T_array, dtype=float)
    S = np.zeros_like(T_array, dtype=float)
    F = np.zeros_like(T_array, dtype=float)

    mask = omega > 1e6  # exclude zero frequency
    omega_m = omega[mask]
    dos_m = dos[mask]

    for i, T in enumerate(T_array):
        if T < 0.1:
            # Zero-point energy
            zpe = np.sum(0.5 * hbar * omega_m * dos_m * df_THz) * N_A
            F[i] = zpe / 1000.0
            continue

        x = hbar * omega_m / (kB * T)
        x = np.clip(x, 0, 500)

        exp_x = np.exp(x)

        # Cv (Eq. 9)
        cv_mode = kB * x ** 2 * exp_x / (exp_x - 1) ** 2
        Cv[i] = np.sum(cv_mode * dos_m * df_THz) * N_A

        # Free energy (Eq. 10)
        zpe = 0.5 * hbar * omega_m
        thermal = np.where(x < 500, kB * T * np.log(1 - np.exp(-x)), -hbar * omega_m)
        F[i] = np.sum((zpe + thermal) * dos_m * df_THz) * N_A / 1000.0

        # Entropy (Eq. 11)
        s_mode = np.where(x < 500,
                          kB * (x / (np.exp(x) - 1) - np.log(1 - np.exp(-x))),
                          0.0)
        S[i] = np.sum(s_mode * dos_m * df_THz) * N_A

    return {'Cv': Cv, 'S': S, 'F': F, 'T': T_array}


def phonon_free_energy_per_atom(freq_THz, dos, T):
    """Compute phonon Helmholtz free energy per atom in eV."""
    omega = freq_THz * 1e12 * 2 * np.pi
    df_THz = freq_THz[1] - freq_THz[0]
    mask = omega > 1e6
    omega_m = omega[mask]
    dos_m = dos[mask]

    if T < 0.1:
        return np.sum(0.5 * hbar * omega_m * dos_m * df_THz) / eV

    x = hbar * omega_m / (kB * T)
    x = np.clip(x, 0, 500)
    zpe = 0.5 * hbar * omega_m
    thermal = np.where(x < 500, kB * T * np.log(1 - np.exp(-x)), -hbar * omega_m)
    return np.sum((zpe + thermal) * dos_m * df_THz) / eV


# ---------------------------------------------------------------------------
# Equation of State (Vinet)
# ---------------------------------------------------------------------------

def vinet_energy(V, E0, V0, B0, B0p):
    """Vinet equation of state: E(V) in eV, V in m^3/atom."""
    x = (V / V0) ** (1.0 / 3.0)
    eta = 1.5 * (B0p - 1)
    E = E0 + (4 * V0 * B0 / (eta ** 2)) * (
        1 - (1 + eta * (1 - x)) * np.exp(-eta * (1 - x))
    )
    return E / eV  # Convert to eV


# ---------------------------------------------------------------------------
# QHA: Volume-dependent properties
# ---------------------------------------------------------------------------

def qha_properties(fc_base, a0, mass, volumes_per_atom, T_array,
                   E0=0, V0=None, B0=None, B0p=None, fc_power=7.0,
                   n_grid=20, n_bins=150):
    """Compute QHA properties.

    Parameters:
        fc_base: force constants at equilibrium
        a0: equilibrium lattice constant
        mass: atomic mass
        volumes_per_atom: array of volumes in m^3/atom
        T_array: temperatures in K
        E0, V0, B0, B0p: Vinet EOS parameters
        fc_power: power law for force constant scaling

    Returns dict with:
        freq_X, freq_L: phonon frequencies at X, L vs volume (THz)
        F_total: total free energy (eV/atom) shape (n_T, n_V)
        V_eq: equilibrium volume at each T
        beta: thermal expansion coefficient
        Cp: heat capacity at constant pressure
    """
    if V0 is None:
        V0 = AL_V0
    if B0 is None:
        B0 = AL_B0
    if B0p is None:
        B0p = AL_B0P

    n_V = len(volumes_per_atom)
    n_T = len(T_array)

    # Compute phonon DOS and free energy at each volume
    freq_X = np.zeros((n_V, 3))
    freq_L = np.zeros((n_V, 3))
    F_ph = np.zeros((n_V, n_T))

    for iv, V in enumerate(volumes_per_atom):
        a = (4 * V) ** (1.0 / 3.0)
        shells = generate_fcc_shells(a, n_shells=4)
        fc = scale_force_constants(fc_base, a, a0, power=fc_power)

        # Frequencies at X and L
        q_X = np.array([1, 0, 0], dtype=float) * 2 * np.pi / a
        q_L = np.array([0.5, 0.5, 0.5]) * 2 * np.pi / a
        freq_X[iv] = phonon_frequencies_Hz(q_X, shells, fc, mass) * 1e-12
        freq_L[iv] = phonon_frequencies_Hz(q_L, shells, fc, mass) * 1e-12

        # DOS
        freq_grid, dos = compute_dos(shells, fc, mass, a, n_grid=n_grid,
                                     n_bins=n_bins, sigma_THz=0.15)

        for it, T in enumerate(T_array):
            F_ph[iv, it] = phonon_free_energy_per_atom(freq_grid, dos, T)

    # Electronic energy (Vinet EOS)
    E_el = np.array([vinet_energy(V, E0, V0, B0, B0p) for V in volumes_per_atom])

    # Total free energy F_total(V, T) = E_el(V) + F_ph(V, T)
    F_total = E_el[:, np.newaxis] + F_ph  # shape (n_V, n_T)

    # Find equilibrium volume at each temperature
    V_eq = np.zeros(n_T)
    for it in range(n_T):
        # Fit polynomial and find minimum
        coeffs = np.polyfit(volumes_per_atom, F_total[:, it], 4)
        p = np.poly1d(coeffs)
        dp = p.deriv()
        # Find root of dp in volume range
        try:
            res = minimize_scalar(p, bounds=(volumes_per_atom[0], volumes_per_atom[-1]),
                                  method='bounded')
            V_eq[it] = res.x
        except Exception:
            V_eq[it] = V0

    # Thermal expansion coefficient
    beta_T = np.zeros(n_T)
    dT = T_array[1] - T_array[0] if len(T_array) > 1 else 1.0
    for it in range(1, n_T - 1):
        beta_T[it] = (V_eq[it + 1] - V_eq[it - 1]) / (2 * dT * V_eq[it])
    if n_T > 2:
        beta_T[0] = (V_eq[1] - V_eq[0]) / (dT * V_eq[0])
        beta_T[-1] = (V_eq[-1] - V_eq[-2]) / (dT * V_eq[-1])

    return {
        'volumes': volumes_per_atom,
        'freq_X': freq_X,
        'freq_L': freq_L,
        'E_el': E_el,
        'F_ph': F_ph,
        'F_total': F_total,
        'V_eq': V_eq,
        'beta': beta_T,
        'T': T_array,
    }
