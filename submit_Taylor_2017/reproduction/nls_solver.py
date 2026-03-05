"""
Core numerical solvers for the Nonlinear Schrodinger Equation with saturation nonlinearity.

PDE: i * dpsi/dt + (1/2) * d^2psi/dx^2 + |psi|^2 * psi / (1 + S * |psi|^2) = 0

Methods implemented:
  1. Split-Step Fourier (SSF)
  2. Central Finite Difference (CFD)

Reference: Taylor (2017), "A comparison between the Split Step Fourier and
Finite-Difference method in analysing the soliton collision of a type of
Nonlinear Schrodinger equation found in the context of optical pulses"
"""

import numpy as np


def soliton(x, S, v=0.0, t=0.0):
    """
    Exact soliton solution of the NLS with saturation nonlinearity.

    psi(x,t) = [2*sqrt(2)*exp(sqrt(2)*x)] / [1 + (3/2 - 2*S)*exp(2*sqrt(2)*x)] * exp(i*t + i*v*x)

    Parameters
    ----------
    x : array_like
        Spatial coordinates.
    S : float
        Saturation parameter.
    v : float
        Velocity parameter (phase gradient).
    t : float
        Time.

    Returns
    -------
    psi : ndarray
        Complex wavefunction values.
    """
    a = np.sqrt(2.0)
    B = 1.5 - 2.0 * S
    numerator = 2.0 * a * np.exp(a * x)
    denominator = 1.0 + B * np.exp(2.0 * a * x)
    # Avoid division by zero for large |x| or bad S
    with np.errstate(over='ignore', invalid='ignore'):
        f = np.where(np.abs(denominator) > 1e-300, numerator / denominator, 0.0)
    psi = f * np.exp(1j * t + 1j * v * x)
    return psi


def init_two_soliton(x1_off, v1, x2_off, v2, L, N, S):
    """
    Initialize two-soliton configuration as superposition.

    Parameters
    ----------
    x1_off, x2_off : float
        Center positions of soliton 1 and 2.
    v1, v2 : float
        Velocities of soliton 1 and 2.
    L : float
        Domain length [0, L].
    N : int
        Number of spatial grid points.
    S : float
        Saturation parameter.

    Returns
    -------
    psi : ndarray
        Initial complex wavefunction.
    x : ndarray
        Spatial grid.
    """
    h = L / N
    x = np.arange(N) * h
    psi1 = soliton(x - x1_off, S, v=v1, t=0.0)
    psi2 = soliton(x - x2_off, S, v=v2, t=0.0)
    return psi1 + psi2, x


def init_one_soliton(x_off, v, L, N, S):
    """Initialize a single soliton."""
    h = L / N
    x = np.arange(N) * h
    psi = soliton(x - x_off, S, v=v, t=0.0)
    return psi, x


def compute_norm(psi, dx):
    """Compute the integral of |psi|^2 using trapezoidal rule."""
    return np.trapz(np.abs(psi)**2, dx=dx)


# ============================================================
# Split-Step Fourier Method
# ============================================================

def splitstep_advance(psi, S, tau, L, N):
    """
    Advance psi by one time step using the Split-Step Fourier method.

    Steps:
      1. Nonlinear half-step: psi *= exp(i * |psi|^2 / (1 + S*|psi|^2) * tau)
      2. Linear step in Fourier space: multiply by exp(-i * k^2/2 * tau)
      3. (Strang splitting uses half-steps; here we use Lie splitting as in the paper)

    Parameters
    ----------
    psi : ndarray
        Current wavefunction.
    S : float
        Saturation parameter.
    tau : float
        Time step.
    L : float
        Domain length.
    N : int
        Number of grid points.

    Returns
    -------
    psi_new : ndarray
        Updated wavefunction.
    """
    # Nonlinear part
    abs_psi_sq = np.abs(psi)**2
    nonlinear_phase = abs_psi_sq / (1.0 + S * abs_psi_sq) * tau
    psi = psi * np.exp(1j * nonlinear_phase)

    # Linear part in Fourier space
    psi_hat = np.fft.fftshift(np.fft.fft(psi))
    # Wavenumbers (shifted)
    k = np.arange(N) - N // 2
    k_phys = 2.0 * np.pi * k / L
    linear_phase = -0.5 * k_phys**2 * tau
    psi_hat *= np.exp(1j * linear_phase)
    psi = np.fft.ifft(np.fft.fftshift(psi_hat))

    return psi


def run_splitstep(psi0, S, T, tau, L, N, store_every=1):
    """
    Run Split-Step simulation.

    Parameters
    ----------
    psi0 : ndarray
        Initial wavefunction.
    S : float
        Saturation parameter.
    T : float
        Total simulation time.
    tau : float
        Time step.
    L : float
        Domain length.
    N : int
        Number of grid points.
    store_every : int
        Store solution every this many steps.

    Returns
    -------
    psi_evolution : ndarray
        Shape (n_stored, N), |psi| at stored time steps.
    times : ndarray
        Time values at stored steps.
    """
    n_steps = int(T / tau)
    stored = []
    times = []
    psi = psi0.copy()

    for step in range(n_steps):
        if step % store_every == 0:
            stored.append(np.abs(psi).copy())
            times.append(step * tau)
        psi = splitstep_advance(psi, S, tau, L, N)

    # Store final state
    stored.append(np.abs(psi).copy())
    times.append(n_steps * tau)

    return np.array(stored), np.array(times)


# ============================================================
# Finite Difference Method (Central Difference)
# ============================================================

def fd_forward(psi, S, tau, h):
    """
    Forward Difference step (used for first time step only).

    psi_{j,k+1} = tau * i * (0.5 * psi_xx + A_{j,k} * psi_{j,k}) + psi_{j,k}
    """
    N = len(psi)
    psi_new = np.zeros(N, dtype=np.complex128)

    for j in range(N):
        # Periodic boundary conditions
        jm1 = (j - 1) % N
        jp1 = (j + 1) % N
        psi_xx = (psi[jm1] - 2.0 * psi[j] + psi[jp1]) / (h * h)
        abs_sq = np.abs(psi[j])**2
        A = abs_sq / (1.0 + S * abs_sq)
        psi_new[j] = 1j * tau * (0.5 * psi_xx + A * psi[j]) + psi[j]

    return psi_new


def fd_central(psi_curr, psi_prev, S, tau, h):
    """
    Central Difference step.

    psi_{j,k+1} = 2*tau*i*(0.5*psi_xx + A_{j,k}*psi_{j,k}) + psi_{j,k-1}
    """
    N = len(psi_curr)
    psi_new = np.zeros(N, dtype=np.complex128)

    for j in range(N):
        jm1 = (j - 1) % N
        jp1 = (j + 1) % N
        psi_xx = (psi_curr[jm1] - 2.0 * psi_curr[j] + psi_curr[jp1]) / (h * h)
        abs_sq = np.abs(psi_curr[j])**2
        A = abs_sq / (1.0 + S * abs_sq)
        psi_new[j] = 2.0 * 1j * tau * (0.5 * psi_xx + A * psi_curr[j]) + psi_prev[j]

    return psi_new


def run_finitediff(psi0, S, T, tau, L, N, store_every=1):
    """
    Run Finite Difference (Central Difference) simulation.

    Parameters
    ----------
    psi0 : ndarray
        Initial wavefunction.
    S : float
        Saturation parameter.
    T : float
        Total simulation time.
    tau : float
        Time step.
    L : float
        Domain length.
    N : int
        Number of grid points.
    store_every : int
        Store solution every this many steps.

    Returns
    -------
    psi_evolution : ndarray
        Shape (n_stored, N), |psi| at stored time steps.
    times : ndarray
        Time values at stored steps.
    """
    h = L / N
    n_steps = int(T / tau)
    stored = []
    times = []

    psi_prev = psi0.copy()
    psi_curr = fd_forward(psi_prev, S, tau, h)

    stored.append(np.abs(psi_prev).copy())
    times.append(0.0)

    for step in range(1, n_steps):
        if step % store_every == 0:
            stored.append(np.abs(psi_curr).copy())
            times.append(step * tau)
        psi_next = fd_central(psi_curr, psi_prev, S, tau, h)
        psi_prev = psi_curr
        psi_curr = psi_next

    stored.append(np.abs(psi_curr).copy())
    times.append(n_steps * tau)

    return np.array(stored), np.array(times)
