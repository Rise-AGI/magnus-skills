"""
Core implementation of the Split-Step Fourier (SSF) method for the
Nonlinear Schrodinger Equation (NLSE):

    i A_t + sigma * d/2 * A_xx + |A|^2 * A = 0

as described in Zoldi et al. (1997), arXiv:physics/9711012.

Also implements the parallel FFT decomposition via 2D matrix factorization
(Eq. 2-4 of the paper).
"""

import numpy as np
import time


def serial_ssf_step(A, dt, k_freq, sigma, d):
    """One step of serial Split-Step Fourier method.

    Steps (from paper Section 2):
    1. Nonlinear step: A1 = exp(dt * N) * A  where N = i|A|^2
    2. Forward FFT: A2 = FFT(A1)
    3. Linear step:  A3 = exp(-i * sigma * d/2 * k^2 * dt) * A2
    4. Backward FFT: A = IFFT(A3)
    """
    # Step 1: Nonlinear step
    A1 = A * np.exp(1j * np.abs(A)**2 * dt)
    # Step 2-3-4: Forward FFT, linear step, backward FFT
    A2 = np.fft.fft(A1)
    A3 = A2 * np.exp(-1j * sigma * d / 2.0 * k_freq**2 * dt)
    return np.fft.ifft(A3)


def parallel_fft_2d(A, M0, M1):
    """Parallel 1D FFT via 2D matrix decomposition (Eq. 2-4 of paper).

    The 1D array A of size N = M0*M1 is written as a 2D matrix a[j,k]
    of size M0 x M1, then:
      Step 1: Independent M1-size FFTs on rows
      Step 2: Multiply by twiddle factors E[j,k] = exp(-2*pi*i*j*k/N)
      Step 3: Independent M0-size FFTs on columns

    Result is stored in rows: F[M1*k1+k0] in row k1, column k0.
    """
    N = M0 * M1
    # Reshape into M0 x M1 matrix (columns: A[0..M0-1], A[M0..2M0-1], ...)
    # a[n0, n1] = A[M0*n1 + n0]
    a = A.reshape(M1, M0).T.copy()  # a[n0, n1]

    # Step 1: M0 independent M1-size FFTs on rows (fixed n0)
    for n0 in range(M0):
        a[n0, :] = np.fft.fft(a[n0, :])

    # Step 2: Multiply by twiddle factors E[n0, k0] = exp(-2*pi*i*n0*k0/N)
    n0_arr = np.arange(M0).reshape(-1, 1)
    k0_arr = np.arange(M1).reshape(1, -1)
    E = np.exp(-2j * np.pi * n0_arr * k0_arr / N)
    a = a * E

    # Step 3: M1 independent M0-size FFTs on columns (fixed k0)
    for k0 in range(M1):
        a[:, k0] = np.fft.fft(a[:, k0])

    # Result: F[M1*k1 + k0] is in a[k1, k0]
    F = np.zeros(N, dtype=complex)
    for k1 in range(M0):
        for k0 in range(M1):
            F[M1 * k1 + k0] = a[k1, k0]
    return F


def parallel_bft_2d(F, M0, M1):
    """Inverse of parallel_fft_2d: backward FFT via 2D decomposition.

    Reverses Steps 1-3 of parallel_fft_2d with conjugate twiddle factors.
    """
    N = M0 * M1
    # Unpack: F[M1*k1 + k0] -> a[k1, k0]
    a = np.zeros((M0, M1), dtype=complex)
    for k1 in range(M0):
        for k0 in range(M1):
            a[k1, k0] = F[M1 * k1 + k0]

    # Step 3 inverse: M1 independent M0-size IFFTs on columns
    for k0 in range(M1):
        a[:, k0] = np.fft.ifft(a[:, k0])

    # Step 2 inverse: Multiply by conjugate twiddle factors
    n0_arr = np.arange(M0).reshape(-1, 1)
    k0_arr = np.arange(M1).reshape(1, -1)
    E_conj = np.exp(2j * np.pi * n0_arr * k0_arr / N)
    a = a * E_conj

    # Step 1 inverse: M0 independent M1-size IFFTs on rows
    for n0 in range(M0):
        a[n0, :] = np.fft.ifft(a[n0, :])

    # Reconstruct: A[M0*n1 + n0] = a[n0, n1]
    A = a.T.reshape(N)
    return A


def parallel_ssf_step(A, dt, k_freq_transposed, sigma, d, M0, M1):
    """One step of parallel SSF (Section 3.1 of paper).

    Steps 1-8 from the paper:
    1. Nonlinear step
    2. Row-FFT
    3. Multiply by E
    4. Column-FFT
    5. Linear step (transposed linear operator)
    6. Column-BFT
    7. Multiply by E*
    8. Row-BFT
    """
    N = M0 * M1
    # Step 1: Nonlinear step
    A1 = A * np.exp(1j * np.abs(A)**2 * dt)

    # Steps 2-4: Forward FFT via 2D decomposition
    A_freq = parallel_fft_2d(A1, M0, M1)

    # Step 5: Linear step with transposed operator
    A_freq = A_freq * np.exp(-1j * sigma * d / 2.0 * k_freq_transposed**2 * dt)

    # Steps 6-8: Backward FFT via 2D decomposition
    return parallel_bft_2d(A_freq, M0, M1)


def soliton_initial(x, amplitude=1.0, velocity=0.0, x0=0.0):
    """Bright soliton solution of NLSE: A(x,0) = a * sech(a*x) * exp(i*v*x).

    For sigma=1, d=1, G=0, the exact solution is:
    A(x,t) = a * sech(a*(x - v*t)) * exp(i*(v*x + (a^2 - v^2/2)*t))
    """
    a = amplitude
    return a / np.cosh(a * (x - x0)) * np.exp(1j * velocity * x)


def run_ssf_simulation(N, L, T_total, n_steps, sigma=1.0, d=1.0,
                       initial_func=None, use_parallel=False):
    """Run SSF simulation of NLSE.

    Args:
        N: Number of spatial grid points (must be power of 2 for parallel)
        L: Spatial domain half-width [-L, L]
        T_total: Total time to simulate
        n_steps: Number of time steps
        sigma: Dispersion sign (+1 anomalous, -1 normal)
        d: Dispersion coefficient
        initial_func: callable(x) -> A(x, 0)
        use_parallel: If True, use 2D decomposed FFT

    Returns:
        x, t_array, A_history (snapshots)
    """
    dx = 2 * L / N
    x = np.linspace(-L, L - dx, N)
    dt = T_total / n_steps

    # Frequency array for linear operator
    k = np.fft.fftfreq(N, d=dx) * 2 * np.pi

    if initial_func is None:
        initial_func = soliton_initial
    A = initial_func(x).astype(complex)

    # For parallel method
    M0 = M1 = None
    k_transposed = None
    if use_parallel:
        M0 = int(np.sqrt(N))
        M1 = N // M0
        assert M0 * M1 == N, f"N={N} must be a perfect square for M0=M1"
        # Build transposed frequency array
        k_transposed = np.zeros(N, dtype=float)
        for k1 in range(M0):
            for k0 in range(M1):
                k_transposed[M1 * k1 + k0] = k[M1 * k1 + k0]

    n_snapshots = min(100, n_steps)
    snapshot_interval = max(1, n_steps // n_snapshots)
    A_history = [np.abs(A).copy()]
    t_array = [0.0]

    for step in range(1, n_steps + 1):
        if use_parallel:
            A = parallel_ssf_step(A, dt, k_transposed, sigma, d, M0, M1)
        else:
            A = serial_ssf_step(A, dt, k, sigma, d)

        if step % snapshot_interval == 0 or step == n_steps:
            A_history.append(np.abs(A).copy())
            t_array.append(step * dt)

    return x, np.array(t_array), np.array(A_history)


def time_fft_methods(N_values, n_repeats=10):
    """Time serial vs parallel FFT for different array sizes.

    Returns dict with N -> (serial_time, parallel_time).
    """
    results = {}
    for N in N_values:
        M = int(np.sqrt(N))
        if M * M != N:
            continue

        A = np.random.randn(N) + 1j * np.random.randn(N)

        # Serial FFT timing
        t0 = time.perf_counter()
        for _ in range(n_repeats):
            np.fft.fft(A)
        serial_time = (time.perf_counter() - t0) / n_repeats

        # Parallel (2D decomposed) FFT timing
        t0 = time.perf_counter()
        for _ in range(n_repeats):
            parallel_fft_2d(A, M, M)
        parallel_time = (time.perf_counter() - t0) / n_repeats

        results[N] = (serial_time, parallel_time)

    return results


def theoretical_speedup(K, P, xi=1.0, f=0.001):
    """Theoretical speedup from Eq. (6) of the paper.

    SU = 2P / (1 + xi/K + f * 2^K / (P*K))

    where N = 2^(2K), P = number of processors.
    """
    return 2 * P / (1.0 + xi / K + f * 2**K / (P * K))
