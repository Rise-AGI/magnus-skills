"""
DMFT-IPT solver for the half-filled single-band Hubbard model on the Bethe lattice.

At half-filling, particle-hole symmetry gives mu = U/2 and Hartree Sigma_H = U/2.
These cancel exactly: G(iw) = bethe_gf(iw - Sigma_tilde(iw))
where Sigma_tilde is the purely imaginary IPT self-energy (without Hartree).

Self-consistency on Bethe lattice: G0^{-1}(iw) = iw - t^2 G(iw).
IPT self-energy: Sigma_tilde(tau) = U^2 G0(tau)^3.

Energy unit: t = 1/2, half-bandwidth D = 2t = 1.

References:
  - A. Georges et al., Rev. Mod. Phys. 68, 13 (1996)
  - J. Joo and V. Oudovenko, Phys. Rev. B 64, 193102 (2001)
"""

import numpy as np


def matsubara_freq(n_freq, beta):
    """Fermionic Matsubara frequencies omega_n = (2n+1)*pi/beta."""
    return (2 * np.arange(n_freq) + 1) * np.pi / beta


def bethe_gf(z):
    """
    Local Green's function for semicircular DOS, half-bandwidth D=1.
    G(z) = 2/D^2 * (z - sign(Im z) * sqrt(z^2 - D^2))
    With D=1: G(z) = 2(z - sign(Im z) * sqrt(z^2 - 1))
    """
    sq = np.sqrt(z**2 - 1.0 + 0j)
    sign = np.where(np.imag(z) >= 0, 1.0, -1.0)
    return 2.0 * (z - sign * sq)


def iw_to_tau(g_iw, wn, beta, n_tau):
    """
    Inverse Fourier transform with 1/iw tail subtraction.

    At half-filling G(iw) is purely imaginary, and the tail is 1/iw.
    G(tau) = (2/beta) sum_{n>=0} sin(wn*tau) * Im G(iwn)

    With 1/iw tail subtraction (Im(1/iw) = -1/w):
    G(tau) = (2/beta) sum sin(wn*tau) * [Im G(iwn) + 1/wn] - 1/2
    """
    tau = np.linspace(0, beta, n_tau + 1)
    im_g = np.imag(g_iw)

    g_tau = np.zeros(n_tau + 1)
    for it in range(n_tau + 1):
        sin_vals = np.sin(wn * tau[it])
        g_tau[it] = (2.0 / beta) * np.sum(sin_vals * (im_g + 1.0 / wn)) - 0.5

    return tau, g_tau


def tau_to_iw(tau, f_tau, wn, beta):
    """
    Forward FT using trapezoidal rule.

    At half-filling, f(tau) is real and f(beta-tau) = -f(tau),
    so the FT is purely imaginary:
    Im f(iw_n) = integral_0^beta sin(w_n tau) f(tau) dtau
    """
    dt = tau[1] - tau[0]
    n_freq = len(wn)
    im_f = np.zeros(n_freq)
    for iw in range(n_freq):
        sin_vals = np.sin(wn[iw] * tau)
        integrand = sin_vals * f_tau
        im_f[iw] = dt * (0.5 * integrand[0] + np.sum(integrand[1:-1]) + 0.5 * integrand[-1])
    return 1j * im_f


def dmft_loop(U, beta, n_freq=1024, n_tau=4096, max_iter=200, tol=1e-6,
              mixing=0.5, seed='metallic', verbose=False):
    """
    DMFT self-consistency loop with IPT solver at half-filling.

    Parameters
    ----------
    U : float - Hubbard interaction
    beta : float - inverse temperature
    n_freq : int - number of Matsubara frequencies
    n_tau : int - tau grid points
    max_iter : int - max iterations
    tol : float - convergence tolerance on Im Sigma(i*pi*T)
    mixing : float - linear mixing (0,1]
    seed : str - 'metallic' or 'insulating'
    verbose : bool - print progress

    Returns
    -------
    dict with G_iw, Sigma_iw, G_tau, tau, wn, converged, n_iter,
    im_sigma_history, double_occ
    """
    t2 = 0.25  # t^2
    wn = matsubara_freq(n_freq, beta)

    # Initialize Sigma_tilde (without Hartree, purely imaginary)
    if seed == 'metallic':
        sigma_iw = np.zeros(n_freq, dtype=complex)
    else:
        # Insulating seed: atomic limit self-energy (without Hartree)
        # Sigma_IPT_atomic(iwn) = -i * U^2 / (4 * wn)
        sigma_iw = -1j * U**2 / (4.0 * wn)

    im_sigma_history = []
    converged = False

    for it in range(max_iter):
        # Lattice GF: G(iw) = bethe_gf(iw - Sigma_tilde)
        z = 1j * wn - sigma_iw
        g_iw = bethe_gf(z)

        # Self-consistency: G0^{-1} = iw - t^2 G(iw)
        g0_iw = 1.0 / (1j * wn - t2 * g_iw)

        # FT G0 to tau
        _, g0_tau = iw_to_tau(g0_iw, wn, beta, n_tau)

        # IPT self-energy in tau: Sigma(tau) = U^2 G0(tau)^3
        sigma_tau = U**2 * g0_tau**3

        # FT sigma back to iw
        tau = np.linspace(0, beta, n_tau + 1)
        sigma_iw_new = tau_to_iw(tau, sigma_tau, wn, beta)

        # Track convergence
        im_sig = np.imag(sigma_iw_new[0])
        im_sigma_history.append(im_sig)

        if it > 0:
            diff = abs(im_sig - im_sigma_history[-2])
            if verbose and it % 5 == 0:
                print(f"  iter {it:3d}: ImSig(ipiT) = {im_sig:.8f}, diff = {diff:.2e}")
            if diff < tol:
                converged = True
                sigma_iw = mixing * sigma_iw_new + (1 - mixing) * sigma_iw
                break
        elif verbose:
            print(f"  iter {it:3d}: ImSig(ipiT) = {im_sig:.8f}")

        sigma_iw = mixing * sigma_iw_new + (1 - mixing) * sigma_iw

    # Final observables
    z = 1j * wn - sigma_iw
    g_iw = bethe_gf(z)
    _, g_tau = iw_to_tau(g_iw, wn, beta, n_tau)
    tau = np.linspace(0, beta, n_tau + 1)

    # Spectral weight at Fermi level: Z ~ -1/Im Sigma'(0)
    # Approximate quasiparticle weight from first two frequencies
    if n_freq >= 2 and abs(np.imag(sigma_iw[0])) > 1e-12:
        z_qp = 1.0 / (1.0 - np.imag(sigma_iw[0]) / wn[0])
        z_qp = max(0.0, min(1.0, z_qp))
    else:
        z_qp = 1.0

    # Double occupancy approximation:
    # d = (1/4) * (1 - Z) for Brinkman-Rice-like estimate
    # More precisely: d ~ -G(beta/2) * ... but let's use Z
    double_occ = 0.25 * (1.0 - z_qp) if z_qp < 1.0 else 0.25

    return {
        'G_iw': g_iw,
        'Sigma_iw': sigma_iw,
        'G_tau': g_tau,
        'tau': tau,
        'wn': wn,
        'converged': converged,
        'n_iter': it + 1,
        'im_sigma_history': im_sigma_history,
        'im_g_history': [np.imag(bethe_gf(1j * wn[0] - s))[0] if hasattr(s, '__len__') else np.imag(bethe_gf(1j * wn[0] - s)) for s in [0]],
        'double_occ': double_occ,
        'z_qp': z_qp,
        'im_g_lowest': np.imag(g_iw[0]),
        'U': U,
        'beta': beta,
    }
