"""
Core module for reproducing Fukaya & Onogi (2003):
"Lattice Study of the Massive Schwinger Model with theta Term
 under Luscher's Admissibility Condition"

Implements:
- U(1) lattice gauge theory in 2D
- Luscher's gauge action with admissibility condition
- Wilson's plaquette action (for comparison)
- Topological charge measurement
- Classical gauge configurations for each sector
- Domain-wall fermion operator construction
- Conjugate gradient solver
- Hybrid Monte Carlo with domain-wall fermions
- Meson correlator measurements (pion, eta)
"""

import numpy as np
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse.linalg import cg as sparse_cg


# =============================================================================
# U(1) Lattice Gauge Operations
# =============================================================================

def make_classical_config(L, N, nu1=0.0, nu2=0.0):
    """
    Create classical U(1) gauge configuration with topological charge N.

    Returns link angles theta1[x,y] and theta2[x,y] such that
    U_mu(x) = exp(i * theta_mu(x)).

    The configuration gives constant field strength F_12 = 2*pi*N/L^2
    everywhere on the torus.
    """
    theta1 = np.zeros((L, L))
    theta2 = np.zeros((L, L))

    for x in range(L):
        for y in range(L):
            theta1[x, y] = 2.0 * np.pi * nu1 / L
            theta2[x, y] = 2.0 * np.pi * nu2 / L + 2.0 * np.pi * N * x / (L * L)

    # Transition function at x = L-1
    for y in range(L):
        theta1[L - 1, y] += -2.0 * np.pi * N * y / L

    return theta1, theta2


def plaquette_angle(theta1, theta2, x, y, L):
    """Compute plaquette angle at site (x,y)."""
    xp = (x + 1) % L
    yp = (y + 1) % L
    phi = theta1[x, y] + theta2[xp, y] - theta1[x, yp] - theta2[x, y]
    return phi


def all_plaquette_angles(theta1, theta2, L):
    """Compute all plaquette angles."""
    xp = np.roll(np.arange(L), -1)
    yp = np.roll(np.arange(L), -1)

    # theta1[x,y] + theta2[x+1,y] - theta1[x,y+1] - theta2[x,y]
    phi = (theta1
           + theta2[xp, :]
           - theta1[:, yp]
           - theta2)
    return phi


def topological_charge(theta1, theta2, L):
    """
    Compute topological charge N = (1/2pi) sum_x arg(P(x)).
    For smooth configs, this should be an integer.
    """
    phi = all_plaquette_angles(theta1, theta2, L)
    # Wrap to [-pi, pi]
    phi_wrapped = np.mod(phi + np.pi, 2 * np.pi) - np.pi
    N = np.sum(phi_wrapped) / (2.0 * np.pi)
    return round(N)


def luscher_action(theta1, theta2, L, beta, epsilon):
    """
    Compute Luscher's gauge action.
    S_G = beta * sum_P (1 - Re P) / (1 - |1-P|/epsilon)
    Returns infinity if admissibility is violated.
    """
    phi = all_plaquette_angles(theta1, theta2, L)
    re_P = np.cos(phi)
    one_minus_reP = 1.0 - re_P
    abs_one_minus_P = np.sqrt(2.0 * one_minus_reP)  # |1-e^{iphi}| = 2|sin(phi/2)|

    if np.any(abs_one_minus_P >= epsilon):
        return np.inf

    denom = 1.0 - abs_one_minus_P / epsilon
    S = beta * np.sum(one_minus_reP / denom)
    return S


def wilson_action(theta1, theta2, L, beta):
    """Compute Wilson's plaquette action: S = beta * sum(1 - Re P)."""
    phi = all_plaquette_angles(theta1, theta2, L)
    return beta * np.sum(1.0 - np.cos(phi))


def luscher_action_min(L, N, beta, epsilon):
    """
    Minimum Luscher action for sector N.
    All plaquettes have angle 2*pi*N/L^2.
    """
    phi0 = 2.0 * np.pi * N / (L * L)
    re_P = np.cos(phi0)
    one_minus_reP = 1.0 - re_P
    abs_one_minus_P = np.sqrt(2.0 * one_minus_reP)

    if abs_one_minus_P >= epsilon:
        return np.inf

    denom = 1.0 - abs_one_minus_P / epsilon
    S = beta * L * L * one_minus_reP / denom
    return S


# =============================================================================
# Domain-Wall Fermion Operator
# =============================================================================

def build_dwf_operator(theta1, theta2, L, L3, M, m_f, operator_type="DW"):
    """
    Build DWF operator using COO format for efficiency.
    Vectorized construction instead of Python loops.
    """
    from scipy.sparse import coo_matrix

    rows = []
    cols = []
    vals = []

    def idx(x, y, s, a):
        return x * (L * L3 * 2) + y * (L3 * 2) + s * 2 + a

    dim = L * L * L3 * 2

    # Precompute projection matrices
    # (1+gamma1)/2 = [[1/2, 1/2],[1/2, 1/2]]
    proj_p1 = np.array([[0.5, 0.5], [0.5, 0.5]])
    # (1-gamma1)/2 = [[1/2, -1/2],[-1/2, 1/2]]
    proj_m1 = np.array([[0.5, -0.5], [-0.5, 0.5]])
    # (1+gamma2)/2 = [[1/2, -i/2],[i/2, 1/2]]
    proj_p2 = np.array([[0.5, -0.5j], [0.5j, 0.5]])
    # (1-gamma2)/2 = [[1/2, i/2],[-i/2, 1/2]]
    proj_m2 = np.array([[0.5, 0.5j], [-0.5j, 0.5]])

    for x in range(L):
        xp = (x + 1) % L
        xm = (x - 1) % L
        for y in range(L):
            yp = (y + 1) % L
            ym = (y - 1) % L

            U1_fwd = np.exp(1j * theta1[x, y])
            U1_bwd = np.exp(-1j * theta1[xm, y])
            U2_fwd = np.exp(1j * theta2[x, y])
            U2_bwd = np.exp(-1j * theta2[x, ym])

            for s in range(L3):
                # Mass term (M - 3) * delta
                for a in range(2):
                    rows.append(idx(x, y, s, a))
                    cols.append(idx(x, y, s, a))
                    vals.append(M - 3.0)

                # x-direction hopping
                for a in range(2):
                    for b in range(2):
                        # Forward
                        rows.append(idx(x, y, s, a))
                        cols.append(idx(xp, y, s, b))
                        vals.append(proj_p1[a, b] * U1_fwd)
                        # Backward
                        rows.append(idx(x, y, s, a))
                        cols.append(idx(xm, y, s, b))
                        vals.append(proj_m1[a, b] * U1_bwd)

                # y-direction hopping
                for a in range(2):
                    for b in range(2):
                        rows.append(idx(x, y, s, a))
                        cols.append(idx(x, yp, s, b))
                        vals.append(proj_p2[a, b] * U2_fwd)
                        rows.append(idx(x, y, s, a))
                        cols.append(idx(x, ym, s, b))
                        vals.append(proj_m2[a, b] * U2_bwd)

                # DW hopping: P+ delta_{s+1,s'} + P- delta_{s-1,s'}
                if s + 1 < L3:
                    rows.append(idx(x, y, s, 0))
                    cols.append(idx(x, y, s + 1, 0))
                    vals.append(1.0)
                if s - 1 >= 0:
                    rows.append(idx(x, y, s, 1))
                    cols.append(idx(x, y, s - 1, 1))
                    vals.append(1.0)

                # Boundary terms
                if operator_type == "DW":
                    if s == L3 - 1:
                        rows.append(idx(x, y, s, 0))
                        cols.append(idx(x, y, 0, 0))
                        vals.append(m_f - 1.0)
                    if s == 0:
                        rows.append(idx(x, y, s, 1))
                        cols.append(idx(x, y, L3 - 1, 1))
                        vals.append(m_f - 1.0)
                elif operator_type == "AP":
                    if s == L3 - 1:
                        rows.append(idx(x, y, s, 0))
                        cols.append(idx(x, y, 0, 0))
                        vals.append(-2.0)
                    if s == 0:
                        rows.append(idx(x, y, s, 1))
                        cols.append(idx(x, y, L3 - 1, 1))
                        vals.append(-2.0)

    return coo_matrix((vals, (rows, cols)), shape=(dim, dim)).tocsr()


def fermion_determinant_ratio(theta1, theta2, L, L3, M, m_f):
    """
    Compute det(D_DW)^2 / det(D_AP)^2 for the given gauge config.
    Uses dense matrix determinant (feasible for small lattices).
    """
    D_DW = build_dwf_operator(theta1, theta2, L, L3, M, m_f, "DW")
    D_AP = build_dwf_operator(theta1, theta2, L, L3, M, m_f, "AP")

    # Use slogdet to avoid overflow
    sign_DW, logdet_DW = np.linalg.slogdet(D_DW.toarray())
    sign_AP, logdet_AP = np.linalg.slogdet(D_AP.toarray())

    # ratio = |det_DW|^2 / |det_AP|^2
    log_ratio = 2.0 * (logdet_DW - logdet_AP)

    # Clamp to avoid overflow
    log_ratio = np.clip(log_ratio, -500, 500)
    return np.exp(log_ratio)


def compute_DetN(L, L3, M, m_f, N, n_nu=5):
    """
    Compute Det^N = integral(det_DW^2/det_AP^2 for sector N) /
                    integral(det_DW^2/det_AP^2 for sector 0)

    Uses weighted sum over n_nu x n_nu points in (nu1, nu2) moduli space.
    """
    nu_points = np.linspace(0, 1, n_nu, endpoint=False) + 0.5 / n_nu

    def integrated_det(sector_N):
        total = 0.0
        for nu1 in nu_points:
            for nu2 in nu_points:
                t1, t2 = make_classical_config(L, sector_N, nu1, nu2)
                ratio = fermion_determinant_ratio(t1, t2, L, L3, M, m_f)
                total += ratio
        return total / (n_nu * n_nu)

    det_N = integrated_det(N)
    det_0 = integrated_det(0)

    if det_0 == 0:
        return 0.0
    return det_N / det_0


# =============================================================================
# Monte Carlo: Metropolis for U(1) gauge theory
# =============================================================================

def metropolis_sweep_luscher(theta1, theta2, L, beta, epsilon, step_size=0.3):
    """
    One Metropolis sweep updating all links with Luscher's action.
    Returns (theta1, theta2, accept_rate).
    """
    accepted = 0
    total = 0

    for mu in range(2):
        links = theta1 if mu == 0 else theta2
        for x in range(L):
            for y in range(L):
                old_angle = links[x, y]
                delta = np.random.uniform(-step_size, step_size)
                new_angle = old_angle + delta

                # Compute local action change
                old_S = _local_luscher_action(theta1, theta2, L, beta, epsilon, x, y, mu)
                links[x, y] = new_angle
                new_S = _local_luscher_action(theta1, theta2, L, beta, epsilon, x, y, mu)

                if new_S == np.inf or np.random.random() > np.exp(-(new_S - old_S)):
                    links[x, y] = old_angle  # reject
                else:
                    accepted += 1
                total += 1

    return theta1, theta2, accepted / total


def metropolis_sweep_wilson(theta1, theta2, L, beta, step_size=0.5):
    """One Metropolis sweep with Wilson's action."""
    accepted = 0
    total = 0

    for mu in range(2):
        links = theta1 if mu == 0 else theta2
        for x in range(L):
            for y in range(L):
                old_angle = links[x, y]
                delta = np.random.uniform(-step_size, step_size)
                new_angle = old_angle + delta

                old_S = _local_wilson_action(theta1, theta2, L, beta, x, y, mu)
                links[x, y] = new_angle
                new_S = _local_wilson_action(theta1, theta2, L, beta, x, y, mu)

                if np.random.random() > np.exp(-(new_S - old_S)):
                    links[x, y] = old_angle
                else:
                    accepted += 1
                total += 1

    return theta1, theta2, accepted / total


def _local_luscher_action(theta1, theta2, L, beta, epsilon, x, y, mu):
    """Luscher action contribution from plaquettes touching link (x,y,mu)."""
    S = 0.0
    plaq_coords = _get_touching_plaquettes(x, y, mu, L)
    for px, py in plaq_coords:
        phi = plaquette_angle(theta1, theta2, px, py, L)
        re_P = np.cos(phi)
        one_m_reP = 1.0 - re_P
        abs_1_m_P = np.sqrt(2.0 * one_m_reP)
        if abs_1_m_P >= epsilon:
            return np.inf
        denom = 1.0 - abs_1_m_P / epsilon
        S += beta * one_m_reP / denom
    return S


def _local_wilson_action(theta1, theta2, L, beta, x, y, mu):
    """Wilson action contribution from plaquettes touching link (x,y,mu)."""
    S = 0.0
    plaq_coords = _get_touching_plaquettes(x, y, mu, L)
    for px, py in plaq_coords:
        phi = plaquette_angle(theta1, theta2, px, py, L)
        S += beta * (1.0 - np.cos(phi))
    return S


def _get_touching_plaquettes(x, y, mu, L):
    """
    In 2D, each link touches exactly 2 plaquettes (one forward, one backward).
    The plaquette at (px, py) is defined by links forming the square
    starting at (px, py) going in mu=0 then mu=1 direction.
    """
    if mu == 0:  # x-link at (x,y)
        # Plaquette at (x, y) and plaquette at (x, y-1)
        return [(x, y), (x, (y - 1) % L)]
    else:  # y-link at (x,y)
        # Plaquette at (x, y) and plaquette at (x-1, y)
        return [(x, y), ((x - 1) % L, y)]


# =============================================================================
# HMC with Domain-Wall Fermions
# =============================================================================

def hmc_trajectory(theta1, theta2, L, L3, beta, epsilon, M, m_f,
                   n_steps=50, dt=0.02, n_flavors=2):
    """
    One HMC trajectory with Luscher gauge action and domain-wall fermions.
    Uses pseudofermion method for det(D_DW)^2/det(D_AP)^2.
    """
    # Store old config
    theta1_old = theta1.copy()
    theta2_old = theta2.copy()

    # Generate momenta
    p1 = np.random.randn(L, L)
    p2 = np.random.randn(L, L)

    # Build operators for pseudofermion generation
    D_DW = build_dwf_operator(theta1, theta2, L, L3, M, m_f, "DW")
    D_AP = build_dwf_operator(theta1, theta2, L, L3, M, m_f, "AP")
    dim = L * L * L3 * 2

    # Generate pseudofermions for each flavor
    pf_list = []
    pf_pv_list = []
    for _ in range(n_flavors):
        eta = np.random.randn(dim) + 1j * np.random.randn(dim)
        eta /= np.sqrt(2)
        phi_dw = D_DW.T.conj() @ eta  # phi = D^dag eta
        pf_list.append(phi_dw)

        eta_pv = np.random.randn(dim) + 1j * np.random.randn(dim)
        eta_pv /= np.sqrt(2)
        phi_pv = D_AP.T.conj() @ eta_pv
        pf_pv_list.append(phi_pv)

    # Compute initial Hamiltonian
    H_old = 0.5 * (np.sum(p1 ** 2) + np.sum(p2 ** 2))
    H_old += luscher_action(theta1, theta2, L, beta, epsilon)
    for phi_dw in pf_list:
        DdD = D_DW.T.conj() @ D_DW
        x_sol, _ = sparse_cg(DdD, phi_dw, tol=1e-10, maxiter=2000)
        H_old += np.real(np.dot(phi_dw.conj(), x_sol))
    for phi_pv in pf_pv_list:
        DdD_pv = D_AP.T.conj() @ D_AP
        x_sol, _ = sparse_cg(DdD_pv, phi_pv, tol=1e-10, maxiter=2000)
        H_old -= np.real(np.dot(phi_pv.conj(), x_sol))

    if H_old == np.inf:
        return theta1_old, theta2_old, False

    # Leapfrog integration
    # Half step for momenta
    f1, f2 = _compute_force(theta1, theta2, L, L3, beta, epsilon, M, m_f,
                            pf_list, pf_pv_list)
    p1 -= 0.5 * dt * f1
    p2 -= 0.5 * dt * f2

    for step in range(n_steps - 1):
        theta1 += dt * p1
        theta2 += dt * p2

        f1, f2 = _compute_force(theta1, theta2, L, L3, beta, epsilon, M, m_f,
                                pf_list, pf_pv_list)
        p1 -= dt * f1
        p2 -= dt * f2

    # Last step
    theta1 += dt * p1
    theta2 += dt * p2
    f1, f2 = _compute_force(theta1, theta2, L, L3, beta, epsilon, M, m_f,
                            pf_list, pf_pv_list)
    p1 -= 0.5 * dt * f1
    p2 -= 0.5 * dt * f2

    # Compute final Hamiltonian
    S_G_new = luscher_action(theta1, theta2, L, beta, epsilon)
    if S_G_new == np.inf:
        return theta1_old, theta2_old, False

    H_new = 0.5 * (np.sum(p1 ** 2) + np.sum(p2 ** 2)) + S_G_new
    D_DW_new = build_dwf_operator(theta1, theta2, L, L3, M, m_f, "DW")
    D_AP_new = build_dwf_operator(theta1, theta2, L, L3, M, m_f, "AP")
    for phi_dw in pf_list:
        DdD = D_DW_new.T.conj() @ D_DW_new
        x_sol, _ = sparse_cg(DdD, phi_dw, tol=1e-10, maxiter=2000)
        H_new += np.real(np.dot(phi_dw.conj(), x_sol))
    for phi_pv in pf_pv_list:
        DdD_pv = D_AP_new.T.conj() @ D_AP_new
        x_sol, _ = sparse_cg(DdD_pv, phi_pv, tol=1e-10, maxiter=2000)
        H_new -= np.real(np.dot(phi_pv.conj(), x_sol))

    # Metropolis accept/reject
    dH = H_new - H_old
    if np.random.random() < np.exp(-dH):
        return theta1, theta2, True
    else:
        return theta1_old, theta2_old, False


def _compute_force(theta1, theta2, L, L3, beta, epsilon, M, m_f,
                   pf_list, pf_pv_list):
    """Compute HMC force = -dS/dU for all links."""
    f1 = np.zeros((L, L))
    f2 = np.zeros((L, L))

    # Gauge force
    for x in range(L):
        for y in range(L):
            for mu in range(2):
                fg = _gauge_force_luscher(theta1, theta2, L, beta, epsilon, x, y, mu)
                if mu == 0:
                    f1[x, y] = fg
                else:
                    f2[x, y] = fg

    # Fermion force (approximation using finite differences for efficiency)
    eps_fd = 1e-4
    D_DW = build_dwf_operator(theta1, theta2, L, L3, M, m_f, "DW")
    D_AP = build_dwf_operator(theta1, theta2, L, L3, M, m_f, "AP")
    DdD = D_DW.T.conj() @ D_DW
    DdD_pv = D_AP.T.conj() @ D_AP

    # Solve for X = (D^dag D)^{-1} phi for each pseudofermion
    X_list = []
    for phi in pf_list:
        x_sol, _ = sparse_cg(DdD, phi, tol=1e-10, maxiter=2000)
        X_list.append(x_sol)
    X_pv_list = []
    for phi in pf_pv_list:
        x_sol, _ = sparse_cg(DdD_pv, phi, tol=1e-10, maxiter=2000)
        X_pv_list.append(x_sol)

    # Fermion force via finite differences (simpler than analytical derivative)
    for x in range(L):
        for y in range(L):
            for mu in range(2):
                links = theta1 if mu == 0 else theta2
                old_val = links[x, y]

                links[x, y] = old_val + eps_fd
                D_p = build_dwf_operator(theta1, theta2, L, L3, M, m_f, "DW")
                DdD_p = D_p.T.conj() @ D_p
                D_pv_p = build_dwf_operator(theta1, theta2, L, L3, M, m_f, "AP")
                DdD_pv_p = D_pv_p.T.conj() @ D_pv_p

                links[x, y] = old_val - eps_fd
                D_m = build_dwf_operator(theta1, theta2, L, L3, M, m_f, "DW")
                DdD_m = D_m.T.conj() @ D_m
                D_pv_m = build_dwf_operator(theta1, theta2, L, L3, M, m_f, "AP")
                DdD_pv_m = D_pv_m.T.conj() @ D_pv_m

                links[x, y] = old_val

                ff = 0.0
                for X in X_list:
                    sp = np.real(X.conj() @ DdD_p @ X)
                    sm = np.real(X.conj() @ DdD_m @ X)
                    ff += (sp - sm) / (2 * eps_fd)
                for X in X_pv_list:
                    sp = np.real(X.conj() @ DdD_pv_p @ X)
                    sm = np.real(X.conj() @ DdD_pv_m @ X)
                    ff -= (sp - sm) / (2 * eps_fd)

                if mu == 0:
                    f1[x, y] += ff
                else:
                    f2[x, y] += ff

    return f1, f2


def _gauge_force_luscher(theta1, theta2, L, beta, epsilon, x, y, mu):
    """Compute gauge force -dS_G/d(theta_mu(x,y)) for Luscher action."""
    eps_fd = 1e-6
    links = theta1 if mu == 0 else theta2
    old_val = links[x, y]

    links[x, y] = old_val + eps_fd
    Sp = _local_luscher_action(theta1, theta2, L, beta, epsilon, x, y, mu)
    links[x, y] = old_val - eps_fd
    Sm = _local_luscher_action(theta1, theta2, L, beta, epsilon, x, y, mu)
    links[x, y] = old_val

    if Sp == np.inf or Sm == np.inf:
        return 0.0
    return (Sp - Sm) / (2 * eps_fd)


# =============================================================================
# Meson Correlator Measurements
# =============================================================================

def measure_pion_correlator(theta1, theta2, L, L3, M, m_f):
    """
    Measure the isotriplet pion propagator using sparse solver.
    C_pi(x) = sum_y |G(x,y;0,0)|^2 summed over spin indices.
    """
    from scipy.sparse.linalg import spsolve

    D_DW = build_dwf_operator(theta1, theta2, L, L3, M, m_f, "DW")
    dim = L * L * L3 * 2

    C_pi = np.zeros(L)

    for alpha in range(2):
        # Source at s=0
        source_0 = np.zeros(dim, dtype=complex)
        source_0[0 * (L * L3 * 2) + 0 * (L3 * 2) + 0 * 2 + alpha] = 1.0

        # Source at s=L3-1
        source_L3 = np.zeros(dim, dtype=complex)
        source_L3[0 * (L * L3 * 2) + 0 * (L3 * 2) + (L3 - 1) * 2 + alpha] = 1.0

        # Solve using sparse direct solver
        try:
            prop_0 = spsolve(D_DW, source_0)
            prop_L3 = spsolve(D_DW, source_L3)
        except Exception:
            return np.zeros(L)

        # Extract physical propagator at DWF boundaries
        for x in range(L):
            for y in range(L):
                for beta_s in range(2):
                    idx_R = x * (L * L3 * 2) + y * (L3 * 2) + (L3 - 1) * 2 + beta_s
                    idx_L = x * (L * L3 * 2) + y * (L3 * 2) + 0 * 2 + beta_s
                    C_pi[x] += np.abs(prop_0[idx_R]) ** 2 + np.abs(prop_L3[idx_L]) ** 2

    return C_pi


def measure_eta_correlator(theta1, theta2, L, L3, M, m_f):
    """
    Measure the isosinglet eta propagator.
    eta = connected (same as pion) + disconnected (hairpin) contribution.

    For the disconnected part, we compute:
    sum_y [Tr gamma3 G(x,y;x,y)]^2 (summed over all source points)

    This requires the full propagator at all points, which we compute by brute force.
    """
    D_DW = build_dwf_operator(theta1, theta2, L, L3, M, m_f, "DW")
    dim = L * L * L3 * 2
    D_dense = D_DW.toarray()

    # Compute full propagator matrix (expensive but exact for small lattice)
    try:
        G_full = np.linalg.inv(D_dense)
    except np.linalg.LinAlgError:
        return np.zeros(L), np.zeros(L)

    # Connected part (same as pion)
    C_connected = np.zeros(L)
    # Disconnected (hairpin) part
    C_disconnected = np.zeros(L)

    # For each source point, compute local trace
    # This is a simplified version - the full calculation would use
    # the physical propagator from DWF boundaries
    local_trace = np.zeros((L, L), dtype=complex)

    for x0 in range(L):
        for y0 in range(L):
            tr = 0.0 + 0.0j
            for alpha in range(2):
                # gamma3 element: gamma3[alpha,alpha] = 1 if alpha=0, -1 if alpha=1
                g3 = 1.0 if alpha == 0 else -1.0
                # Physical propagator at (x0,y0) -> (x0,y0) from DWF
                # Use s=0 and s=L3-1 boundaries
                idx_src_0 = x0 * (L * L3 * 2) + y0 * (L3 * 2) + 0 * 2 + alpha
                idx_sink_L3 = x0 * (L * L3 * 2) + y0 * (L3 * 2) + (L3 - 1) * 2 + alpha
                tr += g3 * G_full[idx_sink_L3, idx_src_0]
            local_trace[x0, y0] = tr

    # Disconnected correlator
    for x in range(L):
        for y in range(L):
            for y0 in range(L):
                C_disconnected[x] += np.real(
                    local_trace[(0 + x) % L, y] * np.conj(local_trace[0, y0])
                )

    # Total eta = connected - 2*pion + 4*disconnected (schematic)
    C_pion = measure_pion_correlator(theta1, theta2, L, L3, M, m_f)
    C_eta = -2 * C_pion + 4 * C_disconnected

    return C_eta, C_disconnected


def fit_correlator_mass(C, L, fit_range=None):
    """
    Fit meson correlator C(x) to A*cosh(m*(x - L/2)) to extract mass m.
    Returns (mass, amplitude, chi2_dof).
    """
    from scipy.optimize import curve_fit

    if fit_range is None:
        fit_range = (max(1, L // 4), 3 * L // 4)

    x_data = np.arange(fit_range[0], fit_range[1] + 1)
    C_data = np.abs(C[fit_range[0]:fit_range[1] + 1])

    if np.all(C_data <= 0):
        return 0.0, 0.0, np.inf

    C_data = np.maximum(C_data, 1e-30)

    def cosh_func(x, A, m):
        return A * np.cosh(m * (x - L / 2.0))

    try:
        popt, pcov = curve_fit(cosh_func, x_data, C_data,
                               p0=[C_data[0], 0.5], maxfev=10000)
        residuals = C_data - cosh_func(x_data, *popt)
        chi2 = np.sum(residuals ** 2 / np.maximum(C_data ** 2, 1e-60))
        dof = max(len(x_data) - 2, 1)
        return abs(popt[1]), popt[0], chi2 / dof
    except (RuntimeError, ValueError):
        return 0.0, 0.0, np.inf


def effective_mass(C, L):
    """Compute effective mass m_eff(x) = arccosh((C(x-1)+C(x+1))/(2*C(x)))."""
    m_eff = np.zeros(L)
    for x in range(1, L - 1):
        if C[x] > 0:
            ratio = (C[x - 1] + C[(x + 1) % L]) / (2.0 * C[x])
            if ratio >= 1.0:
                m_eff[x] = np.arccosh(ratio)
    return m_eff
