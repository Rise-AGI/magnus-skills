"""
Core module for reproducing Hunana et al. (2019):
"An Introductory Guide to Fluid Models with Anisotropic Temperatures
Part 2 - Kinetic Theory, Pade Approximants and Landau Fluid Closures"
arXiv: 1901.09360v2

Contains:
- Exact plasma dispersion function Z(zeta) and response function R(zeta)
- Pade approximants R_{n,n'}(zeta) for n=1..8
- Dispersion relation solvers for ion-acoustic and Langmuir modes
"""

import numpy as np
from scipy.special import erfi, wofz

# ---------------------------------------------------------------------------
# Exact plasma dispersion function and response function
# ---------------------------------------------------------------------------

def Z_exact(zeta):
    """Plasma dispersion function Z(zeta) = i*sqrt(pi)*exp(-zeta^2)*[1+erf(i*zeta)]

    Uses the Faddeeva function w(z) = exp(-z^2)*erfc(-iz) for numerical stability:
    Z(zeta) = i*sqrt(pi)*w(zeta)
    """
    return 1j * np.sqrt(np.pi) * wofz(zeta)


def R_exact(zeta):
    """Plasma response function R(zeta) = 1 + zeta*Z(zeta)"""
    return 1.0 + zeta * Z_exact(zeta)


# ---------------------------------------------------------------------------
# Pade approximants of R(zeta)
# From Eqs. (165), (186)-(188), (221)-(225) and higher-order formulas
# ---------------------------------------------------------------------------

def R_1(zeta):
    """1-pole approximant R_1(zeta), Eq. (165)"""
    return 1.0 / (1.0 - 1j * np.sqrt(np.pi) * zeta)


def R_2_0(zeta):
    """2-pole approximant R_{2,0}(zeta), Eq. (175) via R=1+zeta*Z"""
    z = zeta
    sp = np.sqrt(np.pi)
    num = 1j * sp + 2 * z
    den = 1 - 1j * sp * z - 2 * z**2
    return 1.0 + z * num / den


def R_3_0(zeta):
    """3-pole approximant R_{3,0}(zeta), Eq. (186)"""
    z = zeta
    sp = np.sqrt(np.pi)
    pi_val = np.pi
    a1 = -1j * sp * (pi_val - 3) / (4 - pi_val)
    b1 = -1j * sp / (4 - pi_val)
    b2 = -(3 * pi_val - 8) / (4 - pi_val)
    b3 = 2j * sp * (pi_val - 3) / (4 - pi_val)
    num = 1.0 + a1 * z
    den = 1.0 + b1 * z + b2 * z**2 + b3 * z**3
    return num / den


def R_3_1(zeta):
    """3-pole approximant R_{3,1}(zeta), Eq. (187)"""
    z = zeta
    sp = np.sqrt(np.pi)
    pi_val = np.pi
    a1 = -1j * (4 - pi_val) / sp
    b1 = -4j / sp
    b2 = -2.0
    b3 = 2j * (4 - pi_val) / sp
    num = 1.0 + a1 * z
    den = 1.0 + b1 * z + b2 * z**2 + b3 * z**3
    return num / den


def R_3_2(zeta):
    """3-pole approximant R_{3,2}(zeta), Eq. (188)"""
    z = zeta
    sp = np.sqrt(np.pi)
    a1 = -1j * sp / 2
    b1 = -3j * sp / 2
    b2 = -2.0
    b3 = 1j * sp
    num = 1.0 + a1 * z
    den = 1.0 + b1 * z + b2 * z**2 + b3 * z**3
    return num / den


def R_4_0(zeta):
    """4-pole approximant R_{4,0}(zeta), Eq. (221)"""
    z = zeta
    sp = np.sqrt(np.pi)
    pi_val = np.pi
    D = 6 * pi_val**2 - 29 * pi_val + 32
    a1 = 1j * sp / 2 * (12 * pi_val**2 - 67 * pi_val + 92) / D
    a2 = -(9 * pi_val**2 - 69 * pi_val + 128) / (6 * D)
    b1 = -1j * sp / 2 * (9 * pi_val - 28) / D
    b2 = (36 * pi_val**2 - 195 * pi_val + 256) / (6 * D)
    b3 = -1j * sp * (33 * pi_val - 104) / (6 * D)
    b4 = (9 * pi_val**2 - 69 * pi_val + 128) / (3 * D)
    num = 1.0 + a1 * z + a2 * z**2
    den = 1.0 + b1 * z + b2 * z**2 + b3 * z**3 + b4 * z**4
    return num / den


def R_4_1(zeta):
    """4-pole approximant R_{4,1}(zeta), Eq. (222)"""
    z = zeta
    sp = np.sqrt(np.pi)
    pi_val = np.pi
    D = 16 - 5 * pi_val
    a1 = -1j * sp / 3 * (9 * pi_val - 28) / D
    a2 = -(6 * pi_val**2 - 29 * pi_val + 32) / (3 * D)
    b1 = -1j * 2 * sp / 3 * (10 - 3 * pi_val) / D
    b2 = -(21 * pi_val - 64) / (3 * D)
    b3 = 1j * 2 * sp / 3 * (9 * pi_val - 28) / D
    b4 = 2 * (6 * pi_val**2 - 29 * pi_val + 32) / (3 * D)
    num = 1.0 + a1 * z + a2 * z**2
    den = 1.0 + b1 * z + b2 * z**2 + b3 * z**3 + b4 * z**4
    return num / den


def R_4_2(zeta):
    """4-pole approximant R_{4,2}(zeta), Eq. (223)"""
    z = zeta
    sp = np.sqrt(np.pi)
    pi_val = np.pi
    D = 3 * pi_val - 8
    a1 = -1j * sp * (10 - 3 * pi_val) / D
    a2 = -(16 - 5 * pi_val) / D
    b1 = -1j * sp * 2 / D
    b2 = -(32 - 9 * pi_val) / D
    b3 = 1j * sp * 2 * (10 - 3 * pi_val) / D
    b4 = 2 * (16 - 5 * pi_val) / D
    num = 1.0 + a1 * z + a2 * z**2
    den = 1.0 + b1 * z + b2 * z**2 + b3 * z**3 + b4 * z**4
    return num / den


def R_4_3(zeta):
    """4-pole approximant R_{4,3}(zeta), Eq. (224)/(226)
    Hammett & Perkins (1990) approximant."""
    z = zeta
    sp = np.sqrt(np.pi)
    pi_val = np.pi
    num = 4 - 2j * sp * z - (3 * pi_val - 8) * z**2
    den = (4 - 6j * sp * z - (9 * pi_val - 16) * z**2
           + 4j * sp * z**3 + 2 * (3 * pi_val - 8) * z**4)
    return num / den


def R_4_4(zeta):
    """4-pole approximant R_{4,4}(zeta), Eq. (225)"""
    z = zeta
    sp = np.sqrt(np.pi)
    a1 = -1j * sp / 2
    a2 = -2.0 / 3
    b1 = -3j * sp / 2
    b2 = -4.0
    b3 = 1j * sp
    b4 = 4.0 / 3
    # Note: b1 should give the correct sign
    num = 1.0 + a1 * z + a2 * z**2
    den = 1.0 + b1 * z + b2 * z**2 + b3 * z**3 + b4 * z**4
    return num / den


# ---------------------------------------------------------------------------
# 5-pole approximants: constructed by matching power series and asymptotic
# R_{5,n'}(zeta) = (1 + a1*z + a2*z + a3*z^2) / (1 + b1*z + ... + b5*z^5)
# We compute these numerically by solving the matching conditions
# ---------------------------------------------------------------------------

def _compute_5pole_coeffs():
    """Compute all 5-pole Pade approximant coefficients numerically.

    For R_{5,n'}(zeta) = (1 + a1*z + a2*z^2 + a3*z^3) / (1 + b1*z + ... + b5*z^5)

    Power series of R(zeta) for small zeta:
    R(zeta) = 1 + i*sqrt(pi)*z - 2*z^2 - i*sqrt(pi)*z^3 + (4/3)*z^4 + i*sqrt(pi)/2*z^5
              - (8/15)*z^6 - i*sqrt(pi)/6*z^7 + ...

    Asymptotic: R(zeta) ~ -1/(2*z^2) - 3/(4*z^4) - 15/(8*z^6) - ...
    """
    sp = np.sqrt(np.pi)
    # Power series coefficients of R(zeta): R = sum c_k * zeta^k
    c = [1, 1j * sp, -2, -1j * sp, 4.0/3, 1j * sp / 2, -8.0/15, -1j * sp / 6,
         16.0/105, 1j * sp / 24]

    results = {}

    # R_{5,n'} has form (1 + a1*z + a2*z^2 + a3*z^3) / (1 + b1*z + b2*z^2 + b3*z^3 + b4*z^4 + b5*z^5)
    # Asymptotic: leading is a3/b5 * 1/z^2 = -1/2 => b5 = -2*a3
    # Next: (a2*b5 - a3*b4)/b5^2 * 1/z^3 = 0 => b4 = a2*b5/a3 ... nah
    # Let me use the general approach: match power series coefficients

    # For R_{5,0}: match c0..c7 (8 equations for 8 unknowns: a1,a2,a3,b1..b5)
    # plus asymptotic: b5 = -2*a3
    # Total: 9 constraints for 8 unknowns, but asymptotic gives 1 constraint
    # Actually: 3 numerator + 5 denominator = 8 unknowns
    # Asymptotic constraint: b5 = -2*a3 (reduces to 7 unknowns)
    # Power series: match c0..c6 (7 equations, since c0=1 is automatic)
    # This gives R_{5,0}

    # More systematic approach: use the relation
    # (1 + a1*z + a2*z^2 + a3*z^3) = R(z) * (1 + b1*z + ... + b5*z^5)
    # Matching up to order z^(3+5) = z^8

    # For R_{5,n'}: n' additional asymptotic points beyond the basic b5=-2*a3
    # n'=0: b5=-2*a3, match c1..c7
    # n'=1: b5=-2*a3, b4=-2*a2, match c1..c6
    # n'=2: b5=-2*a3, b4=-2*a2, b3=-2*a1+3*a3, match c1..c5
    # etc.

    for n_prime in range(7):  # 0 to 6
        # Number of asymptotic constraints: 1 + n_prime
        # Number of power series constraints: 7 - n_prime (matching c1..c_{7-n_prime})
        # Total free parameters: a1,a2,a3,b1,b2,b3,b4,b5 = 8

        # Build the matching system
        # (1 + a1*z + a2*z^2 + a3*z^3)(1 + c1*z + c2*z^2 + ...) has to equal
        # (1 + b1*z + ... + b5*z^5) when we multiply both sides
        # Actually: num = R * den, so:
        # 1 + a1*z + a2*z^2 + a3*z^3 = (c0 + c1*z + c2*z^2 + ...)(1 + b1*z + ... + b5*z^5)
        # Matching coefficient of z^k for k=0..N

        # This gives: for k <= 3: a_k = sum_{j=0}^{k} c_{k-j} * b_j  (b_0=1, a_0=1)
        # For k > 3: 0 = sum_{j=0}^{min(k,5)} c_{k-j} * b_j

        # Unknowns: a1, a2, a3, b1, b2, b3, b4, b5
        # Let x = [a1, a2, a3, b1, b2, b3, b4, b5]

        # Asymptotic constraints (from the paper's general scheme):
        # b5 = -2*a3
        # b4 = -2*a2 (if n' >= 1)
        # b3 = -2*a1 + 3*a3 (if n' >= 2), from matching 1/z^4 term
        # b2 = 3*a2 - 2 (if n' >= 3), from matching 1/z^5 term
        # b1 = 3*a1 (if n' >= 4)
        # Then fixed a3, a2, a1 for n' >= 5, 6

        # It's cleaner to solve numerically
        # Use: multiply out and match coefficients
        # (1 + sum_{i=1}^3 a_i z^i) = (sum_{k=0}^{M} c_k z^k)(1 + sum_{j=1}^5 b_j z^j)
        # up to order M = 7 - n_prime + 3 = 10 - n_prime

        n_ps = 7 - n_prime  # number of power series terms to match (after z^0)
        n_as = 1 + n_prime  # number of asymptotic constraints

        # Build linear system
        # Variables: [a1, a2, a3, b1, b2, b3, b4, b5] (8 unknowns)
        A_mat = np.zeros((8, 8), dtype=complex)
        rhs = np.zeros(8, dtype=complex)

        eq_idx = 0

        # Power series equations: coefficient of z^k must match
        # For k=1..n_ps+3 (we need enough), but only k=1..max(n_ps,7)
        for k in range(1, n_ps + 4):
            if eq_idx >= 8:
                break
            # Coefficient of z^k in the product (sum c_j z^j)(1 + sum b_i z^i)
            # = c_k + sum_{i=1}^{min(k,5)} b_i * c_{k-i}
            # This should equal: a_k for k <= 3, or 0 for k > 3
            # So: a_k - c_k - sum b_i * c_{k-i} = 0

            row = np.zeros(8, dtype=complex)
            val = c[k] if k < len(c) else 0

            # a_k contribution (for k=1,2,3)
            if 1 <= k <= 3:
                row[k - 1] = 1.0  # a_k

            # -b_i * c_{k-i} contributions
            for i in range(1, min(k, 5) + 1):
                ci = c[k - i] if (k - i) < len(c) else 0
                row[2 + i] = -ci  # b_i is at index 2+i (b1=index 3, ..., b5=index 7)

            A_mat[eq_idx] = row
            rhs[eq_idx] = val
            eq_idx += 1

        # If we have fewer than 8 equations from power series, add asymptotic constraints
        # Asymptotic: b5 = -2*a3 => b5 + 2*a3 = 0
        if eq_idx < 8:
            row = np.zeros(8, dtype=complex)
            row[7] = 1.0  # b5
            row[2] = 2.0  # a3
            A_mat[eq_idx] = row
            rhs[eq_idx] = 0
            eq_idx += 1

        if n_prime >= 1 and eq_idx < 8:
            # b4 = -2*a2
            row = np.zeros(8, dtype=complex)
            row[6] = 1.0  # b4
            row[1] = 2.0  # a2
            A_mat[eq_idx] = row
            rhs[eq_idx] = 0
            eq_idx += 1

        if n_prime >= 2 and eq_idx < 8:
            # From matching 1/z^4: b3 + 2*a1 - 3*a3 = 0
            # (This comes from: (a2*b5 - a3*b4)/b5^2 * ... = 0 and higher terms)
            # Actually let me use: b3 = -2*a1 (from R_{5,2} scheme)
            # Wait, need to be more careful. Let me use the general scheme for 5-pole.
            # For a 5-pole R with form (1+a1*z+a2*z^2+a3*z^3)/(1+b1*z+...+b5*z^5):
            # Asymptotic: a3/b5 = -1/2 => b5=-2*a3 [always]
            # (a2*b5-a3*b4)/b5^2 = 0 => b4=-2*a2 [n'>=1]
            # Next order: involves b3... let me compute
            # (a1*b5^2 - a2*b4*b5 + a3*(b4^2-b3*b5))/b5^3 = 0
            # With b5=-2*a3, b4=-2*a2:
            # (a1*4*a3^2 - a2*(-2*a2)*(-2*a3) + a3*(4*a2^2+2*a3*b3))/(−8a3^3) = 0
            # 4*a1*a3^2 - 4*a2^2*a3 + 4*a2^2*a3 + 2*a3^2*b3 = 0
            # 4*a1*a3^2 + 2*a3^2*b3 = 0
            # b3 = -2*a1
            row = np.zeros(8, dtype=complex)
            row[5] = 1.0  # b3
            row[0] = 2.0  # a1
            A_mat[eq_idx] = row
            rhs[eq_idx] = 0
            eq_idx += 1

        if n_prime >= 3 and eq_idx < 8:
            # Next asymptotic: b2 = 3*a2 - 2
            row = np.zeros(8, dtype=complex)
            row[4] = 1.0   # b2
            row[1] = -3.0   # -3*a2
            A_mat[eq_idx] = row
            rhs[eq_idx] = -2.0
            eq_idx += 1

        if n_prime >= 4 and eq_idx < 8:
            # b1 = 3*a1
            row = np.zeros(8, dtype=complex)
            row[3] = 1.0   # b1
            row[0] = -3.0   # -3*a1
            A_mat[eq_idx] = row
            rhs[eq_idx] = 0
            eq_idx += 1

        if n_prime >= 5 and eq_idx < 8:
            # a3 specific value, from matching deeper asymptotic
            # For R_{4,4}: a2 = -2/3. For 5-pole analog, a3 = ?
            # From the pattern: for R_{5,5}, a3 should follow from
            # the next asymptotic match. Let me use the paper's Table:
            # For R_{5,5}: precision o(zeta), o(zeta^-7)
            # For R_{5,6}: precision o(zeta^0), o(zeta^-8)
            # The matching becomes very specific.
            # Let me just let the linear system handle it.
            pass

        try:
            x = np.linalg.solve(A_mat, rhs)
            a1, a2, a3 = x[0], x[1], x[2]
            b1, b2, b3, b4, b5 = x[3], x[4], x[5], x[6], x[7]
            results[n_prime] = (a1, a2, a3, b1, b2, b3, b4, b5)
        except np.linalg.LinAlgError:
            results[n_prime] = None

    return results


# Cache the 5-pole coefficients
_5pole_cache = None

def _get_5pole():
    global _5pole_cache
    if _5pole_cache is None:
        _5pole_cache = _compute_5pole_coeffs()
    return _5pole_cache


def _R_5_generic(zeta, n_prime):
    """Generic 5-pole approximant R_{5,n'}(zeta)"""
    coeffs = _get_5pole()
    if n_prime not in coeffs or coeffs[n_prime] is None:
        return np.full_like(np.asarray(zeta, dtype=complex), np.nan)
    a1, a2, a3, b1, b2, b3, b4, b5 = coeffs[n_prime]
    z = np.asarray(zeta, dtype=complex)
    num = 1.0 + a1 * z + a2 * z**2 + a3 * z**3
    den = 1.0 + b1 * z + b2 * z**2 + b3 * z**3 + b4 * z**4 + b5 * z**5
    return num / den

def R_5_0(zeta): return _R_5_generic(zeta, 0)
def R_5_1(zeta): return _R_5_generic(zeta, 1)
def R_5_2(zeta): return _R_5_generic(zeta, 2)
def R_5_3(zeta): return _R_5_generic(zeta, 3)
def R_5_4(zeta): return _R_5_generic(zeta, 4)
def R_5_5(zeta): return _R_5_generic(zeta, 5)
def R_5_6(zeta): return _R_5_generic(zeta, 6)


# ---------------------------------------------------------------------------
# Higher-order approximants (6, 7, 8 pole) computed similarly
# ---------------------------------------------------------------------------

def _compute_npole_coeffs(n_poles, n_prime):
    """Compute n-pole Pade approximant coefficients.

    R_{n,n'}(zeta) = (1 + a1*z + ... + a_{n-2}*z^{n-2}) / (1 + b1*z + ... + b_{2n-3}*z^{2n-3})
    Wait, for n-pole: numerator has n-1 unknowns (a1..a_{n-2}), denominator has n unknowns (b1..b_n)?
    No: for n-pole, the denominator has degree n (n unknowns b1..b_n when b0=1),
    and the numerator degree is ceil(n/2)-1 for the form that gives -1/(2z^2) asymptotic.

    Actually from the paper's pattern:
    - 1-pole: num degree 0, den degree 1 => 1 unknown total
    - 2-pole: R_{2,0} uses Z_{2,0} which is 2-pole in Z
    - 3-pole: num degree 1 (a1), den degree 3 (b1,b2,b3) => 4 unknowns
    - 4-pole: num degree 2 (a1,a2), den degree 4 (b1..b4) => 6 unknowns
    - 5-pole: num degree 3 (a1,a2,a3), den degree 5 (b1..b5) => 8 unknowns
    - 6-pole: num degree 3 (a1,a2,a3), den degree 6 (b1..b6) => 9 unknowns
    Wait, that doesn't follow the pattern. Let me re-read.

    From the paper: R has the form
    R_n(zeta) = (a0 + a1*z + ... ) / (1 + b1*z + ... + b_n*z^n)
    with a0 = 1.

    For the asymptotic -1/(2z^2), the leading term is a_{num_deg}/b_n * z^{num_deg-n} = -1/(2z^2)
    So num_deg - n = -2, i.e., num_deg = n - 2.

    So:
    - 3-pole: num degree 1, den degree 3, unknowns: a1 + b1,b2,b3 = 4
    - 4-pole: num degree 2, den degree 4, unknowns: a1,a2 + b1..b4 = 6
    - 5-pole: num degree 3, den degree 5, unknowns: a1,a2,a3 + b1..b5 = 8
    - 6-pole: num degree 4, den degree 6, unknowns: a1..a4 + b1..b6 = 10
    - 7-pole: num degree 5, den degree 7, unknowns: a1..a5 + b1..b7 = 12
    - 8-pole: num degree 6, den degree 8, unknowns: a1..a6 + b1..b8 = 14
    """
    sp = np.sqrt(np.pi)
    # Power series coefficients of R(zeta)
    # R = 1 + i*sp*z - 2*z^2 - i*sp*z^3 + 4/3*z^4 + i*sp/2*z^5 - 8/15*z^6
    #     - i*sp/6*z^7 + 16/105*z^8 + i*sp/24*z^9 - 32/945*z^10 ...
    # General: c_k for R(z) = sum c_k z^k
    # c_0 = 1, c_1 = i*sp
    # For k >= 2: c_k = (-1)^(k/2) * 2^(k/2) / (k/2)! ... no
    # Actually R(z) = 1 + z*Z(z) where Z has power series
    # Z(z) = i*sp * sum_{n=0}^inf (-1)^n z^(2n) / (2n+1)!! * 2^n ... no
    # Z(z) = i*sp - 2z + i*sp*(-z^2) + 4/3*z^3 + ...
    # Let me just compute enough terms numerically
    c = np.zeros(30, dtype=complex)
    c[0] = 1.0
    c[1] = 1j * sp
    # R(z) = 1 + z*Z(z), Z(z) = i*sp*exp(-z^2)*(1+erf(iz))
    # Power series of Z(z): Z(z) = i*sp * sum_{n=0}^inf (2iz)^n / (1*3*5*...*(2n+1))
    # = i*sp * [1 + 2iz/1 + (2iz)^2/(1*3) + (2iz)^3/(1*3*5) + ...]
    # Hmm, let me just use: Z(z) = sum_{k=0}^inf z_k * z^k
    # z_0 = i*sp, z_1 = -2, z_2 = -i*sp, z_3 = 4/3, z_4 = i*sp/2, z_5 = -8/15, ...
    # General: z_k = (-1)^k * 2 * (i*sp) / (k! * ...) ... use recurrence
    # Actually: Z'(z) = -2(1 + z*Z(z)) = -2*R(z)
    # So z_{k+1} * (k+1) = -2 * c_k where c_k = R coefficient
    # And c_k = z_{k-1} for k >= 1 (since R = 1 + z*Z)
    # Wait: R(z) = 1 + z*Z(z) means c_0 = 1, and for k >= 1: c_k = z_{k-1}
    # And Z'(z) = -2*R(z) means (k+1)*z_{k+1} = -2*c_k

    z_coeffs = np.zeros(30, dtype=complex)
    z_coeffs[0] = 1j * sp
    for k in range(29):
        c_k = z_coeffs[k] if k > 0 else 1.0  # c_0 = 1, c_k = z_{k-1} for k>=1
        if k == 0:
            c_k = 1.0
        else:
            c_k = z_coeffs[k - 1]
        z_coeffs[k + 1] = -2.0 * c_k / (k + 1)

    # Now compute c: c_0=1, c_k = z_{k-1} for k >= 1
    c[0] = 1.0
    for k in range(1, 30):
        c[k] = z_coeffs[k - 1]

    num_deg = n_poles - 2  # degree of numerator
    den_deg = n_poles      # degree of denominator

    n_a = num_deg  # number of a coefficients (a1..a_{num_deg})
    n_b = den_deg  # number of b coefficients (b1..b_{den_deg})
    n_unknowns = n_a + n_b

    # Power series equations to match: n_ps = n_unknowns - (1 + n_prime)
    n_as = 1 + n_prime  # asymptotic constraints
    n_ps = n_unknowns - n_as  # power series constraints

    if n_ps < 0:
        return None

    # Build linear system Ax = rhs
    A_mat = np.zeros((n_unknowns, n_unknowns), dtype=complex)
    rhs = np.zeros(n_unknowns, dtype=complex)

    eq_idx = 0

    # Power series: match coefficient of z^k for k=1..n_ps
    # num(z) = R(z) * den(z)
    # coeff of z^k: a_k = sum_{j=0}^{min(k,den_deg)} c_{k-j} * b_j  (a_0=1, b_0=1)
    # a_k - sum_{j=1}^{min(k,den_deg)} c_{k-j} * b_j = c_k  (for k <= num_deg)
    # -sum_{j=1}^{min(k,den_deg)} c_{k-j} * b_j = c_k          (for k > num_deg)

    for k in range(1, n_ps + 1):
        if eq_idx >= n_unknowns:
            break
        row = np.zeros(n_unknowns, dtype=complex)

        # a_k term (if k <= num_deg)
        if 1 <= k <= num_deg:
            row[k - 1] = 1.0

        # -b_j * c_{k-j} terms
        for j in range(1, min(k, den_deg) + 1):
            c_kj = c[k - j] if (k - j) < len(c) else 0
            row[n_a + j - 1] = -c_kj

        A_mat[eq_idx] = row
        rhs[eq_idx] = c[k] if k < len(c) else 0
        eq_idx += 1

    # Asymptotic constraints
    # These come from matching the asymptotic expansion of R for large z.
    # The general pattern:
    # b_{den_deg} = -2 * a_{num_deg}
    # b_{den_deg-1} = -2 * a_{num_deg-1}
    # b_{den_deg-2} = -2*a_{num_deg-2} + 3*a_{num_deg}   (from 1/z^4 matching)
    # b_{den_deg-3} = -2*a_{num_deg-3} + 3*a_{num_deg-1}  (from 1/z^5)
    # ...and so on

    # Actually, the simplest approach: use the general asymptotic expansion
    # and match term by term. The expansion of N(z)/D(z) for large z
    # is a_{num_deg}/b_{den_deg} * 1/z^2 + ...

    # For the first two asymptotic constraints, the pattern is always:
    # b_n = -2*a_{n-2}  and  b_{n-1} = -2*a_{n-3}
    # For higher ones, cross terms appear.

    # Let's use a more robust approach: match the asymptotic series directly
    # R(z) ~ -1/(2z^2) - 3/(4z^4) - 15/(8z^6) - 105/(16z^8) - ...
    # = sum_{m=1}^inf -(2m-1)!!/(2^m * z^{2m})
    asym_coeffs = []  # coefficients of 1/z^{2}, 1/z^{3}, 1/z^{4}, ...
    # R ~ sum_{j=2}^inf r_j / z^j
    # r_2 = -1/2, r_3 = 0, r_4 = -3/4, r_5 = 0, r_6 = -15/8, ...
    r_asym = np.zeros(30, dtype=complex)
    r_asym[2] = -0.5
    r_asym[4] = -0.75
    r_asym[6] = -15.0/8
    r_asym[8] = -105.0/16
    r_asym[10] = -945.0/32
    # Odd terms are 0

    # N(z)/D(z) = (a_{nd}*z^{nd} + ... + 1) / (b_n*z^n + ... + 1)
    # For large z: = a_{nd}/(b_n * z^{n-nd}) * [1 + lower/z + ...] / [1 + lower/z + ...]
    # = a_{nd}/b_n * z^{-(n-nd)} * (1 + ...)

    # n - nd = den_deg - num_deg = 2, so leading is a_{nd}/b_n / z^2

    # To match asymptotic to the required number of terms, we need to expand
    # N(z)/D(z) as a Laurent series for large z and match each coefficient.

    # N(1/w)/D(1/w) where w = 1/z, expand around w=0
    # N(1/w) = 1 + a1/w + a2/w^2 + a3/w^3 + ...
    # D(1/w) = 1 + b1/w + b2/w^2 + ... + b_n/w^n
    # Multiply both by w^n:
    # w^n * N(1/w) = w^n + a1*w^{n-1} + ... + a_{nd}*w^{n-nd}
    # w^n * D(1/w) = w^n + b1*w^{n-1} + ... + b_n
    # R(1/w) = [w^n*N] / [w^n*D] where we just need the ratio

    # Actually N(z)*D(z)^{-1} for large z means N/D = sum r_j z^{-j}
    # => N(z) = D(z) * sum r_j z^{-j}
    # Matching coefficient of z^{n-2-m} for m = 0, 1, 2, ..., n_as-1:
    # a_{nd-m} or 0 = sum over appropriate products

    # This is getting complicated. Let me just use the direct matching approach
    # from the power series side, which is already set up above, and for the
    # asymptotic constraints, use the relationship that comes from
    # N(z) = R_asym(z) * D(z) matched for high powers of z.

    # For z^{nd}: a_{nd} = r_2 * b_n  (from z^{nd} = z^{n-2} so z^{n-2})
    # r_2 = -1/2, so a_{nd} = -b_n/2 => b_n = -2*a_{nd}

    # For z^{nd-1}: a_{nd-1} = r_2*b_{n-1} + r_3*b_n
    # r_3 = 0, so a_{nd-1} = -b_{n-1}/2 => b_{n-1} = -2*a_{nd-1}

    # For z^{nd-2}: a_{nd-2} = r_2*b_{n-2} + r_3*b_{n-1} + r_4*b_n
    # = -b_{n-2}/2 + 0 - 3/4*b_n
    # So b_{n-2} = -2*a_{nd-2} - 3/2*b_n = -2*a_{nd-2} + 3*a_{nd}

    # For z^{nd-3}: a_{nd-3} or 0 = r_2*b_{n-3} + r_3*b_{n-2} + r_4*b_{n-1} + r_5*b_n
    # = -b_{n-3}/2 + 0 - 3/4*b_{n-1} + 0
    # So b_{n-3} = -2*a_{nd-3} - 3/2*b_{n-1} = -2*a_{nd-3} + 3*a_{nd-1}

    # For z^{nd-4}: a_{nd-4} or 0 = r_2*b_{n-4} + r_4*b_{n-2} + r_6*b_n
    # = -b_{n-4}/2 - 3/4*b_{n-2} - 15/8*b_n
    # So b_{n-4} = -2*a_{nd-4} - 3/2*b_{n-2} - 15/4*b_n

    # General pattern for m-th asymptotic constraint:
    # Matching z^{nd-m}: coefficient on left is a_{nd-m} if nd-m >= 1, else 0 if nd-m=0 (which is 1), else doesn't exist
    # On right: sum_{j even} r_{j} * b_{n-m+2-j}  ... no

    # Let me just implement this properly.
    nd = num_deg
    n = den_deg

    for m in range(n_as):
        if eq_idx >= n_unknowns:
            break
        # Match coefficient of z^{nd-m} in N(z) = R_asym(z) * D(z)
        # Left side: a_{nd-m} if 1 <= nd-m <= nd, 1.0 if nd-m == 0, else 0
        # Right side: sum_{j=2}^{...} r_j * b_{n-j-(nd-m-n)} ... this is confusing

        # Better: N(z) = R(z)*D(z). For large z, powers run from z^n (highest in D)
        # The coefficient of z^{nd-m}:
        # Left: lhs = a_{nd-m} if 1<=nd-m<=nd, else (1 if nd-m==0, else 0)
        # Right: sum over p,q where D has z^p term (b_p, p=0..n) and R_asym has z^{-j} term (r_j)
        # such that p - j = nd - m, i.e., j = p - nd + m
        # For p=0..n: j = p - nd + m, and r_j exists for j >= 2
        # So: p >= nd - m + 2

        row = np.zeros(n_unknowns, dtype=complex)
        rhs_val = 0.0 + 0j

        # Left side
        idx_a = nd - m
        if idx_a == 0:
            rhs_val = 1.0  # a_0 = 1
        elif 1 <= idx_a <= nd:
            row[idx_a - 1] = 1.0  # a_{idx_a}

        # Right side: - sum_p r_{p-nd+m} * b_p for relevant p
        for p in range(n + 1):
            j = p - nd + m
            if j < 2 or j >= len(r_asym):
                continue
            r_j = r_asym[j]
            if abs(r_j) < 1e-30:
                continue
            if p == 0:
                rhs_val += r_j  # b_0 = 1
            else:
                row[n_a + p - 1] -= r_j  # -r_j * b_p contribution to row

        A_mat[eq_idx] = row
        rhs[eq_idx] = rhs_val
        eq_idx += 1

    if eq_idx < n_unknowns:
        return None

    try:
        x = np.linalg.solve(A_mat, rhs)
    except np.linalg.LinAlgError:
        return None

    a_coeffs = np.zeros(num_deg + 1, dtype=complex)
    a_coeffs[0] = 1.0
    for i in range(num_deg):
        a_coeffs[i + 1] = x[i]

    b_coeffs = np.zeros(den_deg + 1, dtype=complex)
    b_coeffs[0] = 1.0
    for i in range(den_deg):
        b_coeffs[i + 1] = x[n_a + i]

    return a_coeffs, b_coeffs


# Cache for generic n-pole approximants
_npole_cache = {}

def R_npole(zeta, n_poles, n_prime):
    """Generic n-pole Pade approximant R_{n,n'}(zeta)."""
    key = (n_poles, n_prime)
    if key not in _npole_cache:
        result = _compute_npole_coeffs(n_poles, n_prime)
        _npole_cache[key] = result

    coeffs = _npole_cache[key]
    if coeffs is None:
        return np.full_like(np.asarray(zeta, dtype=complex), np.nan)

    a_coeffs, b_coeffs = coeffs
    z = np.asarray(zeta, dtype=complex)
    num = np.zeros_like(z)
    den = np.zeros_like(z)
    for i, a in enumerate(a_coeffs):
        num = num + a * z**i
    for i, b in enumerate(b_coeffs):
        den = den + b * z**i
    return num / den


# Convenience wrappers for commonly used approximants
def R_6_4(zeta): return R_npole(zeta, 6, 4)
def R_7_5(zeta): return R_npole(zeta, 7, 5)
def R_8_3(zeta): return R_npole(zeta, 8, 3)


# ---------------------------------------------------------------------------
# Dispersion relation solvers
# ---------------------------------------------------------------------------

def solve_ion_acoustic_kinetic(tau_values, n_pts=200):
    """Solve the exact kinetic dispersion relation for the ion-acoustic mode.

    The dispersion relation (Eq. 423):
    tau * R(zeta_p) + R(zeta_e) = 0

    where tau = T_{0,e}/T_{0,p}, zeta_p = omega/(|k_par|*v_th_p),
    zeta_e = omega/(|k_par|*v_th_e) = zeta_p * sqrt(m_e/m_p) / sqrt(tau)

    We solve for zeta_p as a function of tau.
    In the long-wavelength limit, the ion-acoustic mode approaches a constant zeta.
    """
    from scipy.optimize import fsolve

    results = []
    mp_me = 1836.15  # proton-to-electron mass ratio

    for tau in tau_values:
        def equations(x):
            zr, zi = x
            zeta_p = zr + 1j * zi
            zeta_e = zeta_p / np.sqrt(tau * mp_me)
            val = tau * R_exact(zeta_p) + R_exact(zeta_e)
            return [np.real(val), np.imag(val)]

        # Initial guess from asymptotic: zeta_p ~ sqrt(tau/2) for large tau
        # For tau ~ 1: zeta ~ 1.5 - 0.6i
        if tau <= 1:
            z0 = [1.4, -0.6]
        elif tau <= 5:
            z0 = [1.7 + 0.15 * (tau - 2), -0.4 + 0.05 * (tau - 2)]
        elif tau <= 20:
            z0 = [np.sqrt(tau) * 0.9, -0.1]
        else:
            z0 = [np.sqrt(tau / 2) * 1.02, -0.07]

        try:
            sol = fsolve(equations, z0, full_output=True)
            x_sol = sol[0]
            zeta_p = x_sol[0] + 1j * x_sol[1]
            results.append(zeta_p)
        except Exception:
            results.append(np.nan + 1j * np.nan)

    return np.array(results)


def solve_ion_acoustic_fluid(tau_values, R_func):
    """Solve the fluid dispersion relation for the ion-acoustic mode.

    Replace R(zeta) with R_{n,n'}(zeta) in:
    tau * R_{n,n'}(zeta_p) + R_{n,n'}(zeta_e) = 0
    """
    from scipy.optimize import fsolve

    results = []
    mp_me = 1836.15

    for tau in tau_values:
        def equations(x):
            zr, zi = x
            zeta_p = zr + 1j * zi
            zeta_e = zeta_p / np.sqrt(tau * mp_me)
            val = tau * R_func(zeta_p) + R_func(zeta_e)
            return [np.real(val), np.imag(val)]

        if tau <= 1:
            z0 = [1.4, -0.5]
        elif tau <= 5:
            z0 = [1.7 + 0.15 * (tau - 2), -0.3]
        elif tau <= 20:
            z0 = [np.sqrt(tau) * 0.9, -0.1]
        else:
            z0 = [np.sqrt(tau / 2) * 1.02, -0.07]

        try:
            sol = fsolve(equations, z0, full_output=True)
            x_sol = sol[0]
            zeta_p = x_sol[0] + 1j * x_sol[1]
            results.append(zeta_p)
        except Exception:
            results.append(np.nan + 1j * np.nan)

    return np.array(results)


def solve_langmuir_kinetic(k_lambda_D_values):
    """Solve the exact kinetic dispersion relation for the Langmuir mode.

    Eq. (453): 1 + 1/(k^2 * lambda_D^2) * R(zeta_e) = 0
    where zeta_e = omega / (|k|*v_th_e) and k*lambda_D is the input.

    Returns complex zeta_e.
    """
    from scipy.optimize import fsolve

    results = []

    for klD in k_lambda_D_values:
        def equations(x):
            zr, zi = x
            zeta = zr + 1j * zi
            val = 1.0 + R_exact(zeta) / klD**2
            return [np.real(val), np.imag(val)]

        # Initial guess from fluid: omega^2 = omega_pe^2 + 3*k^2*v_th^2
        # zeta = omega/(k*v_th), omega^2/omega_pe^2 = 1 + 3*k^2*lD^2
        # zeta^2 = omega^2/(k^2*v_th^2) = (omega_pe/(k*v_th))^2 * (1+3*k^2*lD^2)
        # = (1/klD)^2 * (1 + 3*klD^2) / 2  (since v_th^2 = 2*T/m)
        # Actually lambda_D^2 = T/(m*omega_pe^2), v_th^2 = 2T/m, so lD^2 = v_th^2/(2*omega_pe^2)
        # So omega_pe^2/(k^2*v_th^2) = 1/(2*klD^2)
        # zeta^2 = (1+3*klD^2)/(2*klD^2)
        zeta_r_guess = np.sqrt((1 + 3 * klD**2) / (2 * klD**2))
        # Damping estimate from asymptotic formula
        zeta_i_guess = -np.sqrt(np.pi / 8) / klD**3 * np.exp(-(1 + 3 * klD**2) / (2 * klD**2))
        zeta_i_guess = max(zeta_i_guess, -zeta_r_guess * 0.5)
        zeta_i_guess = min(zeta_i_guess, -1e-6)

        z0 = [zeta_r_guess, zeta_i_guess]

        try:
            sol = fsolve(equations, z0, full_output=True)
            x_sol = sol[0]
            zeta = x_sol[0] + 1j * x_sol[1]
            results.append(zeta)
        except Exception:
            results.append(np.nan + 1j * np.nan)

    return np.array(results)


def solve_langmuir_fluid(k_lambda_D_values, R_func):
    """Solve the Langmuir dispersion relation with a Pade approximant."""
    from scipy.optimize import fsolve

    results = []

    for klD in k_lambda_D_values:
        def equations(x):
            zr, zi = x
            zeta = zr + 1j * zi
            val = 1.0 + R_func(zeta) / klD**2
            return [np.real(val), np.imag(val)]

        zeta_r_guess = np.sqrt((1 + 3 * klD**2) / (2 * klD**2))
        zeta_i_guess = -0.01 * zeta_r_guess
        z0 = [zeta_r_guess, zeta_i_guess]

        try:
            sol = fsolve(equations, z0, full_output=True)
            x_sol = sol[0]
            zeta = x_sol[0] + 1j * x_sol[1]
            results.append(zeta)
        except Exception:
            results.append(np.nan + 1j * np.nan)

    return np.array(results)


# Registry of all approximant functions
PADE_REGISTRY = {
    (1, 0): R_1,
    (2, 0): R_2_0,
    (3, 0): R_3_0,
    (3, 1): R_3_1,
    (3, 2): R_3_2,
    (4, 0): R_4_0,
    (4, 1): R_4_1,
    (4, 2): R_4_2,
    (4, 3): R_4_3,
    (4, 4): R_4_4,
    (5, 0): R_5_0,
    (5, 1): R_5_1,
    (5, 2): R_5_2,
    (5, 3): R_5_3,
    (5, 4): R_5_4,
    (5, 5): R_5_5,
    (5, 6): R_5_6,
    (6, 4): R_6_4,
    (7, 5): R_7_5,
    (8, 3): R_8_3,
}
