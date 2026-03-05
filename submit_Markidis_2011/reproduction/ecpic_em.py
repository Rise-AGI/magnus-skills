"""
1D3V Electromagnetic PIC code using Boris push + FDTD.

Used for Weibel instability simulation (Figures 10, 11).

Normalization follows the paper (same as electrostatic appendix):
  - omega_pe = 1, c = 1, m_e = 1, |e| absorbed into Q
  - Q = WP^2 / (QM * N/L) absorbs the 4*pi factor
  - Ampere: dE/dt = c*curl_B - J (no explicit 4*pi since it's in Q)
  - Faraday: dB/dt = -c*curl_E
  - Energy: KE = 0.5*|Q|*sum(v^2), FE = 0.5*sum(E^2)*dx, ME = 0.5*sum(B^2)*dx
"""

import numpy as np


def gather_1d(field, x, dx, NG):
    """CIC interpolation: grid -> particles."""
    xn = x / dx
    g1 = np.floor(xn).astype(int) % NG
    g2 = (g1 + 1) % NG
    w2 = xn - np.floor(xn)
    w1 = 1.0 - w2
    return w1 * field[g1] + w2 * field[g2]


def scatter_1d(val, x, dx, NG):
    """CIC scatter: particles -> grid."""
    xn = x / dx
    g1 = np.floor(xn).astype(int) % NG
    g2 = (g1 + 1) % NG
    w2 = xn - np.floor(xn)
    w1 = 1.0 - w2
    out = np.zeros(NG)
    np.add.at(out, g1, w1 * val)
    np.add.at(out, g2, w2 * val)
    return out


def boris_push(vx, vy, vz, Epx, Epy, Epz, Bpx, Bpy, Bpz, qmdt2, c=1.0):
    """Boris velocity push. qmdt2 = (q/m)*dt/2."""
    # Half E push
    vxm = vx + qmdt2 * Epx
    vym = vy + qmdt2 * Epy
    vzm = vz + qmdt2 * Epz
    # Rotation
    tx = qmdt2 / c * Bpx
    ty = qmdt2 / c * Bpy
    tz = qmdt2 / c * Bpz
    t2 = tx**2 + ty**2 + tz**2
    sx = 2*tx/(1+t2); sy = 2*ty/(1+t2); sz = 2*tz/(1+t2)
    vpx = vxm + vym*tz - vzm*ty
    vpy = vym + vzm*tx - vxm*tz
    vpz = vzm + vxm*ty - vym*tx
    vxp = vxm + vpy*sz - vpz*sy
    vyp = vym + vpz*sx - vpx*sz
    vzp = vzm + vpx*sy - vpy*sx
    # Half E push
    return vxp + qmdt2*Epx, vyp + qmdt2*Epy, vzp + qmdt2*Epz


def run_weibel(L, DT, NT, NG, N, WP, QM_e, vthx, vthy, mass_ratio=1836.0):
    """Run 1D3V EM-PIC Weibel instability simulation.

    Uses same normalization as paper's ES code:
    Q = WP^2 / (QM * N/L), which absorbs 4*pi.
    """
    dx = L / NG
    c = 1.0

    Ne = N
    Ni = N

    np.random.seed(42)

    # Electrons: bi-Maxwellian
    xe = np.random.uniform(0, L, Ne)
    vex = vthx * np.random.randn(Ne)
    vey = vthy * np.random.randn(Ne)
    vez = np.zeros(Ne)

    # Ions: same T, heavier
    xi = np.random.uniform(0, L, Ni)
    vthi_x = vthx / np.sqrt(mass_ratio)
    vthi_y = vthy / np.sqrt(mass_ratio)
    vix = vthi_x * np.random.randn(Ni)
    viy = vthi_y * np.random.randn(Ni)
    viz = np.zeros(Ni)

    # Charges: Q = WP^2 / (QM * N/L)
    # For electrons: QM_e = -1, so Qe = WP^2 / (-1 * Ne/L) = -WP^2*L/Ne
    Qe = WP**2 / (QM_e * Ne / L)  # negative
    # For ions: QM_i = +1/mass_ratio
    QM_i = np.abs(QM_e) / mass_ratio
    Qi = -Qe  # positive, neutralizing

    # Fields
    Ex = np.zeros(NG)
    Ey = np.zeros(NG)
    Ez = np.zeros(NG)
    Bx = np.zeros(NG)
    By = np.zeros(NG)
    Bz = np.zeros(NG)

    energy_hist = []
    momentum_hist = []
    Bz_k1_hist = []

    for it in range(NT):
        # --- Record observables ---
        KE_e = 0.5 * np.abs(Qe) * np.sum(vex**2 + vey**2 + vez**2)
        KE_i = 0.5 * np.abs(Qi) * np.sum(vix**2 + viy**2 + viz**2)
        FE = 0.5 * np.sum(Ex**2 + Ey**2 + Ez**2) * dx
        ME = 0.5 * np.sum(Bx**2 + By**2 + Bz**2) * dx
        energy_hist.append(KE_e + KE_i + FE + ME)

        px = np.abs(Qe)/np.abs(QM_e) * np.sum(vex) + Qi/QM_i * np.sum(vix)
        momentum_hist.append(px)

        Bz_hat = np.fft.fft(Bz)
        Bz_k1_hist.append(np.abs(Bz_hat[1]) * 2.0 / NG)

        # --- Half B update (Faraday): dB/dt = -c*curl_E ---
        dEz_dx = (np.roll(Ez, -1) - Ez) / dx
        dEy_dx = (np.roll(Ey, -1) - Ey) / dx
        # curl_E in 1D: (0, -dEz/dx, dEy/dx)
        # dB/dt = -c*curl_E -> dBy/dt = c*dEz/dx, dBz/dt = -c*dEy/dx
        By = By + 0.5 * c * DT * dEz_dx
        Bz = Bz - 0.5 * c * DT * dEy_dx

        # --- Push electrons ---
        Epx = gather_1d(Ex, xe, dx, NG)
        Epy = gather_1d(Ey, xe, dx, NG)
        Epz = gather_1d(Ez, xe, dx, NG)
        Bpx = gather_1d(Bx, xe, dx, NG)
        Bpy = gather_1d(By, xe, dx, NG)
        Bpz = gather_1d(Bz, xe, dx, NG)
        vex, vey, vez = boris_push(vex, vey, vez, Epx, Epy, Epz,
                                     Bpx, Bpy, Bpz, QM_e*DT/2, c)
        xe = (xe + vex * DT) % L

        # --- Push ions ---
        Epx = gather_1d(Ex, xi, dx, NG)
        Epy = gather_1d(Ey, xi, dx, NG)
        Epz = gather_1d(Ez, xi, dx, NG)
        Bpx = gather_1d(Bx, xi, dx, NG)
        Bpy = gather_1d(By, xi, dx, NG)
        Bpz = gather_1d(Bz, xi, dx, NG)
        vix, viy, viz = boris_push(vix, viy, viz, Epx, Epy, Epz,
                                     Bpx, Bpy, Bpz, QM_i*DT/2, c)
        xi = (xi + vix * DT) % L

        # --- Current density (no explicit 4*pi, absorbed in Q) ---
        Jx = scatter_1d(Qe*vex, xe, dx, NG)/dx + scatter_1d(Qi*vix, xi, dx, NG)/dx
        Jy = scatter_1d(Qe*vey, xe, dx, NG)/dx + scatter_1d(Qi*viy, xi, dx, NG)/dx
        Jz = scatter_1d(Qe*vez, xe, dx, NG)/dx + scatter_1d(Qi*viz, xi, dx, NG)/dx

        # --- E update (Ampere): dE/dt = c*curl_B - J ---
        # curl_B in 1D: (0, -dBz/dx, dBy/dx)
        dBz_dx = (Bz - np.roll(Bz, 1)) / dx
        dBy_dx = (By - np.roll(By, 1)) / dx
        Ex = Ex - Jx * DT  # no curl_B_x component in 1D
        Ey = Ey + c * DT * (-dBz_dx) - Jy * DT
        Ez = Ez + c * DT * dBy_dx - Jz * DT

        # --- Second half B update ---
        dEz_dx = (np.roll(Ez, -1) - Ez) / dx
        dEy_dx = (np.roll(Ey, -1) - Ey) / dx
        By = By + 0.5 * c * DT * dEz_dx
        Bz = Bz - 0.5 * c * DT * dEy_dx

        if it % 50 == 0:
            print(f"  Step {it}/{NT}, Etot={energy_hist[-1]:.4e}, |Bz_k1|={Bz_k1_hist[-1]:.4e}")

    return {
        'energy': np.array(energy_hist),
        'momentum': np.array(momentum_hist),
        'Bz_k1': np.array(Bz_k1_hist),
    }
