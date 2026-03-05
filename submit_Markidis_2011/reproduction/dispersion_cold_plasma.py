"""
Numerical dispersion relation for the energy-conserving PIC method (Section 5).

For cold plasma (f_0(v) = delta(v)), the dispersion relation is:
  tan(omega*dt/2) * sin(omega*dt/2) = (omega_pe * dt/2)^2

This is compared with the explicit PIC dispersion relation:
  4*sin^2(omega*dt/2) = (omega_pe * dt)^2  ->  exponential growth for omega_pe*dt > 2

Also plots the warm plasma dispersion (Langmuir waves):
  omega^2 = omega_pe^2 + 3*k^2*vth^2
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import brentq


def ec_pic_dispersion_cold(wpe_dt):
    """Solve tan(w*dt/2)*sin(w*dt/2) = (wpe*dt/2)^2 for w*dt.
    Returns all real roots."""
    target = (wpe_dt / 2.0) ** 2
    roots = []
    # Search for roots of f(x) = tan(x)*sin(x) - target = 0 where x = omega*dt/2
    # Branches in (n*pi, (n+1/2)*pi) for each n
    for n in range(20):
        lo = n * np.pi + 0.01
        hi = (n + 0.5) * np.pi - 0.01
        if lo >= hi:
            continue
        try:
            f_lo = np.tan(lo) * np.sin(lo) - target
            f_hi = np.tan(hi) * np.sin(hi) - target
            if f_lo * f_hi < 0:
                x = brentq(lambda x: np.tan(x) * np.sin(x) - target, lo, hi)
                roots.append(2 * x)  # omega * dt
        except Exception:
            pass
    return roots


def explicit_pic_dispersion_cold(wpe_dt):
    """Solve 4*sin^2(omega*dt/2) = (omega_pe*dt)^2.
    Has real solutions only for omega_pe*dt <= 2."""
    if wpe_dt <= 2.0:
        x = np.arcsin(wpe_dt / 2.0)
        return [2 * x]  # omega * dt
    else:
        return []  # No real solution -> exponential growth


# Scan omega_pe * dt
wpe_dt_range = np.linspace(0.1, 8.0, 200)

ec_omega = []
ex_omega = []

for wpe_dt in wpe_dt_range:
    roots_ec = ec_pic_dispersion_cold(wpe_dt)
    if roots_ec:
        ec_omega.append(roots_ec[0] / wpe_dt)  # omega / omega_pe
    else:
        ec_omega.append(np.nan)

    roots_ex = explicit_pic_dispersion_cold(wpe_dt)
    if roots_ex:
        ex_omega.append(roots_ex[0] / wpe_dt)
    else:
        ex_omega.append(np.nan)

ec_omega = np.array(ec_omega)
ex_omega = np.array(ex_omega)

# Save data
data = np.column_stack([wpe_dt_range, ec_omega, ex_omega])
np.savetxt("../data/dispersion_cold_plasma.csv", data, delimiter=",",
           header="wpe_dt,omega_over_wpe_EC,omega_over_wpe_explicit",
           comments="", fmt="%.8e")

# Plot
fig, ax = plt.subplots(figsize=(8, 6))
ax.plot(wpe_dt_range, ec_omega, 'r-', linewidth=2, label='EC-PIC (always stable)')
ax.plot(wpe_dt_range, ex_omega, 'b--', linewidth=2, label='Explicit PIC')
ax.axvline(x=2.0, color='k', linestyle=':', alpha=0.5, label='$\\omega_{pe}\\Delta t = 2$ (explicit stability limit)')
ax.set_xlabel('$\\omega_{pe} \\Delta t$', fontsize=12)
ax.set_ylabel('$\\omega / \\omega_{pe}$', fontsize=12)
ax.set_title('Cold Plasma Dispersion: EC-PIC vs Explicit PIC', fontsize=14)
ax.legend(fontsize=11)
ax.set_xlim([0, 8])
ax.set_ylim([0, 2])
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("../plots/dispersion_cold_plasma.png", dpi=150)
print("Saved dispersion_cold_plasma.png and dispersion_cold_plasma.csv")
