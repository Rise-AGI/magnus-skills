from magnus import submit_job, JobType
from typing import Annotated, Literal, Optional

IpMA = Annotated[float, {"label": "Initial Plasma Current (MA)", "description": "Initial total plasma current in mega-amperes"}]
SeedFrac = Annotated[float, {"label": "RE Seed Fraction", "description": "Initial RE current as fraction of total current"}]
SimTime = Annotated[float, {"label": "Simulation Time (ms)", "description": "Total simulation time in milliseconds"}]

def blueprint(
    i_p_ma: IpMA = 14.5,
    seed_frac: SeedFrac = 0.001,
    sim_time_ms: SimTime = 12.0,
):
    script = f'''
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.constants import e, m_e, epsilon_0, mu_0, c, pi, m_p
from scipy.integrate import solve_ivp

I_p0 = {i_p_ma} * 1e6
B_phi = 4.8; n_e = 5e19; eta = 1.24e-4; a = 2.0; R = 6.2
T_eV = 2.35; Z = 1.0; lnL = 15.0
I_r0 = {seed_frac} * I_p0
tau_R = mu_0 * a**2 / eta
E_c = n_e * e**3 * lnL / (4 * pi * epsilon_0**2 * m_e * c**2)
nu_fp = n_e * e**4 * lnL / (4 * pi * epsilon_0**2 * m_e**2 * c**3)

def rhs(t, y):
    I_total, I_r = y
    I_r = max(I_r, 0); I_total = max(I_total, I_r)
    I_th = I_total - I_r
    E_par = eta * I_th / (pi * a**2)
    dI_th_dt = -I_th / tau_R
    eps_c = abs(E_par) / E_c if E_c > 0 else 0
    gamma_aval = 0
    if eps_c > 1.0 and I_r > 0:
        growth = (eps_c - 1) / lnL
        sqrt_fac = np.sqrt(pi / (3 * (Z + 5)))
        denom_arg = max(eps_c**2 + 3, 1.0)
        bracket = max(1 - 1.0/eps_c + 4*pi*(Z+1)**2/(3*(Z+5)*denom_arg), 0.01)
        gamma_aval = min(nu_fp * growth * sqrt_fac * bracket**(-0.5), 500.0)
    dI_r_dt = gamma_aval * I_r
    if I_r + dI_r_dt * 1e-6 > I_total:
        dI_r_dt = max(0, (I_total - I_r) * 100)
    return [dI_th_dt + dI_r_dt, dI_r_dt]

t_final = {sim_time_ms} * 1e-3
sol = solve_ivp(rhs, (0, t_final), [I_p0, I_r0],
                method="RK45", max_step=1e-6, rtol=1e-8, atol=1e-10)
times = sol.t * 1e3
I_total = sol.y[0] * 1e-6; I_r = sol.y[1] * 1e-6

# Baseline (no REs)
sol2 = solve_ivp(lambda t, y: [-y[0]/tau_R], (0, t_final), [I_p0],
                 method="RK45", max_step=1e-6, rtol=1e-8)
t2 = sol2.t * 1e3; I_base = sol2.y[0] * 1e-6

t_common = np.linspace(0, {sim_time_ms}, 500)
I_total_i = np.interp(t_common, times, I_total)
I_r_i = np.interp(t_common, times, I_r)
I_base_i = np.interp(t_common, t2, I_base)

np.savetxt("fig6a_vde_currents.csv",
           np.column_stack([t_common, I_total_i, I_r_i, I_total_i-I_r_i, I_base_i]),
           delimiter=",", header="time_ms,I_total_RE_MA,I_r_MA,I_th_RE_MA,I_total_baseline_MA",
           comments="", fmt="%.8e")
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(t_common, I_total_i, "b-", lw=2, label="I_total (with REs)")
ax.plot(t_common, I_r_i, "r-", lw=2, label="I_r (RE)")
ax.plot(t_common, I_total_i-I_r_i, "g--", lw=2, label="I_th (with REs)")
ax.plot(t_common, I_base_i, "k:", lw=2, label="I_total (no REs)")
ax.set_xlabel("Time [ms]"); ax.set_ylabel("Current [MA]")
ax.set_title("ITER VDE current evolution"); ax.legend(); ax.grid(True, alpha=0.3)
ax.set_xlim(0, {sim_time_ms}); ax.set_ylim(0, {i_p_ma}*1.1)
fig.savefig("fig6a_vde_currents.png", dpi=150); plt.close()
print(f"With REs: I_total={{I_total_i[-1]:.2f}}, I_r={{I_r_i[-1]:.2f}} MA")
print(f"Baseline: I_total={{I_base_i[-1]:.2f}} MA")
'''

    submit_job(
        task_name="Bandaru2018: ITER VDE with Runaway Electrons",
        description=f"VDE simulation (I_p={i_p_ma} MA, seed={seed_frac}, t={sim_time_ms} ms)",
        namespace="Rise-AGI",
        repo_name="magnus-skills",
        commit_sha="HEAD",
        entry_command=f'python3 -m pip install --break-system-packages numpy scipy matplotlib && python3 -c "{script}"',
        container_image="docker://git.pku.edu.cn/2200011363/knowledge-distiller:latest",
    )
