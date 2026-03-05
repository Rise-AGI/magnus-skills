from magnus import submit_job, JobType
from typing import Annotated, Literal, Optional

ThermalQuenchTime = Annotated[float, {"label": "Thermal Quench Time (ms)", "description": "Timescale of exponential temperature decay during thermal quench"}]
FinalTemperature = Annotated[float, {"label": "Final Temperature (eV)", "description": "Post-quench equilibrium temperature in eV"}]
SeedFraction = Annotated[float, {"label": "RE Seed Fraction", "description": "Initial runaway electron current as fraction of total current"}]

def blueprint(
    tau_tq_ms: ThermalQuenchTime = 5.0,
    t_final_ev: FinalTemperature = 5.0,
    seed_frac: SeedFraction = 0.001,
):
    script = f'''
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.constants import e, m_e, epsilon_0, mu_0, c, pi, m_p
from scipy.integrate import solve_ivp

class PlasmaParams:
    def __init__(self, R=10.0, a=1.0, B_phi0=1.0, T0_keV=1.7, n0=1e20,
                 Z_eff=1.0, eta0=1.1e-7, lnLambda=15.0, I_p_MA=0.67):
        self.R = R; self.a = a; self.B_phi0 = B_phi0
        self.T0 = T0_keV * 1e3 * e; self.n0 = n0
        self.Z_eff = Z_eff; self.eta0 = eta0; self.lnLambda = lnLambda
        self.I_p = I_p_MA * 1e6
    def resistivity(self, T_J):
        return self.eta0 * max(T_J / self.T0, 1e-10)**(-1.5)
    def dreicer_field(self, T_J, n_e):
        return n_e * e**3 * self.lnLambda / (4 * pi * epsilon_0**2 * T_J)
    def critical_field(self, n_e):
        return n_e * e**3 * self.lnLambda / (4 * pi * epsilon_0**2 * m_e * c**2)
    def collision_freq(self, T_J, n_e):
        v_te = np.sqrt(2 * T_J / m_e)
        return n_e * e**4 * self.lnLambda / (4 * pi * epsilon_0**2 * m_e**2 * v_te**3)

params = PlasmaParams()
T_final = {t_final_ev} * e
tau_TQ = {tau_tq_ms} * 1e-3
cross_section = pi * params.a**2
L = mu_0 * params.R * (np.log(8 * params.R / params.a) - 1.5)
I_r_seed = {seed_frac} * params.I_p

def temperature(t):
    return T_final + (params.T0 - T_final) * np.exp(-t / tau_TQ)

def rhs(t, y):
    I_th, I_r = max(y[0], 0), max(y[1], 0)
    T = temperature(t)
    eta = params.resistivity(T)
    E_par = eta * I_th / cross_section
    R_eff = eta * 2 * pi * params.R / cross_section
    dI_th_dt = -R_eff * I_th / L
    E_c = params.critical_field(params.n0)
    eps_c = abs(E_par) / E_c if E_c > 0 else 0
    gamma_aval = 0
    if eps_c > 1.0 and I_r > 0:
        Z = params.Z_eff
        nu_fp = params.n0 * e**4 * params.lnLambda / (4 * pi * epsilon_0**2 * m_e**2 * c**3)
        growth = (eps_c - 1) / params.lnLambda
        sqrt_fac = np.sqrt(pi / (3 * (Z + 5)))
        denom_arg = max(eps_c**2 + 3, 1.0)
        bracket = max(1 - 1.0/eps_c + 4*pi*(Z+1)**2/(3*(Z+5)*denom_arg), 0.01)
        gamma_aval = nu_fp * growth * sqrt_fac * bracket**(-0.5)
    dI_r_dt = gamma_aval * I_r
    E_D = params.dreicer_field(T, params.n0)
    eps_d = abs(E_par) / E_D if E_D > 0 else 0
    S_dreicer = 0
    if eps_d >= 0.01:
        Z = params.Z_eff
        nu_ee = params.collision_freq(T, params.n0)
        power = -3.0/16.0*(1+Z)
        exp_arg = -0.25/eps_d - np.sqrt(1+Z)/np.sqrt(eps_d)
        if exp_arg > -500:
            correction = T/(m_e*c**2)*(eps_d**(-2)/8.0 + (2.0/3.0)*np.sqrt(1+Z)*eps_d**(-1.5))
            S_dreicer = (0.21+0.11*Z)*params.n0*nu_ee*eps_d**power*np.exp(exp_arg)*abs(correction)
    vol = cross_section * 2 * pi * params.R
    dI_r_dreicer = e * c * S_dreicer * vol * 0.5
    dI_r_dt += dI_r_dreicer
    dI_th_dt -= dI_r_dreicer
    return [dI_th_dt, dI_r_dt]

sol = solve_ivp(rhs, (0, 70e-3), [params.I_p - I_r_seed, I_r_seed],
                method="Radau", max_step=3.5e-4, rtol=1e-8, atol=[1.0, 1e-5], dense_output=True)
t_eval = np.linspace(0, 70e-3, 500)
y = sol.sol(t_eval)
times = t_eval * 1e3
I_th = np.maximum(y[0], 0) * 1e-6
I_r = np.maximum(y[1], 0) * 1e-6
I_total = I_th + I_r
np.savetxt("fig3a_current_quench.csv",
           np.column_stack([times, I_total, I_r, I_th]),
           delimiter=",", header="time_ms,I_total_MA,I_re_MA,I_thermal_MA",
           comments="", fmt="%.8e")
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(times, I_total, "b-", lw=2, label="I_total")
ax.plot(times, I_r, "r--", lw=2, label="I_r (RE)")
ax.plot(times, I_th, "g:", lw=2, label="I_th (thermal)")
ax.set_xlabel("Time [ms]"); ax.set_ylabel("Current [MA]")
ax.set_title("Current quench with RE generation"); ax.legend(); ax.grid(True, alpha=0.3)
fig.savefig("fig3a_current_quench.png", dpi=150); plt.close()
print(f"Final: I_total={{I_total[-1]:.4f}}, I_r={{I_r[-1]:.4f}}, I_th={{I_th[-1]:.4f}} MA")
'''

    submit_job(
        task_name="Bandaru2018: Current Quench with RE Generation",
        description=f"1D current quench simulation (tau_TQ={tau_tq_ms}ms, T_f={t_final_ev}eV, seed={seed_frac})",
        namespace="Rise-AGI",
        repo_name="magnus-skills",
        commit_sha="HEAD",
        entry_command=f'python3 -m pip install --break-system-packages numpy scipy matplotlib && python3 -c "{script}"',
        container_image="docker://git.pku.edu.cn/2200011363/knowledge-distiller:latest",
    )
