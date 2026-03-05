"""
Blueprint: Runaway electron momentum evolution (Wang et al. 2016, Fig. 4).

Solves the relaxation (gyro-center) equations for momentum dynamics of
runaway electrons in tokamak fields with synchrotron + curvature radiation.
Outputs momentum trajectory and energy evolution.
"""

ElectricField = Annotated[float, {
    "placeholder": "0.2",
    "description": "Loop electric field E_l [V/m]",
}]

MagneticField = Annotated[float, {
    "placeholder": "2.0",
    "description": "Toroidal magnetic field B0 [T]",
}]

MajorRadius = Annotated[float, {
    "placeholder": "1.7",
    "description": "Major radius R0 [m]",
}]

SafetyFactor = Annotated[float, {
    "placeholder": "2.0",
    "description": "Safety factor q",
}]

SimTime = Annotated[float, {
    "placeholder": "3.5",
    "description": "Simulation time [s]",
}]


def blueprint(
    electric_field: ElectricField = 0.2,
    magnetic_field: MagneticField = 2.0,
    major_radius: MajorRadius = 1.7,
    safety_factor: SafetyFactor = 2.0,
    sim_time: SimTime = 3.5,
):
    description = """## Runaway Electron Momentum Evolution
Solves relaxation equations for runaway electron dynamics in tokamak fields.
Wang, Qin & Liu, Phys. Plasmas 23, 062505 (2016), Figure 4."""

    entry_command = f"""
set -e
pip install numpy scipy matplotlib --quiet
cd /workspace

cat > runaway_momentum.py << 'PYEOF'
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import json, sys

E_CHARGE = 1.602176634e-19
M_E = 9.1093837015e-31
C_LIGHT = 2.99792458e8
EPSILON_0 = 8.8541878128e-12
PI = np.pi
M0C = M_E * C_LIGHT
M0C2_MEV = 0.51099895
C_SYNC = E_CHARGE**4 / (6*PI*EPSILON_0*M_E**4*C_LIGHT**3)
TAU_S_COEFF = C_SYNC * M_E
TAU_C_COEFF = E_CHARGE**2 / (6*PI*EPSILON_0*C_LIGHT*M_E)

E_l = {electric_field}
B0 = {magnetic_field}
R0 = {major_radius}
q = {safety_factor}
t_max = {sim_time}
r_orbit = 0.1
p_par0, p_perp0 = 5.0, 1.0

kappa_sq = 1.0/R0**2
if r_orbit > 0:
    kh = r_orbit / (r_orbit**2 + q**2*R0**2)
    kappa_sq += kh**2
R_eff = 1.0 / np.sqrt(kappa_sq)

alpha = E_CHARGE * E_l / M0C
tau_s = TAU_S_COEFF * B0**2
tau_c = TAU_C_COEFF / R_eff**2

def rhs(t, y):
    pp, pperp = y
    pperp = max(abs(pperp), 1e-30)
    p2 = pp**2 + pperp**2
    if p2 < 1e-60:
        return [alpha, 0.0]
    gamma = np.sqrt(1.0 + p2)
    D = (tau_s*pperp**2 + tau_c*pp**4)*gamma/p2
    return [alpha - D*pp, -D*pperp]

sol = solve_ivp(rhs, [0, t_max], [p_par0, p_perp0],
                t_eval=np.linspace(0, t_max, 10000),
                method='RK45', rtol=1e-10, atol=1e-15, max_step=t_max/100)

p_par = sol.y[0]
p_perp = np.abs(sol.y[1])
gamma = np.sqrt(1 + p_par**2 + p_perp**2)
energy = (gamma - 1) * M0C2_MEV

gamma4_an = 6*PI*EPSILON_0*E_l*R_eff**2/E_CHARGE
gamma_an = gamma4_an**0.25
E_an = (gamma_an - 1)*M0C2_MEV

fig, axes = plt.subplots(1, 2, figsize=(12, 5))
ax = axes[0]
ax.plot(p_par, p_perp, 'b-', linewidth=0.5)
ax.set_xlabel(r'$p_\\parallel / m_0 c$', fontsize=12)
ax.set_ylabel(r'$p_\\perp / m_0 c$', fontsize=12)
ax.set_title('(a) Momentum trajectory', fontsize=12)

ax = axes[1]
ax.plot(sol.t, energy, 'r-', linewidth=1)
ax.set_xlabel('Time [s]', fontsize=12)
ax.set_ylabel('Kinetic Energy [MeV]', fontsize=12)
ax.set_title('(b) Energy evolution', fontsize=12)
ax.axhline(y=energy[-1], color='k', linestyle='--', alpha=0.5,
           label=f'$E_{{max}}$ = {{energy[-1]:.1f}} MeV')
ax.legend(fontsize=10)
plt.tight_layout()
plt.savefig('momentum_evolution.png', dpi=150, bbox_inches='tight')
plt.close()

np.savetxt('momentum_data.csv',
           np.column_stack([sol.t, p_par, p_perp, gamma, energy]),
           delimiter=',',
           header='time,p_parallel,p_perpendicular,gamma,energy_MeV',
           comments='', fmt='%.8e')

result = {{
    'E_max_numerical': float(np.max(energy)),
    'E_max_analytical': float(E_an),
    'gamma_max': float(np.max(gamma)),
    'max_p_parallel': float(np.max(p_par)),
    'parameters': {{'E_l': E_l, 'B0': B0, 'R0': R0, 'q': q}}
}}
print(json.dumps(result, indent=2))
with open('/workspace/result.json', 'w') as f:
    json.dump(result, f, indent=2)
PYEOF

python3 runaway_momentum.py

cat result.json > "$MAGNUS_RESULT"
"""

    submit_job(
        task_name="[Blueprint] Runaway Electron Momentum (Wang 2016)",
        description=description,
        namespace="Rise-AGI",
        repo_name="magnus-skills",
        entry_command=entry_command,
        container_image="docker://pytorch/pytorch:2.5.1-cuda12.4-cudnn9-runtime",
        gpu_count=0,
        gpu_type="cpu",
        job_type=JobType.A2,
    )
