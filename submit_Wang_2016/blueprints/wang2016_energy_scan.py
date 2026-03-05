"""
Blueprint: Energy limit parametric scan (Wang et al. 2016, Figs. 9/10/12).

Computes the synchrotron energy limit and balance time as a function of
a swept parameter (electric field, magnetic field, major radius, or
safety factor) for runaway electrons in tokamak fields.
"""

SweepParam = Annotated[str, {
    "placeholder": "electric_field",
    "description": "Parameter to sweep: electric_field, major_radius, or safety_factor",
}]

SweepMin = Annotated[float, {
    "placeholder": "0.05",
    "description": "Minimum value of swept parameter",
}]

SweepMax = Annotated[float, {
    "placeholder": "5.0",
    "description": "Maximum value of swept parameter",
}]

NumPoints = Annotated[int, {
    "placeholder": "40",
    "description": "Number of sweep points",
}]

ElectricField = Annotated[float, {
    "placeholder": "0.2",
    "description": "Loop electric field E_l [V/m] (fixed when not swept)",
}]

MagneticField = Annotated[float, {
    "placeholder": "2.0",
    "description": "Toroidal magnetic field B0 [T]",
}]

MajorRadius = Annotated[float, {
    "placeholder": "1.7",
    "description": "Major radius R0 [m] (fixed when not swept)",
}]

SafetyFactor = Annotated[float, {
    "placeholder": "2.0",
    "description": "Safety factor q (fixed when not swept)",
}]


def blueprint(
    sweep_param: SweepParam = "electric_field",
    sweep_min: SweepMin = 0.05,
    sweep_max: SweepMax = 5.0,
    num_points: NumPoints = 40,
    electric_field: ElectricField = 0.2,
    magnetic_field: MagneticField = 2.0,
    major_radius: MajorRadius = 1.7,
    safety_factor: SafetyFactor = 2.0,
):
    description = f"""## Runaway Electron Energy Limit Scan
Sweeping **{sweep_param}** from {sweep_min} to {sweep_max} ({num_points} points).
Wang, Qin & Liu, Phys. Plasmas 23, 062505 (2016)."""

    entry_command = f"""
set -e
pip install numpy scipy matplotlib --quiet
cd /workspace

cat > energy_scan.py << 'PYEOF'
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import json

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

def find_energy_limit(E_l, B0, R0, q, r_orbit=0.1, t_max=8.0):
    kappa_sq = 1.0/R0**2
    if r_orbit > 0:
        kh = r_orbit / (r_orbit**2 + q**2*R0**2)
        kappa_sq += kh**2
    R_eff = 1.0/np.sqrt(kappa_sq)
    alpha = E_CHARGE*E_l/M0C
    tau_s = TAU_S_COEFF*B0**2
    tau_c = TAU_C_COEFF/R_eff**2

    def rhs(t, y):
        pp, pperp = y
        pperp = max(abs(pperp), 1e-30)
        p2 = pp**2 + pperp**2
        if p2 < 1e-60: return [alpha, 0.0]
        gamma = np.sqrt(1.0+p2)
        D = (tau_s*pperp**2 + tau_c*pp**4)*gamma/p2
        return [alpha - D*pp, -D*pperp]

    sol = solve_ivp(rhs, [0, t_max], [5.0, 1.0],
                    t_eval=np.linspace(0, t_max, 5000),
                    method='RK45', rtol=1e-10, atol=1e-15, max_step=t_max/100)
    gamma = np.sqrt(1 + sol.y[0]**2 + sol.y[1]**2)
    energy = (gamma-1)*M0C2_MEV
    E_max = float(np.max(energy))
    threshold = 0.9*E_max
    idx = np.searchsorted(energy, threshold)
    t_blc = float(sol.t[min(idx, len(sol.t)-1)])
    return E_max, t_blc

sweep_param = "{sweep_param}"
sweep_vals = np.linspace({sweep_min}, {sweep_max}, {num_points})
E_l_base = {electric_field}
B0_base = {magnetic_field}
R0_base = {major_radius}
q_base = {safety_factor}

E_max_list, t_blc_list = [], []
for val in sweep_vals:
    E_l = val if sweep_param == "electric_field" else E_l_base
    R0 = val if sweep_param == "major_radius" else R0_base
    q = val if sweep_param == "safety_factor" else q_base
    E_max, t_blc = find_energy_limit(E_l, B0_base, R0, q)
    E_max_list.append(E_max)
    t_blc_list.append(t_blc)
    print(f"  {{val:.3f}} -> E_max={{E_max:.1f}} MeV, t_blc={{t_blc:.3f}} s")

E_max_arr = np.array(E_max_list)
t_blc_arr = np.array(t_blc_list)

labels = {{"electric_field": "E_l [V/m]", "major_radius": "R0 [m]", "safety_factor": "q"}}
xlabel = labels.get(sweep_param, sweep_param)

fig, ax1 = plt.subplots(figsize=(8, 5))
ax1.plot(sweep_vals, E_max_arr, 'b-', linewidth=2, label='$E_{{max}}$')
ax1.set_xlabel(xlabel, fontsize=12)
ax1.set_ylabel('Maximum Energy [MeV]', fontsize=12, color='blue')
ax1.tick_params(axis='y', labelcolor='blue')
ax2 = ax1.twinx()
ax2.plot(sweep_vals, t_blc_arr, 'r--', linewidth=2, label='$t_{{blc}}$')
ax2.set_ylabel('Balance Time [s]', fontsize=12, color='red')
ax2.tick_params(axis='y', labelcolor='red')
fig.legend(loc='upper left', bbox_to_anchor=(0.15, 0.95))
plt.title(f'Energy Limit vs {{xlabel}}', fontsize=12)
plt.tight_layout()
plt.savefig('energy_scan.png', dpi=150, bbox_inches='tight')
plt.close()

np.savetxt('energy_scan.csv',
           np.column_stack([sweep_vals, E_max_arr, t_blc_arr]),
           delimiter=',', header=f'{{sweep_param}},E_max_MeV,t_balance_s',
           comments='', fmt='%.8e')

result = {{
    'sweep_param': sweep_param,
    'E_max_range': [float(E_max_arr.min()), float(E_max_arr.max())],
    't_blc_range': [float(t_blc_arr.min()), float(t_blc_arr.max())],
    'num_points': {num_points},
}}
print(json.dumps(result, indent=2))
with open('/workspace/result.json', 'w') as f:
    json.dump(result, f, indent=2)
PYEOF

python3 energy_scan.py

cat result.json > "$MAGNUS_RESULT"
"""

    submit_job(
        task_name=f"[Blueprint] Energy Scan: {sweep_param} (Wang 2016)",
        description=description,
        namespace="Rise-AGI",
        repo_name="magnus-skills",
        entry_command=entry_command,
        container_image="docker://pytorch/pytorch:2.5.1-cuda12.4-cudnn9-runtime",
        gpu_count=0,
        gpu_type="cpu",
        job_type=JobType.A2,
    )
