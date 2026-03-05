"""
Blueprint: Analytical energy limit calculator (Wang et al. 2016).

Quick analytical estimate of the runaway electron energy limit using
the curvature-radiation balance formula. No ODE integration needed.
Returns gamma_max and E_max for given tokamak parameters.
"""

ElectricField = Annotated[float, {
    "placeholder": "0.2",
    "description": "Loop electric field E_l [V/m]",
}]

MajorRadius = Annotated[float, {
    "placeholder": "1.7",
    "description": "Major radius R0 [m]",
}]

SafetyFactor = Annotated[float, {
    "placeholder": "2.0",
    "description": "Safety factor q",
}]

OrbitRadius = Annotated[float, {
    "placeholder": "0.1",
    "description": "Orbit minor radius r [m]",
}]


def blueprint(
    electric_field: ElectricField = 0.2,
    major_radius: MajorRadius = 1.7,
    safety_factor: SafetyFactor = 2.0,
    orbit_radius: OrbitRadius = 0.1,
):
    description = """## Analytical Runaway Electron Energy Limit
Computes the energy limit from curvature-radiation balance:
gamma^4 = 6*pi*eps0 * E_l * R_eff^2 / e
Wang, Qin & Liu, Phys. Plasmas 23, 062505 (2016)."""

    entry_command = f"""
set -e
python3 -c "
import math, json

E_CHARGE = 1.602176634e-19
EPSILON_0 = 8.8541878128e-12
M0C2_MEV = 0.51099895
PI = math.pi

E_l = {electric_field}
R0 = {major_radius}
q = {safety_factor}
r = {orbit_radius}

kappa_sq = 1.0/R0**2
if r > 0:
    kh = r / (r**2 + q**2*R0**2)
    kappa_sq += kh**2
R_eff = 1.0 / math.sqrt(kappa_sq)

gamma4 = 6*PI*EPSILON_0*E_l*R_eff**2 / E_CHARGE
gamma_max = gamma4**0.25
E_max = (gamma_max - 1)*M0C2_MEV

result = {{
    'R_eff_m': round(R_eff, 4),
    'gamma_max': round(gamma_max, 2),
    'E_max_MeV': round(E_max, 2),
    'parameters': {{'E_l': E_l, 'R0': R0, 'q': q, 'r': r}},
}}
print(json.dumps(result, indent=2))
with open('/tmp/result.json', 'w') as f:
    json.dump(result, f, indent=2)
"

cat /tmp/result.json > "$MAGNUS_RESULT"
"""

    submit_job(
        task_name="[Blueprint] Analytical Energy Limit (Wang 2016)",
        description=description,
        namespace="Rise-AGI",
        repo_name="magnus-skills",
        entry_command=entry_command,
        container_image="docker://pytorch/pytorch:2.5.1-cuda12.4-cudnn9-runtime",
        gpu_count=0,
        gpu_type="cpu",
        job_type=JobType.A2,
    )
