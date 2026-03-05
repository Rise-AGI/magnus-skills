from typing import Annotated, Literal, Optional

EpsIReal = Annotated[float, {"label": "Re[eps_i]", "description": "Real part of metal dielectric constant"}]
EpsIImag = Annotated[float, {"label": "Im[eps_i]", "description": "Imaginary part (loss) of metal dielectric constant"}]
OmegaRc = Annotated[float, {"label": "omega*R/c", "description": "Normalized radius parameter (typically 100-1000)"}]
ThetaDeg = Annotated[float, {"label": "Bend angle (deg)", "description": "Bend angle in degrees"}]

def blueprint(
    eps_i_real: EpsIReal = -15.0,
    eps_i_imag: EpsIImag = 0.5,
    omega_R_over_c: OmegaRc = 800.0,
    theta_deg: ThetaDeg = 90.0,
):
    submit_job(
        task_name="Hasegawa 2004: SPP Scattering Coefficients",
        repo_name="magnus-skills",
        branch="hasegawa2004",
        commit_sha="45d86ae3fc57664f302b5697db69f10727d954c6",
        description=(
            "Compute T, R, radiation loss P, upper bound Tu, coupling efficiency, "
            "and field profiles for SPP propagation around a metallic bend. "
            "From Hasegawa et al., Appl. Phys. Lett. 84, 1835 (2004)."
        ),
        entry_command=(
            "python3 -m pip install --break-system-packages numpy scipy matplotlib mpmath && "
            "cd reproduction && "
            "python3 fig1_coefficients.py"
        ),
        container_image="docker://git.pku.edu.cn/2200011363/knowledge-distiller:latest",
    )
