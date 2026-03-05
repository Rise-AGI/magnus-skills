from typing import Annotated, Literal, Optional

WavelengthA = Annotated[str, {"label": "Wavelengths (nm)", "description": "Comma-separated wavelengths, e.g. 500,600,700"}]
NpointsR = Annotated[int, {"label": "R grid points", "description": "Number of bend radius points (50-300)"}]
ThetaDeg = Annotated[float, {"label": "Bend angle (degrees)", "description": "Bend angle in degrees"}]

def blueprint(
    wavelengths: WavelengthA = "500,600,700",
    n_points_R: NpointsR = 120,
    theta_deg: ThetaDeg = 90.0,
):
    submit_job(
        task_name="Hasegawa 2004: Upper Bound Transmittance Tu vs R",
        repo_name="magnus-skills",
        branch="hasegawa2004",
        commit_sha="45d86ae3fc57664f302b5697db69f10727d954c6",
        description=(
            "Compute SPP upper bound transmittance Tu as a function of bend radius R "
            "for silver-air interface. Reproduces Figure 2 from Hasegawa et al., "
            "Appl. Phys. Lett. 84, 1835 (2004)."
        ),
        entry_command=(
            "python3 -m pip install --break-system-packages numpy scipy matplotlib mpmath && "
            "cd reproduction && "
            "python3 fig2_transmittance_bound.py"
        ),
        container_image="docker://git.pku.edu.cn/2200011363/knowledge-distiller:latest",
    )
