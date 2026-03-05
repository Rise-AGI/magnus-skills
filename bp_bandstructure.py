from typing import Annotated, Literal, Optional
from magnus import submit_job

EpsRod = Annotated[float, {"label": "Rod Dielectric Constant", "description": "Dielectric constant of the high-index layer (default: 13.0 for GaAs)"}]
LFrac = Annotated[float, {"label": "Layer Fraction l/P", "description": "Fraction of period occupied by the high-index layer (default: 0.2)"}]
OmegaMax = Annotated[float, {"label": "Max Frequency", "description": "Maximum normalized frequency omega*P/(2*pi*c) (default: 0.6)"}]
NPoints = Annotated[int, {"label": "Frequency Points", "description": "Number of frequency sampling points (default: 3000)"}]

def blueprint(
    eps_rod: EpsRod = 13.0,
    l_frac: LFrac = 0.2,
    omega_max: OmegaMax = 0.6,
    n_points: NPoints = 3000,
):
    submit_job(
        task_name="Huttunen2003: 1D PhC Band Structure",
        description="Compute band structure of 1D photonic crystal using transfer matrix method. "
                    "Produces band diagram and identifies photonic band gaps.",
        repo_name="magnus-skills",
        branch="huttunen2003",
        entry_command=(
            f"pip install numpy scipy matplotlib && "
            f"cd submit_Huttunen_2003/reproduction && "
            f"MPLBACKEND=Agg python3 fig1_bandstructure.py"
        ),
        container_image="docker://git.pku.edu.cn/2200011363/knowledge-distiller:latest",
    )
