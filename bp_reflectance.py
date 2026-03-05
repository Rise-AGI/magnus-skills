from typing import Annotated, Literal, Optional
from magnus import submit_job

NPeriods = Annotated[str, {"label": "Number of Periods", "description": "Comma-separated list of period counts (default: 3,5,10,20)"}]
OmegaMax = Annotated[float, {"label": "Max Frequency", "description": "Maximum normalized frequency (default: 0.6)"}]

def blueprint(
    n_periods: NPeriods = "3,5,10,20",
    omega_max: OmegaMax = 0.6,
):
    submit_job(
        task_name="Huttunen2003: PhC Reflectance Spectrum",
        description="Compute reflectance spectrum of finite 1D photonic crystal stacks "
                    "using transfer matrix method. Shows band gap formation.",
        repo_name="magnus-skills",
        branch="huttunen2003",
        entry_command=(
            f"pip install numpy scipy matplotlib && "
            f"cd submit_Huttunen_2003/reproduction && "
            f"MPLBACKEND=Agg python3 fig2_reflectance_spectrum.py"
        ),
        container_image="docker://git.pku.edu.cn/2200011363/knowledge-distiller:latest",
    )
