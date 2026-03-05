from typing import Annotated, Literal, Optional
from magnus import submit_job

HFrac = Annotated[float, {"label": "Slab Height h/P", "description": "Slab height as fraction of period (default: 0.5)"}]
EpsBMax = Annotated[float, {"label": "Max Boundary Eps", "description": "Maximum boundary dielectric constant to scan (default: 12.0)"}]

def blueprint(
    h_frac: HFrac = 0.5,
    eps_b_max: EpsBMax = 12.0,
):
    submit_job(
        task_name="Huttunen2003: Guided Mode Analysis",
        description="Analyze slab waveguide guided modes as function of boundary material. "
                    "Shows how boundary material controls mode propagation.",
        repo_name="magnus-skills",
        branch="huttunen2003",
        entry_command=(
            f"pip install numpy scipy matplotlib && "
            f"cd submit_Huttunen_2003/reproduction && "
            f"MPLBACKEND=Agg python3 fig3_guided_modes.py"
        ),
        container_image="docker://pytorch/pytorch:2.5.1-cuda12.4-cudnn9-runtime",
    )
