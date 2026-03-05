from typing import Annotated, Literal, Optional

ContainerImage = "docker://pytorch/pytorch:2.5.1-cuda12.4-cudnn9-runtime"

Alpha = Annotated[float, {"label": "Beam Density Ratio", "description": "alpha = nb/np"}]
GammaB = Annotated[float, {"label": "Beam Lorentz Factor", "description": "Relativistic Lorentz factor of electron beam"}]
OmegaB = Annotated[float, {"label": "Magnetic Field", "description": "Normalized cyclotron frequency"}]

def blueprint(
    alpha: Alpha = 0.1,
    gamma_b: GammaB = 2.0,
    omega_b: OmegaB = 3.0,
):
    submit_job(
        task_name="Bret2009: Instability Hierarchy Map",
        description=f"Compute analytical instability hierarchy in (alpha, gamma_b) parameter space for various magnetic field strengths. Uses Tables 1 & 2 growth rate expressions.",
        entry_command=f"pip install numpy scipy matplotlib && cd /magnus/workspace/repository/submit_Bret_2009/reproduction && python3 fig5_hierarchy.py",
        container_image=ContainerImage,
    )
