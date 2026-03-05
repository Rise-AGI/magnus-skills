from typing import Annotated, Literal, Optional

ContainerImage = "docker://pytorch/pytorch:2.5.1-cuda12.4-cudnn9-runtime"

Alpha = Annotated[float, {"label": "Beam Density Ratio", "description": "alpha = nb/np, beam-to-plasma density ratio"}]
GammaB = Annotated[float, {"label": "Beam Lorentz Factor", "description": "Relativistic Lorentz factor of electron beam"}]
OmegaB = Annotated[float, {"label": "Magnetic Field", "description": "OmegaB = omega_b/omega_p, normalized cyclotron frequency"}]
MassRatio = Annotated[str, {"label": "Mass Ratio", "description": "Electron-to-proton mass ratio: '1/1836' (physical) or '1/100' (PIC)"}]

def blueprint(
    alpha: Alpha = 0.1,
    gamma_b: GammaB = 20.0,
    omega_b: OmegaB = 1.0,
    mass_ratio: MassRatio = "1/1836",
):
    submit_job(
        task_name="Bret2009: Parallel Growth Rates",
        description=f"Compute growth rates of electrostatic (Two-Stream, Buneman) and electromagnetic modes for flow-aligned wave vectors. alpha={alpha}, gamma_b={gamma_b}, OmegaB={omega_b}, R={mass_ratio}",
        entry_command=f"pip install numpy scipy matplotlib && cd /magnus/workspace/repository/submit_Bret_2009/reproduction && python3 fig1_parallel_growth.py",
        container_image=ContainerImage,
    )
