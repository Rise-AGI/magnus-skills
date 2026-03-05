from typing import Annotated, Literal, Optional

NumParticles = Annotated[int, {"label": "Number of Particles", "description": "Computational particles per species"}]
NumSteps = Annotated[int, {"label": "Time Steps", "description": "Simulation time steps"}]
Anisotropy = Annotated[float, {"label": "Temperature Anisotropy", "description": "a = vthy^2/vthx^2 - 1"}]

def blueprint(
    num_particles: NumParticles = 50000,
    num_steps: NumSteps = 2000,
    anisotropy: Anisotropy = 15.0,
):
    submit_job(
        task_name="EM-PIC: Weibel Instability",
        description="Simulate the Weibel instability with electromagnetic PIC. "
                    "Bi-Maxwellian electron distribution with temperature anisotropy drives "
                    "magnetic field growth. Compares Bz growth rate with linear theory "
                    "(gamma=0.22*wpe for a=15). "
                    "From Markidis & Lapenta, J. Comput. Phys. 230 (2011) 7037-7052.",
        container_image="docker://pytorch/pytorch:2.5.1-cuda12.4-cudnn9-runtime",
        entry_command=(
            "pip install numpy scipy matplotlib && "
            "git clone https://github.com/Rise-AGI/magnus-skills.git repo && "
            "cd repo && git checkout markidis2011 && "
            "cd submit_Markidis_2011/reproduction && "
            "MPLBACKEND=Agg python3 fig10_11_weibel.py"
        ),
    )
