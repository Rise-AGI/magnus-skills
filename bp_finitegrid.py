from typing import Annotated, Literal, Optional

NumParticles = Annotated[int, {"label": "Number of Particles", "description": "Computational particles (higher = less noise)"}]
NumSteps = Annotated[int, {"label": "Time Steps", "description": "Number of simulation cycles"}]

def blueprint(
    num_particles: NumParticles = 10000,
    num_steps: NumSteps = 200,
):
    submit_job(
        task_name="EC-PIC: Finite Grid Instability Test",
        description="Test energy-conserving PIC against finite grid instability. "
                    "Maxwellian plasma with dx >> lambda_D: EC-PIC remains stable while "
                    "explicit PIC exhibits numerical heating (~2% energy increase). "
                    "From Markidis & Lapenta, J. Comput. Phys. 230 (2011) 7037-7052.",
        container_image="docker://pytorch/pytorch:2.5.1-cuda12.4-cudnn9-runtime",
        entry_command=(
            "pip install numpy scipy matplotlib && "
            "git clone https://github.com/Rise-AGI/magnus-skills.git repo && "
            "cd repo && git checkout markidis2011 && "
            "cd submit_Markidis_2011/reproduction && "
            "MPLBACKEND=Agg python3 fig7_finite_grid_energy.py"
        ),
    )
