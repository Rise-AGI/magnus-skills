from typing import Annotated, Literal, Optional

Tolerance = Annotated[str, {"label": "Solver Tolerance", "description": "JFNK solver tolerance (e.g. 1e-7)"}]
NumParticles = Annotated[int, {"label": "Number of Particles", "description": "Number of computational particles (higher = less noise, slower)"}]
NumSteps = Annotated[int, {"label": "Time Steps", "description": "Number of simulation time steps"}]

def blueprint(
    tolerance: Tolerance = "1e-7",
    num_particles: NumParticles = 20000,
    num_steps: NumSteps = 300,
):
    submit_job(
        task_name="EC-PIC: Two-Stream Instability",
        description="Simulate two-stream instability with energy-conserving PIC. "
                    "Compares EC-PIC growth rate with linear theory (gamma=0.35355*wpe) "
                    "and energy conservation with explicit PIC. "
                    "From Markidis & Lapenta, J. Comput. Phys. 230 (2011) 7037-7052.",
        container_image="docker://pytorch/pytorch:2.5.1-cuda12.4-cudnn9-runtime",
        entry_command=(
            "pip install numpy scipy matplotlib && "
            "git clone https://github.com/Rise-AGI/magnus-skills.git repo && "
            "cd repo && git checkout markidis2011 && "
            "cd submit_Markidis_2011/reproduction && "
            f"MPLBACKEND=Agg python3 fig8_twostream_growth.py && "
            f"MPLBACKEND=Agg python3 fig9_twostream_energy.py"
        ),
    )
