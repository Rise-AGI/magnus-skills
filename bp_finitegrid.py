from typing import Annotated

NumParticles = Annotated[int, {"label": "Particles", "description": "Computational particles"}]

def blueprint(num_particles: NumParticles = 2000):
    submit_job(
        task_name="EC-PIC: Finite Grid Instability",
        description="Test EC-PIC against finite grid instability in Maxwellian plasma",
        repo_name="magnus-skills",
        namespace="Rise-AGI",
        branch="markidis2011",
        commit_sha="008166e4eef080edad8463a17fa7732a8b176355",
        container_image="docker://pytorch/pytorch:2.5.1-cuda12.4-cudnn9-runtime",
        entry_command="pip install numpy scipy matplotlib && cd submit_Markidis_2011/reproduction && MPLBACKEND=Agg python3 fig7_finite_grid_energy.py",
    )
