from typing import Annotated

Anisotropy = Annotated[float, {"label": "Anisotropy", "description": "Temperature anisotropy a = vthy^2/vthx^2 - 1"}]

def blueprint(anisotropy: Anisotropy = 15.0):
    submit_job(
        task_name="EM-PIC: Weibel Instability",
        description="Weibel instability simulation with electromagnetic PIC",
        repo_name="magnus-skills",
        namespace="Rise-AGI",
        branch="markidis2011",
        commit_sha="008166e4eef080edad8463a17fa7732a8b176355",
        container_image="docker://pytorch/pytorch:2.5.1-cuda12.4-cudnn9-runtime",
        entry_command="pip install numpy scipy matplotlib && cd submit_Markidis_2011/reproduction && MPLBACKEND=Agg python3 fig10_11_weibel.py",
    )
