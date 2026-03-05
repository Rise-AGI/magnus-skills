from typing import Annotated

Tolerance = Annotated[str, {"label": "Solver Tolerance", "description": "JFNK solver tolerance"}]

def blueprint(tolerance: Tolerance = "1e-7"):
    submit_job(
        task_name="EC-PIC: Two-Stream Instability",
        description="Two-stream instability with energy-conserving PIC",
        repo_name="magnus-skills",
        namespace="Rise-AGI",
        branch="markidis2011",
        commit_sha="8f13d2bb33f429ea85692d36d30541a27be576b1",
        container_image="docker://pytorch/pytorch:2.5.1-cuda12.4-cudnn9-runtime",
        entry_command="pip install numpy scipy matplotlib && cd submit_Markidis_2011/reproduction && MPLBACKEND=Agg python3 fig8_twostream_growth.py && MPLBACKEND=Agg python3 fig9_twostream_energy.py",
    )
