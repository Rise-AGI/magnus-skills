from typing import Annotated

ContainerImage = "docker://pytorch/pytorch:2.5.1-cuda12.4-cudnn9-runtime"

def blueprint():
    submit_job(
        task_name="Bret2009: Run All Reproductions",
        description="Run all figure reproduction scripts for Bret (2009): parallel growth rates (Fig 1), instability hierarchy maps (Fig 5).",
        entry_command="pip install numpy scipy matplotlib && cd /magnus/workspace/repository/submit_Bret_2009/reproduction && python3 run_all.py",
        container_image=ContainerImage,
    )
