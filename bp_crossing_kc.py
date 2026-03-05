from magnus import submit_job, JobType
from typing import Annotated, Literal, Optional

def blueprint():
    submit_job(
        task_name="Ferrenberg2018: Cumulant Crossing K_c",
        description="Plot K_c from cumulant crossing technique vs L_min with 1 and 2 "
                    "correction terms. Reproduces Fig 6 of Ferrenberg et al. (2018).",
        repo_name="magnus-skills",
        namespace="Rise-AGI",
        branch="ferrenberg2018",
        commit_sha="9432b9250c71470e62c7b71eaa9c420b1b02969b",
        entry_command=(
            "cd submit_Ferrenberg_2018 && python3 -m pip install numpy scipy matplotlib && "
            "cd reproduction && export MPLBACKEND=Agg && "
            "python3 fig6_kc_crossing.py"
        ),
        container_image="docker://pytorch/pytorch:2.5.1-cuda12.4-cudnn9-runtime",
    )
