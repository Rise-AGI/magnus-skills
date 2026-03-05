from magnus import submit_job, JobType
from typing import Annotated, Literal, Optional

def blueprint():
    submit_job(
        task_name="Ferrenberg2018: Cumulant Crossing K_c",
        description="Plot K_c from cumulant crossing technique vs L_min with 1 and 2 "
                    "correction terms. Reproduces Fig 6 of Ferrenberg et al. (2018).",
        entry_command=(
            "pip install numpy scipy matplotlib && "
            "cd reproduction && export MPLBACKEND=Agg && "
            "python3 fig6_kc_crossing.py"
        ),
        container_image="docker://git.pku.edu.cn/2200011363/knowledge-distiller:latest",
    )
