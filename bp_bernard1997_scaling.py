from typing import Annotated, Literal

BetaChoice = Annotated[
    Literal["6.40", "6.60", "6.80", "7.20", "all"],
    {"label": "Beta value", "description": "Which coupling constant beta to analyze (or all)"}
]

def blueprint(beta: BetaChoice = "all"):
    submit_job(
        task_name="Bernard1997: Scaling Analysis",
        description="Compute Tc scaling ratios (Figs 12-16) for improved Wilson lattice QCD",
        container_image="docker://git.pku.edu.cn/2200011363/knowledge-distiller:latest",
        repo_name="sampleproject",
        branch="main",
        commit_sha="621e4974ca25ce531773def586ba3ed8e736b3fc",
        entry_command=(
            "pip install numpy scipy matplotlib && "
            "cd /magnus/workspace/repository/submit_Bernard_1997/reproduction && "
            "MPLBACKEND=Agg python3 fig12_tc_over_mv.py && "
            "MPLBACKEND=Agg python3 fig13_tc_scaling.py && "
            "MPLBACKEND=Agg python3 fig14_tc_vs_a.py && "
            "MPLBACKEND=Agg python3 fig15_r0_sigma.py && "
            "MPLBACKEND=Agg python3 fig16_mv_r0.py && "
            "echo 'Scaling analysis complete.'"
        ),
    )
