from typing import Annotated, Literal

Action = Annotated[
    Literal["phase_diagram", "scaling_plots", "all"],
    {"label": "Computation", "description": "Which figures to compute: phase_diagram (Fig 9), scaling_plots (Figs 12-16), or all"}
]

def blueprint(action: Action = "all"):
    submit_job(
        task_name="Bernard1997: Phase Diagram & All Figures",
        description="Compute phase diagram and scaling figures from MILC lattice QCD paper",
        container_image="docker://git.pku.edu.cn/2200011363/knowledge-distiller:latest",
        repo_name="sampleproject",
        branch="main",
        commit_sha="621e4974ca25ce531773def586ba3ed8e736b3fc",
        entry_command=(
            "pip install numpy scipy matplotlib && "
            "cd /magnus/workspace/repository/submit_Bernard_1997/reproduction && "
            "MPLBACKEND=Agg python3 fig9_phase_diagram.py && "
            "MPLBACKEND=Agg python3 fig11_polyakov_vs_amps2.py && "
            "MPLBACKEND=Agg python3 fig12_tc_over_mv.py && "
            "MPLBACKEND=Agg python3 fig13_tc_scaling.py && "
            "MPLBACKEND=Agg python3 fig14_tc_vs_a.py && "
            "MPLBACKEND=Agg python3 fig15_r0_sigma.py && "
            "MPLBACKEND=Agg python3 fig16_mv_r0.py && "
            "echo 'All figures computed successfully.'"
        ),
    )
