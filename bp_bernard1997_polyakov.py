def blueprint():
    submit_job(
        task_name="Bernard1997: Polyakov Loop Comparison",
        description="Compute Polyakov loop vs (aMPS)^2 comparing improved and unimproved Wilson actions (Fig 11)",
        container_image="docker://git.pku.edu.cn/2200011363/knowledge-distiller:latest",
        repo_name="sampleproject",
        branch="main",
        commit_sha="621e4974ca25ce531773def586ba3ed8e736b3fc",
        entry_command=(
            "pip install numpy scipy matplotlib && "
            "cd /magnus/workspace/repository/submit_Bernard_1997/reproduction && "
            "MPLBACKEND=Agg python3 fig11_polyakov_vs_amps2.py && "
            "echo 'Polyakov loop comparison complete.'"
        ),
    )
