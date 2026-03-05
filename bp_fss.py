from typing import Annotated, Literal, Optional

ContainerImage = "docker://git.pku.edu.cn/2200011363/knowledge-distiller:latest"

McSteps = Annotated[int, {"label": "MC Steps (L=32)", "description": "Number of measurement sweeps for L=32. Smaller lattices use proportionally more."}]
Seed = Annotated[int, {"label": "Random Seed", "description": "RNG seed for reproducibility."}]

def blueprint(n_mc: McSteps = 4000, seed: Seed = 42):
    submit_job(
        task_name="Kotze 2008: Finite Size Scaling (Figs 19-21)",
        description="Finite size scaling analysis for 2D Ising model: critical exponents, reduced-unit collapse, and Binder cumulant intersection for Tc determination. Lattices L=4,8,16,32.",
        container_image=ContainerImage,
        repo_name="magnus-skills",
        branch="kotze2008",
        entry_command="python3 -m pip install --break-system-packages numpy scipy matplotlib && cd submit_Kotze_2008/reproduction && MPLBACKEND=Agg python3 fig_fss.py",
    )
