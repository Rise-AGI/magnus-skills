from typing import Annotated, Literal, Optional

ContainerImage = "docker://git.pku.edu.cn/2200011363/knowledge-distiller:latest"

McSteps = Annotated[int, {"label": "MC Steps (L=16)", "description": "Number of measurement sweeps for L=16 lattice. Other sizes scale proportionally."}]
EqSteps = Annotated[int, {"label": "Equilibration Sweeps", "description": "Number of equilibration sweeps per temperature."}]
Seed = Annotated[int, {"label": "Random Seed", "description": "RNG seed for reproducibility."}]

def blueprint(n_mc: McSteps = 5000, n_eq: EqSteps = 2000, seed: Seed = 42):
    submit_job(
        task_name="Kotze 2008: Main Observables (Figs 7-9, 14)",
        description="2D Ising Monte Carlo for L=2,4,8,16. Produces energy, heat capacity, magnetization, susceptibility vs temperature.",
        container_image=ContainerImage,
        repo_name="magnus-skills",
        branch="kotze2008",
        entry_command="python3 -m pip install --break-system-packages numpy scipy matplotlib && cd submit_Kotze_2008/reproduction && MPLBACKEND=Agg python3 fig_main_observables.py",
    )
