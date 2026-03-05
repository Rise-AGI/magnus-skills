from typing import Annotated, Literal, Optional

ContainerImage = "docker://git.pku.edu.cn/2200011363/knowledge-distiller:latest"

McSteps = Annotated[int, {"label": "MC Steps (L=2)", "description": "Number of measurement sweeps for L=2 lattice."}]
Seed = Annotated[int, {"label": "Random Seed", "description": "RNG seed for reproducibility."}]

def blueprint(n_mc: McSteps = 100000, seed: Seed = 42):
    submit_job(
        task_name="Kotze 2008: Exact Comparison (Fig 17)",
        description="Compare Monte Carlo simulation of 2x2 Ising lattice against exact analytical solutions for heat capacity and susceptibility.",
        container_image=ContainerImage,
        entry_command=f"pip install numpy scipy matplotlib && cd /tmp && git clone https://github.com/Rise-AGI/magnus-skills.git -b kotze2008 repo && cd repo/submit_Kotze_2008/reproduction && python3 fig17_exact_comparison.py",
    )
