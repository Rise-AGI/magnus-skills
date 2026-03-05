from typing import Annotated, Literal, Optional

ContainerImage = "docker://git.pku.edu.cn/2200011363/knowledge-distiller:latest"

McSteps = Annotated[int, {"label": "MC Steps (L=32)", "description": "Number of measurement sweeps for L=32. Smaller lattices use proportionally more."}]
Seed = Annotated[int, {"label": "Random Seed", "description": "RNG seed for reproducibility."}]

def blueprint(n_mc: McSteps = 4000, seed: Seed = 42):
    submit_job(
        task_name="Kotze 2008: Finite Size Scaling (Figs 19-21)",
        description="Finite size scaling analysis for 2D Ising model: critical exponents, reduced-unit collapse, and Binder cumulant intersection for Tc determination. Lattices L=4,8,16,32.",
        container_image=ContainerImage,
        entry_command=f"pip install numpy scipy matplotlib && cd /tmp && git clone https://github.com/Rise-AGI/magnus-skills.git -b kotze2008 repo && cd repo/submit_Kotze_2008/reproduction && python3 fig_fss.py",
    )
