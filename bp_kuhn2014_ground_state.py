from magnus import submit_job
from typing import Annotated, Literal, Optional

ModelType = Annotated[Literal["cqed", "zd", "both"], {"label": "Model Type", "description": "Which gauge model: truncated cQED, Zd, or both"}]
MaxD = Annotated[int, {"label": "Max Link Dimension", "description": "Maximum link Hilbert space dimension d (odd: 3,5,7,9)"}]
MaxN = Annotated[int, {"label": "Max System Size", "description": "Maximum number of lattice sites N (even)"}]

def blueprint(model: ModelType = "both", max_d: MaxD = 7, max_n: MaxN = 8):
    submit_job(
        task_name="Kuhn2014: Ground State Properties (Figs 1-2)",
        description="Compute ground state energy density and energy gap for the Schwinger model with finite-dimensional link variables using exact diagonalization. Reproduces Figures 1 and 2 from Kuhn, Cirac, Banuls (2014).",
        entry_command="pip install numpy scipy matplotlib && cd reproduction && python3 fig1_energy_density.py && python3 fig2_energy_gap.py",
        container_image="docker://pytorch/pytorch:2.5.1-cuda12.4-cudnn9-runtime",
    )
