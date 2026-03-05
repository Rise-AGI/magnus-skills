from magnus import submit_job
from typing import Annotated, Literal, Optional

ModelType = Annotated[Literal["cqed", "zd", "both"], {"label": "Model Type", "description": "Which gauge model: truncated cQED, Zd, or both"}]
MaxLambda = Annotated[float, {"label": "Max Noise Strength", "description": "Maximum noise strength parameter lambda"}]

def blueprint(model: ModelType = "both", max_lambda: MaxLambda = 0.003):
    submit_job(
        task_name="Kuhn2014: Noise Effects (Figs 5-6)",
        description="Study effect of gauge invariance breaking noise on adiabatic preparation. Computes penalty energy and overlap vs noise strength. Reproduces Figures 5-6 from Kuhn, Cirac, Banuls (2014).",
        entry_command="pip install numpy scipy matplotlib && cd reproduction && python3 fig5_noise_cqed.py && python3 fig6_noise_zd.py",
        container_image="docker://pytorch/pytorch:2.5.1-cuda12.4-cudnn9-runtime",
    )
