from magnus import submit_job
from typing import Annotated, Literal, Optional

ModelType = Annotated[Literal["cqed", "zd", "both"], {"label": "Model Type", "description": "Which gauge model: truncated cQED, Zd, or both"}]
FinalX = Annotated[float, {"label": "Final x", "description": "Final value of dimensionless parameter x = 1/(ag)^2"}]
MaxT = Annotated[float, {"label": "Max Evolution Time", "description": "Maximum total adiabatic evolution time T"}]

def blueprint(model: ModelType = "both", final_x: FinalX = 10.0, max_T: MaxT = 100.0):
    submit_job(
        task_name="Kuhn2014: Adiabatic Preparation (Figs 3-4)",
        description="Simulate adiabatic preparation of the Schwinger model ground state using cubic ramp x(t) = xF*(t/T)^3. Computes overlap with exact ground state. Reproduces Figures 3-4 from Kuhn, Cirac, Banuls (2014).",
        entry_command="pip install numpy scipy matplotlib && cd reproduction && python3 fig3_adiabatic_cqed.py && python3 fig4_adiabatic_zd.py",
        container_image="docker://pytorch/pytorch:2.5.1-cuda12.4-cudnn9-runtime",
    )
