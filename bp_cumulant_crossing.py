from magnus import submit_job, JobType
from typing import Annotated, Literal, Optional

LatticeSize = Annotated[int, {"label": "Max Lattice Size", "description": "Largest lattice L to simulate (default 32). Larger = slower."}]
NMeasure = Annotated[int, {"label": "Measurements", "description": "Number of Wolff steps for measurement (default 50000)"}]

def blueprint(max_L: LatticeSize = 32, n_measure: NMeasure = 50000):
    submit_job(
        task_name="Ferrenberg2018: Binder Cumulant Crossing",
        description=f"Wolff cluster MC simulation for 3D Ising model up to L={max_L}. "
                    f"Computes U_4 vs K via histogram reweighting and U_4 vs L "
                    f"self-consistency check. Reproduces Figs 5, 7 of Ferrenberg et al. (2018).",
        entry_command=(
            "pip install numpy scipy matplotlib && "
            "cd reproduction && export MPLBACKEND=Agg && "
            "python3 fig5_cumulant_crossing.py && "
            "python3 fig7_u4_consistency.py"
        ),
        container_image="docker://git.pku.edu.cn/2200011363/knowledge-distiller:latest",
    )
