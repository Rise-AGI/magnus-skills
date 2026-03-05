from typing import Annotated, Literal

T_hop = Annotated[float, {"label": "Hopping energy t (eV)", "description": "Nearest-neighbor hopping parameter in eV"}]
T_prime_ratio = Annotated[float, {"label": "t'/t ratio", "description": "Next-nearest-neighbor hopping ratio t'/t"}]

def blueprint(t: T_hop = 2.7, tp_ratio: T_prime_ratio = 0.2):
    submit_job(
        task_name="CastroNeto2009: Band Structure & DOS",
        description="Compute graphene tight-binding band structure (Fig 3) and density of states (Fig 5) from Castro Neto et al. Rev. Mod. Phys. 81, 109 (2009).",
        repo_name="magnus-skills",
        namespace="Rise-AGI",
        branch="castroneto2009",
        entry_command=(
            "pip install numpy scipy matplotlib && "
            "cd submit_CastroNeto_2009/reproduction && "
            "MPLBACKEND=Agg python3 fig3_band_structure.py && "
            "MPLBACKEND=Agg python3 fig5_dos.py && "
            "echo 'Band structure and DOS computation complete'"
        ),
        container_image="docker://pytorch/pytorch:2.5.1-cuda12.4-cudnn9-runtime",
    )
