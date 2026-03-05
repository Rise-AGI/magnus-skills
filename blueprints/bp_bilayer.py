from typing import Annotated, Literal

V_bias = Annotated[float, {"label": "Interlayer bias V (eV)", "description": "Half the potential difference between layers in eV"}]

def blueprint(v_bias: V_bias = 0.1):
    submit_job(
        task_name="CastroNeto2009: Bilayer Band Structure",
        description="Compute bilayer graphene band structure with/without interlayer bias (Figs 10-11) from Castro Neto et al. Rev. Mod. Phys. 81, 109 (2009).",
        repo_name="magnus-skills",
        namespace="Rise-AGI",
        branch="castroneto2009",
        entry_command=(
            "pip install numpy scipy matplotlib && "
            "cd submit_CastroNeto_2009/reproduction && "
            "MPLBACKEND=Agg python3 fig10_11_bilayer.py && "
            "echo 'Bilayer band structure computation complete'"
        ),
        container_image="docker://pytorch/pytorch:2.5.1-cuda12.4-cudnn9-runtime",
    )
