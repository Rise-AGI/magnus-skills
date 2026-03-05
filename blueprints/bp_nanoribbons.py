from typing import Annotated, Literal

N_width = Annotated[int, {"label": "Ribbon width N", "description": "Number of unit cells across the nanoribbon width"}]

def blueprint(n_width: N_width = 200):
    submit_job(
        task_name="CastroNeto2009: Nanoribbon Spectra",
        description="Compute energy spectra of zigzag and armchair graphene nanoribbons (Fig 17) from Castro Neto et al. Rev. Mod. Phys. 81, 109 (2009).",
        entry_command=(
            "pip install numpy scipy matplotlib && "
            "cd /tmp && "
            "git clone https://github.com/Rise-AGI/magnus-skills.git -b castroneto2009 --depth 1 && "
            "cd magnus-skills/submit_CastroNeto_2009/reproduction && "
            "MPLBACKEND=Agg python3 fig17_nanoribbons.py && "
            "echo 'Nanoribbon spectra computation complete'"
        ),
        container_image="docker://pytorch/pytorch:2.5.1-cuda12.4-cudnn9-runtime",
    )
