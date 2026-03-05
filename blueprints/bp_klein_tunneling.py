from typing import Annotated, Literal

E_barrier = Annotated[float, {"label": "Electron energy (eV)", "description": "Energy of incident Dirac electrons in eV"}]
V0_barrier = Annotated[float, {"label": "Barrier height V0 (eV)", "description": "Square potential barrier height in eV"}]

def blueprint(energy: E_barrier = 0.080, v0: V0_barrier = 0.200):
    submit_job(
        task_name="CastroNeto2009: Klein Tunneling",
        description="Compute Klein tunneling transmission T(phi) for Dirac electrons through a square barrier (Fig 7) from Castro Neto et al. Rev. Mod. Phys. 81, 109 (2009).",
        entry_command=(
            "pip install numpy scipy matplotlib && "
            "cd /tmp && "
            "git clone https://github.com/Rise-AGI/magnus-skills.git -b castroneto2009 --depth 1 && "
            "cd magnus-skills/submit_CastroNeto_2009/reproduction && "
            "MPLBACKEND=Agg python3 fig7_klein_tunneling.py && "
            "echo 'Klein tunneling computation complete'"
        ),
        container_image="docker://pytorch/pytorch:2.5.1-cuda12.4-cudnn9-runtime",
    )
