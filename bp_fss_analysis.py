from magnus import submit_job, JobType
from typing import Annotated, Literal, Optional

LMin = Annotated[int, {"label": "Minimum Lattice Size", "description": "Smallest L included in the FSS fit (default 16)"}]
NumCorrections = Annotated[Literal["1", "2", "3"], {"label": "Number of Corrections", "description": "Number of fixed correction exponents: 1 (omega1=0.83), 2 (+omega2=4), or 3 (+omega_nu=1.6)"}]

def blueprint(l_min: LMin = 16, num_corrections: NumCorrections = "3"):
    submit_job(
        task_name="Ferrenberg2018: FSS Analysis (nu and Kc)",
        description=f"Finite-size scaling analysis of 3D Ising critical exponent nu and coupling Kc. "
                    f"Plots nu and Kc vs L_min with {num_corrections} correction exponent(s), "
                    f"starting from L_min={l_min}. Reproduces Figs 2-4 of Ferrenberg et al. (2018).",
        repo_name="magnus-skills",
        namespace="Rise-AGI",
        branch="ferrenberg2018",
        commit_sha="9432b9250c71470e62c7b71eaa9c420b1b02969b",
        entry_command=(
            "cd submit_Ferrenberg_2018 && python3 -m pip install numpy scipy matplotlib && "
            "cd reproduction && export MPLBACKEND=Agg && "
            "python3 fig2_nu_scaling.py && "
            "python3 fig3_kc_one_correction.py && "
            "python3 fig4_kc_scaling.py"
        ),
        container_image="docker://pytorch/pytorch:2.5.1-cuda12.4-cudnn9-runtime",
    )
