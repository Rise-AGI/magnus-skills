from typing import Annotated, Literal, Optional

ScanParam = Annotated[str, {"label": "Scan Parameter", "description": "Which parameter to vary: B (magnetic field), ne (density), or T (temperature)"}]

def blueprint(
    scan_param: ScanParam = "B",
):
    submit_job(
        task_name="Fulop2006: Parameter Scan (Figs 4-6)",
        description="Compute growth rate vs propagation angle for varying plasma parameters (B, n_e, or T). Generates Figs 4, 5, and 6 of Fulop et al. 2006.",
        container_image="docker://git.pku.edu.cn/2200011363/knowledge-distiller:latest",
        entry_command="pip install numpy scipy matplotlib && cd /magnus/workspace/repository/submit_Fulop_2006/reproduction && MPLBACKEND=Agg python3 fig4_angle_B.py && MPLBACKEND=Agg python3 fig5_angle_density.py && MPLBACKEND=Agg python3 fig6_angle_temperature.py && echo 'All done'",
    )
