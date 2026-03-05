from typing import Annotated, Literal, Optional

Angle = Annotated[float, {"label": "Propagation Angle (deg)", "description": "Wave propagation angle theta_k in degrees"}]
Density = Annotated[float, {"label": "Electron Density (m^-3)", "description": "Background electron density"}]
MagField = Annotated[float, {"label": "Magnetic Field (T)", "description": "Toroidal magnetic field strength in Tesla"}]
Temperature = Annotated[float, {"label": "Temperature (eV)", "description": "Electron temperature in eV"}]

def blueprint(
    theta_k: Angle = 85.0,
    n_e: Density = 5e19,
    b_t: MagField = 2.0,
    t_ev: Temperature = 10.0,
):
    submit_job(
        task_name="Fulop2006: Dispersion & Growth Rate (Figs 2-3)",
        description="Compute magnetosonic-whistler dispersion relation and growth rate vs wave number at fixed propagation angle",
        container_image="docker://git.pku.edu.cn/2200011363/knowledge-distiller:latest",
        entry_command="pip install numpy scipy matplotlib && cd /magnus/workspace/repository/submit_Fulop_2006/reproduction && MPLBACKEND=Agg python3 fig2_dispersion.py && MPLBACKEND=Agg python3 fig3_growth_rate.py && echo 'All done'",
    )
