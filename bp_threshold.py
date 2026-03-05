from typing import Annotated, Literal, Optional

MagneticField = Annotated[float, {"label": "Magnetic Field Range Max (T)", "description": "Upper limit of B_T range in Tesla for stability threshold plot"}]
RunawayFraction1 = Annotated[float, {"label": "Runaway Fraction 1", "description": "First runaway fraction n_r/n_e"}]
RunawayFraction2 = Annotated[float, {"label": "Runaway Fraction 2", "description": "Second runaway fraction n_r/n_e"}]
RunawayFraction3 = Annotated[float, {"label": "Runaway Fraction 3", "description": "Third runaway fraction n_r/n_e"}]

def blueprint(
    b_max: MagneticField = 4.0,
    nr_ratio_1: RunawayFraction1 = 5e-4,
    nr_ratio_2: RunawayFraction2 = 1e-3,
    nr_ratio_3: RunawayFraction3 = 5e-3,
):
    submit_job(
        task_name="Fulop2006: Stability Threshold (Fig 1)",
        description="Compute stability threshold T_crit vs B_T from Eq.(25) of Fulop et al. 2006 for three runaway fractions",
        container_image="docker://git.pku.edu.cn/2200011363/knowledge-distiller:latest",
        entry_command=f"pip install numpy scipy matplotlib && python3 -c \"\nimport numpy as np\nimport matplotlib; matplotlib.use('Agg')\nimport matplotlib.pyplot as plt\nZ = 1\nB_T = np.linspace(1, {b_max}, 200)\nfor ratio in [{nr_ratio_1}, {nr_ratio_2}, {nr_ratio_3}]:\n    T_crit = (Z**2 * B_T / (20 * ratio))**(2.0/3.0)\n    plt.plot(B_T, T_crit, linewidth=2, label=f'n_r/n_e = {{ratio}}')\nplt.xlabel('B [T]'); plt.ylabel('T_eV')\nplt.legend(); plt.grid(True, alpha=0.3)\nplt.savefig('fig1_threshold.png', dpi=300)\nprint('Done: fig1_threshold.png')\n\"",
    )
