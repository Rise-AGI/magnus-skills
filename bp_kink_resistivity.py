from magnus import submit_job, JobType
from typing import Annotated

REFrac1 = Annotated[float, {"label": "RE Fraction 1", "description": "First RE current fraction (0 to 1)"}]
REFrac2 = Annotated[float, {"label": "RE Fraction 2", "description": "Second RE current fraction (0 to 1)"}]
REFrac3 = Annotated[float, {"label": "RE Fraction 3", "description": "Third RE current fraction (0 to 1)"}]

def blueprint(
    f_re_1: REFrac1 = 0.0,
    f_re_2: REFrac2 = 0.5,
    f_re_3: REFrac3 = 1.0,
):
    script = (
        'import numpy as np\n'
        'import matplotlib; matplotlib.use("Agg")\n'
        'import matplotlib.pyplot as plt\n'
        'S_inv = np.logspace(-7, -2, 200)\n'
        'fractions = [' + str(f_re_1) + ', ' + str(f_re_2) + ', ' + str(f_re_3) + ']\n'
        'labels = ["I_r/I_p = " + str(f) for f in fractions]\n'
        'colors = ["blue", "red", "green"]\n'
        'styles = ["-", "--", ":"]\n'
        'fig, ax = plt.subplots(figsize=(8, 6))\n'
        'all_data = [S_inv]\n'
        'col_names = ["S_inv"]\n'
        'for f_RE, label, color, style in zip(fractions, labels, colors, styles):\n'
        '    S_eff = np.maximum(S_inv * (1 - f_RE), 1e-20)\n'
        '    gamma = np.where(S_eff < 1e-3, 0.5 * S_eff**(1./3.), 0.3 * S_eff**(1./3.) + 0.1 * S_eff**(2./3.))\n'
        '    ax.loglog(S_inv, gamma, color=color, ls=style, lw=2, label=label)\n'
        '    all_data.append(gamma)\n'
        '    col_names.append("gamma_fRE_" + str(f_RE))\n'
        'ax.set_xlabel("S^{-1}"); ax.set_ylabel("gamma * tau_A")\n'
        'ax.set_title("Kink mode growth rate vs resistivity"); ax.legend(); ax.grid(True, alpha=0.3, which="both")\n'
        'ax.set_xlim(1e-7, 1e-2)\n'
        'fig.savefig("fig4b_kink_scaling.png", dpi=150); plt.close()\n'
        'np.savetxt("fig4b_kink_scaling.csv", np.column_stack(all_data), delimiter=",",'
        ' header=",".join(col_names), comments="", fmt="%.8e")\n'
        'psi_N = np.linspace(0, 2.0, 500)\n'
        'eta_ratio = np.ones_like(psi_N)\n'
        'for i, pn in enumerate(psi_N):\n'
        '    if pn <= 0.8: eta_ratio[i] = 1.0\n'
        '    elif pn <= 1.0: eta_ratio[i] = 1.0 + 2.0*((pn-0.8)/0.2)**2\n'
        '    elif pn <= 1.5: eta_ratio[i] = 3.0\n'
        '    else: eta_ratio[i] = 3.0 * np.exp(5*(pn-1.5))\n'
        'np.savetxt("fig5_resistivity.csv", np.column_stack([psi_N, eta_ratio]),'
        ' delimiter=",", header="psi_N,eta_over_eta_axis", comments="", fmt="%.8e")\n'
        'fig2, ax2 = plt.subplots(figsize=(8, 5))\n'
        'ax2.semilogy(psi_N, eta_ratio, "b-", lw=2)\n'
        'ax2.axvline(x=1.0, color="gray", ls="--", alpha=0.5, label="LCFS")\n'
        'ax2.set_xlabel("psi_N"); ax2.set_ylabel("eta/eta_axis")\n'
        'ax2.set_title("Resistivity profile"); ax2.legend(); ax2.grid(True, alpha=0.3)\n'
        'fig2.savefig("fig5_resistivity.png", dpi=150); plt.close()\n'
        'print("Done: fig4b_kink_scaling and fig5_resistivity generated")\n'
    )

    cmd = "python3 -m pip install --break-system-packages numpy matplotlib && python3 -c '" + script.replace("'", "'\\''") + "'"

    submit_job(
        task_name="Bandaru2018: Kink Scaling & Resistivity Profile",
        description="Kink mode scaling (Fig 4b) and resistivity profile (Fig 5)",
        namespace="Rise-AGI",
        repo_name="magnus-skills",
        commit_sha="HEAD",
        entry_command=cmd,
        container_image="docker://git.pku.edu.cn/2200011363/knowledge-distiller:latest",
    )
