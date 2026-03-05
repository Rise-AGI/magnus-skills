from typing import Annotated, Literal, Optional
from magnus import submit_job

ContainerImage = "docker://git.pku.edu.cn/2200011363/knowledge-distiller:latest"

LMax = Annotated[int, {"label": "Max Angular Momentum l", "description": "Maximum angular momentum quantum number"}]
NFreq = Annotated[int, {"label": "Frequency Points", "description": "Number of nqR values to scan"}]
ROverA = Annotated[float, {"label": "R/a Ratio", "description": "Cylinder radius / lattice constant"}]
RefractiveIndex = Annotated[float, {"label": "Refractive Index n", "description": "Refractive index of cylinders"}]

def blueprint(
    r_over_a: ROverA = 0.35,
    refractive_index: RefractiveIndex = 4.0,
    l_max: LMax = 2,
    n_freq: NFreq = 300,
):
    script = """
import numpy as np
from scipy.special import jv, jvp, hankel1, h1vp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

R_over_a = """ + str(r_over_a) + """
n_refr = """ + str(refractive_index) + """
l_max = """ + str(l_max) + """
n_freq = """ + str(n_freq) + """
a = 1.0
R = R_over_a * a

def mie_coeff(l, q):
    nqR = n_refr*q*R; qR = q*R
    num = jvp(l,nqR)*jv(l,qR) - n_refr*jvp(l,qR)*jv(l,nqR)
    den = n_refr*jv(l,nqR)*h1vp(l,qR) - hankel1(l,qR)*jvp(l,nqR)
    return num / den

nqR_array = np.linspace(0.3, 8.0, n_freq)
b_single = {}
for l in range(l_max+1):
    b_single[l] = np.zeros(n_freq)

for idx, nqR in enumerate(nqR_array):
    q = nqR / (n_refr * R)
    for l in range(l_max+1):
        s_l = mie_coeff(l, q)
        b_single[l][idx] = abs(s_l)

fig, ax = plt.subplots(figsize=(8, 5))
colors = ["b", "r", "g", "m", "c"]
for l in range(l_max+1):
    ax.plot(nqR_array, b_single[l], colors[l%5], label="l="+str(l))
ax.set_xlabel("nqR"); ax.set_ylabel("|s_l|")
ax.set_title("Mie resonances")
ax.legend(); ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("mie_resonances.png", dpi=200)
data = np.column_stack([nqR_array] + [b_single[l] for l in range(l_max+1)])
np.savetxt("mie_resonances.csv", data, delimiter=",", fmt="%.8e")
print("Done")
"""
    submit_job(
        task_name="Mie Resonances: R/a=" + str(r_over_a) + ", n=" + str(refractive_index),
        description="Compute single-cylinder Mie scattering coefficients vs frequency",
        entry_command="pip install numpy scipy matplotlib && cat << 'PYEOF' > run.py\n" + script + "\nPYEOF\npython3 run.py",
        container_image=ContainerImage,
    )
