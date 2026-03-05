from typing import Annotated, Literal, Optional
from magnus import submit_job

ContainerImage = "docker://git.pku.edu.cn/2200011363/knowledge-distiller:latest"

LMax = Annotated[int, {"label": "Max Angular Momentum l", "description": "Maximum angular momentum quantum number for Mie expansion"}]
NArray = Annotated[int, {"label": "Array Size N", "description": "NxN square array of cylinders"}]
NFreq = Annotated[int, {"label": "Frequency Points", "description": "Number of frequency points to compute"}]
ROverA = Annotated[float, {"label": "R/a Ratio", "description": "Cylinder radius / lattice constant"}]
RefractiveIndex = Annotated[float, {"label": "Refractive Index n", "description": "Refractive index of cylinders"}]
KFactor = Annotated[float, {"label": "k/q Factor", "description": "Ratio of Bloch wavevector to free-space wavevector"}]

SCRIPT = r'''
import numpy as np
from scipy.special import jv, jvp, hankel1, h1vp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys, json

params = json.loads(sys.argv[1])
R_over_a = params["R_over_a"]
n_refr = params["n_refr"]
N_array = params["N_array"]
l_max = params["l_max"]
n_freq = params["n_freq"]
k_factor = params["k_factor"]

a = 1.0
R = R_over_a * a

def mie_scattering_coeff(l, n, q, R):
    nqR = n * q * R; qR = q * R
    num = jvp(l, nqR) * jv(l, qR) - n * jvp(l, qR) * jv(l, nqR)
    den = n * jv(l, nqR) * h1vp(l, qR) - hankel1(l, qR) * jvp(l, nqR)
    return num / den

def lattice_sums_vectorized(q, k_x, N, a, l_max, damping=0.0):
    half_N = N // 2
    n1 = np.arange(-half_N, half_N + 1)
    N1, N2 = np.meshgrid(n1, n1)
    N1, N2 = N1.ravel(), N2.ravel()
    mask = ~((N1 == 0) & (N2 == 0))
    X, Y = N1[mask] * a, N2[mask] * a
    R_arr = np.sqrt(X**2 + Y**2)
    phi_arr = np.arctan2(Y, X)
    qR_arr = q * R_arr
    phase_arr = np.exp(1j * k_x * X) * np.exp(-damping * R_arr)
    sums = {}
    for m in range(-2 * l_max, 2 * l_max + 1):
        sums[m] = np.sum(np.exp(1j * m * phi_arr) * hankel1(m, qR_arr) * phase_arr)
    return sums

def solve_bloch(q, k_x):
    ls = list(range(-l_max, l_max + 1)); n_l = len(ls)
    s = {l: mie_scattering_coeff(l, n_refr, q, R) for l in ls}
    lsums = lattice_sums_vectorized(q, k_x, N_array, a, l_max, 0.005 / a)
    T_mat = np.eye(n_l, dtype=complex)
    rhs = np.zeros(n_l, dtype=complex)
    for i, l in enumerate(ls):
        rhs[i] = s[l] * (1j ** l)
        for j, lp in enumerate(ls):
            T_mat[i, j] -= s[l] * lsums.get(l - lp, 0)
    return dict(zip(ls, np.linalg.solve(T_mat, rhs)))

nqR_array = np.linspace(0.5, 8.0, n_freq)
sum_imag = np.zeros(n_freq)
b0_sq = np.zeros(n_freq)

for idx, nqR in enumerate(nqR_array):
    q = nqR / (n_refr * R)
    k_x = k_factor * q
    try:
        beta = solve_bloch(q, k_x)
        for l in range(-l_max, l_max + 1):
            sum_imag[idx] += beta[l].imag
        b0_sq[idx] = abs(sum(beta[l] for l in range(-l_max, l_max + 1)))**2
    except Exception:
        pass

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 10))
ax1.plot(nqR_array, sum_imag, "b-")
ax1.set_xlabel("nqR"); ax1.set_ylabel("Sum Im(b_0l)")
ax1.set_title(f"Scattered field: R/a={R_over_a}, n={n_refr}, {N_array}x{N_array}")
ax2.semilogy(nqR_array, b0_sq, "b-")
ax2.set_xlabel("nqR"); ax2.set_ylabel("|b_0|^2")
ax2.set_title("Scattering amplitude")
plt.tight_layout()
plt.savefig("mie_scattering_results.png", dpi=200)

np.savetxt("mie_scattering_data.csv",
    np.column_stack([nqR_array, sum_imag, b0_sq]),
    header="nqR,sum_Im_b0l,b0_squared", delimiter=",", fmt="%.8e")
print("Done. Output: mie_scattering_results.png, mie_scattering_data.csv")
'''

def blueprint(
    r_over_a: ROverA = 0.1,
    refractive_index: RefractiveIndex = 10.0,
    n_array: NArray = 101,
    l_max: LMax = 1,
    n_freq: NFreq = 200,
    k_factor: KFactor = 1.01,
):
    import json
    params_json = json.dumps({
        "R_over_a": r_over_a, "n_refr": refractive_index,
        "N_array": n_array, "l_max": l_max, "n_freq": n_freq,
        "k_factor": k_factor,
    })
    submit_job(
        task_name=f"Mie Scattering: R/a={r_over_a}, n={refractive_index}, {n_array}x{n_array}",
        description="Compute Mie scattering coefficients for an array of dielectric cylinders using the multiple scattering approach (Juarez et al.)",
        entry_command=f"pip install numpy scipy matplotlib && python3 -c \"import sys; open('run.py','w').write(sys.argv[1])\" '{SCRIPT}' && python3 run.py '{params_json}'",
        container_image=ContainerImage,
    )
