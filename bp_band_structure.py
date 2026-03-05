from typing import Annotated, Literal, Optional
from magnus import submit_job

ContainerImage = "docker://git.pku.edu.cn/2200011363/knowledge-distiller:latest"

NFreq = Annotated[int, {"label": "Frequency Points", "description": "Number of nqR values (vertical axis)"}]
NKvec = Annotated[int, {"label": "Wavevector Points", "description": "Number of ka values (horizontal axis)"}]
NArray = Annotated[int, {"label": "Array Size N", "description": "NxN cylinder array for scattering computation"}]
LMax = Annotated[int, {"label": "Max Angular Momentum l", "description": "Max l for Mie expansion"}]
ROverA = Annotated[float, {"label": "R/a Ratio", "description": "Cylinder radius / lattice constant"}]
RefractiveIndex = Annotated[float, {"label": "Refractive Index n", "description": "Refractive index of cylinders"}]
Eta = Annotated[float, {"label": "Smoothing eta", "description": "Broadening parameter for visualization"}]

SCRIPT = r'''
import numpy as np
from scipy.special import jv, jvp, hankel1, h1vp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys, json

p = json.loads(sys.argv[1])
a = 1.0; R = p["R_over_a"] * a; n_refr = p["n_refr"]
l_max = p["l_max"]; N_array = p["N_array"]; eta = p["eta"]

def mie_coeff(l, q):
    nqR = n_refr * q * R; qR = q * R
    num = jvp(l, nqR) * jv(l, qR) - n_refr * jvp(l, qR) * jv(l, nqR)
    den = n_refr * jv(l, nqR) * h1vp(l, qR) - hankel1(l, qR) * jvp(l, nqR)
    return num / den

def lat_sums(q, k_x):
    half_N = N_array // 2
    n1 = np.arange(-half_N, half_N + 1)
    N1, N2 = np.meshgrid(n1, n1); N1, N2 = N1.ravel(), N2.ravel()
    mask = ~((N1 == 0) & (N2 == 0)); X, Y = N1[mask]*a, N2[mask]*a
    R_arr = np.sqrt(X**2+Y**2); phi = np.arctan2(Y,X); qR_a = q*R_arr
    ph = np.exp(1j*k_x*X)*np.exp(-0.003/a*R_arr)
    sums = {}
    for m in range(-2*l_max, 2*l_max+1):
        sums[m] = np.sum(np.exp(1j*m*phi)*hankel1(m,qR_a)*ph)
    return sums

def solve(q, k_x):
    ls = list(range(-l_max,l_max+1)); nl = len(ls)
    s = {l: mie_coeff(l,q) for l in ls}
    S = lat_sums(q, k_x)
    T = np.eye(nl, dtype=complex); rhs = np.zeros(nl, dtype=complex)
    for i,l in enumerate(ls):
        rhs[i] = s[l]*(1j**l)
        for j,lp in enumerate(ls):
            T[i,j] -= s[l]*S.get(l-lp,0)
    return dict(zip(ls, np.linalg.solve(T, rhs)))

nqR_vals = np.linspace(0.5, 7.0, p["n_freq"])
ka_vals = np.linspace(0.05, np.pi, p["n_kvec"])
intensity = np.zeros((len(nqR_vals), len(ka_vals)))

for j, ka in enumerate(ka_vals):
    if j % 10 == 0: print(f"ka {j+1}/{len(ka_vals)}")
    k_x = ka / a
    for i, nqR in enumerate(nqR_vals):
        q = nqR / (n_refr * R)
        try:
            beta = solve(q, k_x)
            b0 = sum(beta[l] for l in range(-l_max,l_max+1))
            b0sq = abs(b0)**2
            intensity[i,j] = eta/(eta**2+1.0/b0sq) if b0sq > 1e-30 else 0
        except Exception: pass

intensity[intensity<=0] = 1e-15
fig, ax = plt.subplots(figsize=(8,6))
ax.pcolormesh(ka_vals, nqR_vals, np.log10(intensity), shading="auto", cmap="hot")
ax.set_xlabel("ka"); ax.set_ylabel("nqR")
ax.set_title(f"Band structure: R/a={p['R_over_a']}, n={n_refr}")
plt.tight_layout()
plt.savefig("band_structure.png", dpi=200)
np.savez("band_structure.npz", nqR=nqR_vals, ka=ka_vals, intensity=intensity)
print("Done. Output: band_structure.png, band_structure.npz")
'''

def blueprint(
    r_over_a: ROverA = 0.35,
    refractive_index: RefractiveIndex = 4.0,
    n_array: NArray = 51,
    l_max: LMax = 2,
    n_freq: NFreq = 80,
    n_kvec: NKvec = 50,
    eta: Eta = 0.1,
):
    import json
    params_json = json.dumps({
        "R_over_a": r_over_a, "n_refr": refractive_index,
        "N_array": n_array, "l_max": l_max,
        "n_freq": n_freq, "n_kvec": n_kvec, "eta": eta,
    })
    submit_job(
        task_name=f"Photonic Bands: R/a={r_over_a}, n={refractive_index}",
        description="Compute photonic band structure via multiple scattering (Juarez et al.)",
        entry_command=f"pip install numpy scipy matplotlib && python3 -c \"import sys; open('run.py','w').write(sys.argv[1])\" '{SCRIPT}' && python3 run.py '{params_json}'",
        container_image=ContainerImage,
    )
