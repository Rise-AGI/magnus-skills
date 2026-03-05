from typing import Annotated, Literal, Optional

ContainerImage = "docker://pytorch/pytorch:2.5.1-cuda12.4-cudnn9-runtime"

LatticeSize = Annotated[int, {"label": "Lattice Size L", "description": "Linear size of L^3 cubic lattice"}]
DisorderW = Annotated[float, {"label": "Disorder W", "description": "Disorder strength"}]
DisorderType = Annotated[Literal["box", "gaussian"], {"label": "Disorder Type", "description": "Distribution type"}]
NumRealizations = Annotated[int, {"label": "Realizations", "description": "Number of disorder realizations"}]
Regime = Annotated[Literal["metallic", "localized", "critical"], {"label": "Regime", "description": "Transport regime"}]

def blueprint(
    L: LatticeSize = 10,
    W: DisorderW = 2.0,
    disorder: DisorderType = "gaussian",
    n_realizations: NumRealizations = 500,
    regime: Regime = "metallic"
):
    submit_job(
        task_name=f"Markos2006: Level Statistics ({regime})",
        description=f"Compute level spacing distribution p(s) for 3D Anderson model ({regime} regime, L={L}, W={W}, {disorder})",
        repo_name="magnus-skills",
        commit_sha="HEAD",
        container_image=ContainerImage,
        entry_command=f"""pip install numpy scipy matplotlib && python3 -c "
import numpy as np
from scipy.linalg import eigvalsh
from scipy import sparse
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

L={L}; W={W}; disorder='{disorder}'; n_real={n_realizations}; regime='{regime}'
rng_base = 42

all_spacings = []
for r in range(n_real):
    rng = np.random.default_rng(rng_base + r)
    N = L**3
    if disorder == 'box':
        diag = W * rng.uniform(-0.5, 0.5, size=N)
    else:
        diag = W * rng.normal(0, 1, size=N)
    row, col, val = [], [], []
    for x in range(L):
        for y in range(L):
            for z in range(L):
                i = (x%L)*L*L + (y%L)*L + (z%L)
                for dx,dy,dz in [(1,0,0),(0,1,0),(0,0,1)]:
                    j = ((x+dx)%L)*L*L + ((y+dy)%L)*L + ((z+dz)%L)
                    row.extend([i,j]); col.extend([j,i]); val.extend([1.0,1.0])
    H = sparse.csr_matrix((val,(row,col)),shape=(N,N))
    H.setdiag(diag)
    eigs = eigvalsh(H.toarray())
    E = np.sort(eigs)
    mask = (E >= -0.5) & (E <= 0.5)
    E_sel = E[mask]
    s = np.diff(E_sel)
    if len(s) > 0:
        s = s / np.mean(s)
        all_spacings.extend(s)
    if (r+1) % 100 == 0:
        print(f'Realization {{r+1}}/{{n_real}}')

all_spacings = np.array(all_spacings)
print(f'Total spacings: {{len(all_spacings)}}')

bins = np.linspace(0, 5, 100)
centers = 0.5*(bins[:-1]+bins[1:])
hist, _ = np.histogram(all_spacings, bins=bins, density=True)
s_th = np.linspace(0, 5, 300)
wigner = (np.pi/2)*s_th*np.exp(-np.pi*s_th**2/4)
poisson = np.exp(-s_th)
np.savetxt(f'level_spacing_{{regime}}.csv', np.column_stack([centers, hist]),
           delimiter=',', header='s,p_s', fmt='%.8e', comments='')

fig, ax = plt.subplots(figsize=(8,6))
ax.bar(centers, hist, width=centers[1]-centers[0], alpha=0.6, label=f'{{disorder}} W={{W}}, L={{L}}')
ax.plot(s_th, wigner, 'k-', lw=2, label='Wigner (GOE)')
ax.plot(s_th, poisson, 'k--', lw=2, label='Poisson')
ax.set_xlabel('s'); ax.set_ylabel('p(s)')
ax.set_title(f'Level spacing: {{regime}} regime')
ax.legend(); ax.set_xlim(0, 4); plt.tight_layout()
plt.savefig(f'level_spacing_{{regime}}.png', dpi=150)
print('Done.')
"
""",
    )
