from typing import Annotated, Literal, Optional

ContainerImage = "docker://pytorch/pytorch:2.5.1-cuda12.4-cudnn9-runtime"

LatticeSize = Annotated[int, {"label": "Lattice Size L", "description": "Linear size of L^3 cubic lattice"}]
DisorderStrength = Annotated[float, {"label": "Disorder W", "description": "Disorder strength parameter"}]
DisorderType = Annotated[Literal["box", "gaussian"], {"label": "Disorder Type", "description": "box or gaussian distribution"}]
NumRealizations = Annotated[int, {"label": "Realizations", "description": "Number of disorder realizations"}]

def blueprint(
    L: LatticeSize = 10,
    W: DisorderStrength = 5.0,
    disorder: DisorderType = "box",
    n_realizations: NumRealizations = 5
):
    submit_job(
        task_name="Markos2006: Density of States",
        description=f"Compute density of states for 3D Anderson model (L={L}, W={W}, {disorder})",
        repo_name="magnus-skills",
        commit_sha="HEAD",
        container_image=ContainerImage,
        entry_command=f"""pip install numpy scipy matplotlib && python3 -c "
import numpy as np
from scipy.linalg import eigvalsh
from scipy import sparse

L={L}; W={W}; disorder='{disorder}'; n_real={n_realizations}
rng = np.random.default_rng(42)
all_eigs = []
for r in range(n_real):
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
    all_eigs.extend(eigs)
    print(f'Realization {{r+1}}/{{n_real}} done')
all_eigs = np.array(all_eigs)
counts, edges = np.histogram(all_eigs, bins=200, range=(-15,15), density=True)
centers = 0.5*(edges[:-1]+edges[1:])
np.savetxt('dos_result.csv', np.column_stack([centers, counts]), delimiter=',', header='energy,dos', fmt='%.8e', comments='')
print(f'Done. {{len(all_eigs)}} eigenvalues computed.')
"
""",
    )
