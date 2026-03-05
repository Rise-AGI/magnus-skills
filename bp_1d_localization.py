from typing import Annotated, Literal, Optional

ContainerImage = "docker://pytorch/pytorch:2.5.1-cuda12.4-cudnn9-runtime"

SystemLength = Annotated[int, {"label": "Chain Length Lz", "description": "Length of the 1D disordered chain"}]
DisorderW = Annotated[float, {"label": "Disorder W", "description": "Disorder strength"}]
Energy = Annotated[float, {"label": "Energy E", "description": "Electron energy"}]
NumSamples = Annotated[int, {"label": "Samples", "description": "Number of disorder realizations"}]

def blueprint(
    Lz: SystemLength = 300,
    W: DisorderW = 1.0,
    E: Energy = 1.0,
    n_samples: NumSamples = 50000
):
    submit_job(
        task_name="Markos2006: 1D Localization",
        description=f"Compute 1D Anderson localization: p(x) distribution and mean <x> vs Lz",
        repo_name="magnus-skills",
        commit_sha="HEAD",
        container_image=ContainerImage,
        entry_command=f"""pip install numpy scipy matplotlib && python3 -c "
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

Lz={Lz}; W={W}; E={E}; n_samples={n_samples}
rng = np.random.default_rng(42)
x_vals = np.zeros(n_samples)
for i in range(n_samples):
    eps = W * rng.uniform(-0.5, 0.5, size=Lz)
    T = np.eye(2)
    for n in range(Lz):
        Tn = np.array([[E - eps[n], -1.0], [1.0, 0.0]])
        T = Tn @ T
    TtT = T.T @ T
    cosh2 = max(1.0, (np.trace(TtT) + 2) / 4.0)
    x_vals[i] = 2 * np.arccosh(np.sqrt(cosh2))
    if (i+1) % 10000 == 0:
        print(f'Sample {{i+1}}/{{n_samples}}')

mean_x = np.mean(x_vals)
var_x = np.var(x_vals)
print(f'<x> = {{mean_x:.4f}}, var(x) = {{var_x:.4f}}')

bins = np.linspace(0, 30, 150)
hist, edges = np.histogram(x_vals, bins=bins, density=True)
centers = 0.5*(edges[:-1]+edges[1:])
gaussian = np.exp(-(centers-mean_x)**2/(2*var_x))/np.sqrt(2*np.pi*var_x)
np.savetxt('px_distribution.csv', np.column_stack([centers, hist, gaussian]),
           delimiter=',', header='x,p_x,gaussian', fmt='%.8e', comments='')

fig, ax = plt.subplots(figsize=(8,6))
ax.plot(centers, hist, 'b-', lw=1.5, label=f'p(x), Lz={{Lz}}')
ax.plot(centers, gaussian, 'r--', lw=1.5, label=f'Gaussian')
ax.set_xlabel('x'); ax.set_ylabel('p(x)')
ax.set_title(f'1D Anderson: Lz={{Lz}}, W={{W}}, E={{E}}')
ax.legend(); plt.tight_layout()
plt.savefig('px_distribution.png', dpi=150)
print('Done.')
"
""",
    )
