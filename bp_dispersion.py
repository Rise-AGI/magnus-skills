from typing import Annotated, Literal

CParam = Annotated[
    float,
    {"label": "Coupling c", "description": "c = mu/(1-mu), controls speed ratio of two fermion species"}
]
E2Pi = Annotated[
    float,
    {"label": "e^2/pi", "description": "Schwinger mass squared parameter (sets mass scale)"}
]
KMax = Annotated[
    float,
    {"label": "k_max", "description": "Maximum momentum k for the plot"}
]

def blueprint(
    c: CParam = 5.0,
    e2_over_pi: E2Pi = 1.0,
    k_max: KMax = 3.0,
):
    submit_job(
        task_name="Melnikov2000: Coupled Dispersion",
        description="Compute eigenvalues of the coupled 2x2 dispersion matrix M (Eq. 93) for the modified SLAC derivative, showing Goldstone and massive modes. Reproduces Fig 3 of Melnikov & Weinstein (2000).",
        container_image="docker://pytorch/pytorch:2.5.1-cuda12.4-cudnn9-runtime",
        entry_command=f"pip install numpy matplotlib && python3 -c \"\nimport numpy as np, os\nos.environ['MPLBACKEND']='Agg'\nimport matplotlib; matplotlib.use('Agg')\nimport matplotlib.pyplot as plt\nc={c}; e2pi={e2_over_pi}; kmax={k_max}\nk = np.linspace(0.001, kmax, 300)\nE1, E2 = [], []\nfor ki in k:\n    M = np.array([[ki**2+e2pi, np.sqrt(c)*e2pi],[np.sqrt(c)*e2pi, c**2*ki**2+c*e2pi]])\n    ev = np.sort(np.linalg.eigvalsh(M))\n    E1.append(np.sqrt(max(ev[0],0))); E2.append(np.sqrt(max(ev[1],0)))\nE1,E2 = np.array(E1),np.array(E2)\nfig,ax = plt.subplots(figsize=(8,6))\nax.plot(k,E1,'b-',lw=2,label='Lower branch')\nax.plot(k,E2,'r-',lw=2,label='Upper branch')\nax.plot(k,np.sqrt(c)*k,'b--',alpha=0.4,label='sqrt(c)*k')\nax.axhline(np.sqrt((1+c)*e2pi),color='r',ls='--',alpha=0.4,label='sqrt((1+c)*e2/pi)')\nax.set_xlabel('k'); ax.set_ylabel('E'); ax.set_title(f'Coupled dispersion c={c}')\nax.legend(); ax.grid(alpha=0.3)\nfig.savefig('dispersion.png',dpi=150)\nnp.savetxt('dispersion.csv',np.column_stack([k,E1,E2]),delimiter=',',header='k,E1,E2',comments='')\nprint('Done')\n\"",
    )
