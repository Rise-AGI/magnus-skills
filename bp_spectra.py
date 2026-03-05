from typing import Annotated, Literal

DerivativeType = Annotated[
    Literal["all", "naive", "wilson", "slac", "modified_slac"],
    {"label": "Derivative Type", "description": "Which fermion derivative to compute (or 'all' for comparison)"}
]
MuParam = Annotated[
    float,
    {"label": "mu parameter", "description": "Parameter mu for modified SLAC derivative (0 < mu < 1)"}
]
WilsonR = Annotated[
    float,
    {"label": "Wilson r", "description": "Wilson parameter r (typically 1.0)"}
]

def blueprint(
    derivative_type: DerivativeType = "all",
    mu: MuParam = 0.5,
    wilson_r: WilsonR = 1.0,
):
    submit_job(
        task_name="Melnikov2000: Energy Spectra",
        description="Compute one-particle energy spectra E_k for various lattice fermion derivatives (naive, Wilson, SLAC, modified SLAC, perfect Wilson). Reproduces Figs 1, 2, 4 of Melnikov & Weinstein (2000).",
        container_image="docker://pytorch/pytorch:2.5.1-cuda12.4-cudnn9-runtime",
        entry_command=f"pip install numpy scipy matplotlib && python3 -c \"\nimport numpy as np\nimport os\nos.environ['MPLBACKEND'] = 'Agg'\nimport matplotlib\nmatplotlib.use('Agg')\nimport matplotlib.pyplot as plt\n\ndef naive(xi): return np.sin(xi), np.zeros_like(xi)\ndef wilson(xi, r={wilson_r}): return np.sin(xi), r*(1-np.cos(xi))\ndef slac(xi): return xi.copy(), np.zeros_like(xi)\ndef mod_slac(xi, mu={mu}): return np.where(xi<mu*np.pi, xi, (mu/(mu-1))*(np.pi-xi)), np.zeros_like(xi)\ndef perfect_wilson(xi): return np.where(xi<np.pi/2, xi, (np.pi/2)*np.sin(np.pi-xi)), np.where(xi<np.pi/2, 0.0, (np.pi/2)*np.cos(np.pi-xi))\ndef E(Z,X): return np.sqrt(Z**2+X**2)\n\nxi = np.linspace(0, np.pi, 500)\nfig, ax = plt.subplots(figsize=(10,7))\nfor name, func in [('Naive',naive),('Wilson',wilson),('SLAC',slac),('Mod SLAC',lambda x: mod_slac(x,{mu})),('Perfect Wilson',perfect_wilson)]:\n    Z,X = func(xi)\n    ax.plot(xi/np.pi, E(Z,X), lw=2, label=name)\nax.set_xlabel('ka/pi'); ax.set_ylabel('E_k (1/a)'); ax.legend(); ax.grid(alpha=0.3)\nax.set_title('One-particle energy spectra: Lattice Schwinger Model')\nfig.savefig('spectra.png', dpi=150); print('Done')\n\"",
    )
