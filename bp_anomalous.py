from typing import Annotated

WilsonR = Annotated[
    float,
    {"label": "Wilson r", "description": "Wilson parameter r (typically 0.5, 1.0, or 2.0)"}
]

def blueprint(
    wilson_r: WilsonR = 1.0,
):
    submit_job(
        task_name="Melnikov2000: Anomalous Commutator",
        description="Compute the anomalous commutator W/(k/pi) for naive and Wilson derivatives numerically, and verify against analytical values. Eq. 83-85 of Melnikov & Weinstein (2000).",
        container_image="docker://pytorch/pytorch:2.5.1-cuda12.4-cudnn9-runtime",
        entry_command=f"pip install numpy scipy matplotlib && python3 -c \"\nimport numpy as np, os\nos.environ['MPLBACKEND']='Agg'\nimport matplotlib; matplotlib.use('Agg')\nimport matplotlib.pyplot as plt\nnpts=50000\nxi=np.linspace(0,np.pi,npts); dxi=xi[1]-xi[0]\ndef compute_W(Z,X):\n    E=np.sqrt(Z**2+X**2); E=np.where(E<1e-15,1e-15,E)\n    ct,st=Z/E,X/E\n    d2Z=np.gradient(np.gradient(Z,dxi),dxi)\n    d2X=np.gradient(np.gradient(X,dxi),dxi)\n    return np.trapezoid(d2Z*ct+d2X*st,xi)\nZ_n,X_n=np.sin(xi),np.zeros(npts)\nW_naive=compute_W(Z_n,X_n)\nprint(f'Naive: W/(k/pi) = {{W_naive:.6f}} (expected -2)')\nZ_w,X_w=np.sin(xi),{wilson_r}*(1-np.cos(xi))\nW_wilson=compute_W(Z_w,X_w)\nprint(f'Wilson r={wilson_r}: W/(k/pi) = {{W_wilson:.6f}} (expected -2)')\nmu_scan=np.linspace(0.05,0.95,50)\nW_an=[-(1+m/(1-m)) for m in mu_scan]\nfig,(a1,a2)=plt.subplots(1,2,figsize=(12,5))\na1.bar(['Naive',f'Wilson r={wilson_r}'],[W_naive,W_wilson],color='steelblue')\na1.axhline(-2,color='r',ls='--'); a1.axhline(-1,color='g',ls=':')\na1.set_ylabel('W/(k/pi)'); a1.set_title('Numerical')\na2.plot(mu_scan,W_an,'b-',lw=2)\na2.axhline(-1,color='g',ls=':',label='Continuum')\na2.set_xlabel('mu'); a2.set_ylabel('W/(k/pi)'); a2.set_title('Mod SLAC (analytical)')\na2.legend(); fig.tight_layout()\nfig.savefig('anomalous.png',dpi=150); print('Done')\n\"",
    )
