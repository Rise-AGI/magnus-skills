# Reproduction Instructions for AI Agent

This document provides complete instructions for an AI agent to reproduce key computational results from the referenced physics paper.

---

## 1. Article Information

| Field | Value |
|-------|-------|
| **Title** | Quantum simulation of the Schwinger model: A study of feasibility |
| **Authors** | Stefan Kuhn, J. Ignacio Cirac, Mari-Carmen Banuls |
| **Journal** | Physical Review A **90**, 042305 (2014) |
| **arXiv** | 1407.4995v2 |
| **DOI** | 10.1103/PhysRevA.90.042305 |
| **Source Markdown** | `kuhn2014.md` |

---

## 2. Main Formulas Used in Reproduction

### 2.1 Kogut-Susskind Lattice Hamiltonian (Eq. 1)
$$
\mathcal{H} = \frac{g^2 a}{2} \sum_n (L_n^z)^2 + m \sum_n (-1)^n \phi_n^\dagger \phi_n - \frac{i}{2a} \sum_n (\phi_n^\dagger L_n^+ \phi_{n+1} - \text{h.c.})
$$
where $g$ is the coupling constant, $m$ is the fermion mass, $a$ is the lattice spacing. The fields $\phi_n$ are single-component fermionic fields, and $L_n^+, L_n^z$ are link operators.

### 2.2 Rescaled Hamiltonian (W)
$$
W = \frac{2\mathcal{H}}{ag^2} = \sum_n (L_n^z)^2 + \mu \sum_n (-1)^n \phi_n^\dagger \phi_n + x \sum_n (\phi_n^\dagger L_n^+ \phi_{n+1} + \text{h.c.})
$$
where $x = 1/(ag)^2$ and $\mu = 2m/(ag^2)$. We study the massless case $\mu = 0$.

### 2.3 Gauss Law Constraint (Eq. 2)
$$
G_n = L_n^z - L_{n-1}^z - \phi_n^\dagger \phi_n + \frac{1}{2}[1 - (-1)^n] = 0 \quad \forall n
$$

### 2.4 Truncated cQED Link Operators (Eq. 3)
$$
L_n^+ = i \frac{a_n^\dagger b_n}{\sqrt{l(l+1)}}, \quad L_n^z = \frac{1}{2}(a_n^\dagger a_n - b_n^\dagger b_n)
$$
with $l = N_0/2$ and dimension $d = N_0 + 1$. The truncation means $L_n^+|J\rangle = 0$ (no wraparound).

### 2.5 Zd Link Operators (Eq. 4)
$$
L_n^+ = \sum_{k=-J}^{J} |\phi_{k+1}\rangle\langle\phi_k|, \quad L_n^z = \sum_{k=-J}^{J} k|\phi_k\rangle\langle\phi_k|
$$
with periodic identification $|\phi_{J+1}\rangle = |\phi_{-J}\rangle$, giving dimension $d = 2J+1$.

### 2.6 Jordan-Wigner Transformation
In the spin formulation: $\phi_n^\dagger \phi_n \to (1 + \sigma_z^{(n)})/2$, and the hopping becomes:
$$
\phi_n^\dagger L_n^+ \phi_{n+1} \to \sigma_+^{(n)} L_n^+ \sigma_-^{(n+1)}
$$

### 2.7 Adiabatic Ramp
$$
x(t) = x_F \cdot (t/T)^3
$$
Cubic ramp from $x=0$ (strong coupling) to $x = x_F$ over total time $T$.

### 2.8 Noise Term
$$
H_{\text{noise}} = \lambda x(t) \sum_n (L_n^+ + L_n^-)
$$
for the Zd model, and $\lambda x(t) \sum_n (a_n^\dagger b_n + b_n^\dagger a_n)$ for the truncated cQED model.

### 2.9 Continuum Limit
The exact energy density for the massless Schwinger model is $\omega_0/g = -1/\pi \approx -0.31831$.

---

## 3. Main Methodology

### 3.1 Exact Diagonalization Approach
The paper uses MPS (Matrix Product States) for large systems (N=50-200). Our reproduction uses exact diagonalization (ED) for smaller systems (N=4-12), which is feasible because:
1. The Gauss law constraint dramatically reduces the Hilbert space
2. For open boundary conditions, once fermion occupations are fixed, link variables are determined
3. The physical Hilbert space grows only as $\binom{N}{N/2}$ (the zero-charge sector)

### 3.2 Physical State Construction
For each fermion configuration in the zero total charge sector:
1. Compute cumulative charge to determine each link's electric field
2. For cQED: reject states where any $|L_n^z| > J$
3. For Zd: take $L_n^z \mod d$ (cyclic)

### 3.3 Hamiltonian Construction
Build the Hamiltonian matrix directly in the physical subspace:
- **Diagonal elements**: electric field energy $\sum_n (L_n^z)^2$
- **Off-diagonal elements**: hopping terms connecting states that differ by one fermion hop

### 3.4 Adiabatic Evolution
Time evolution using matrix exponential with second-order decomposition:
- Split time into $n_{\text{steps}}$ intervals of size $\delta t = T/n_{\text{steps}}$
- At each step, evaluate $H(x(t_{\text{mid}}))$ and apply $e^{-iH\delta t}$

### 3.5 Noise Simulation
For noisy evolution, work in the FULL (unconstrained) Hilbert space since noise breaks Gauss law.

---

## 4. Input Parameters

### 4.1 Physical Parameters
| Parameter | Symbol | Values Used | Description |
|-----------|--------|-------------|-------------|
| System size | $N$ | 4, 6, 8, 10, 12 | Number of lattice sites (even) |
| Link dimension | $d$ | 3, 5, 7 | Dimension of link Hilbert space (odd) |
| Lattice spacing param | $x$ | 0.5 - 25 | $x = 1/(ag)^2$ |
| Final x for adiabatic | $x_F$ | 10 | Target coupling |
| Evolution time | $T$ | 10 - 100 | Total adiabatic time |
| Noise strength | $\lambda$ | $10^{-4}$ - $3\times10^{-3}$ | Relative noise level |
| Mass parameter | $\mu$ | 0 | Massless case |

### 4.2 Numerical Parameters
| Parameter | Value | Notes |
|-----------|-------|-------|
| Time steps (adiabatic) | 200 | For Figs 3-4 |
| Time steps (noisy) | 300 | For Figs 5-6 |
| Bond dimension (paper) | 50-200 | MPS; we use ED instead |

### 4.3 Parameter Variations per Figure
| Figure | Varied Parameter | Values |
|--------|-----------------|--------|
| Fig 1 | $x$, $d$, model | $x \in \{1,4,9,16,25\}$; $d \in \{3,5,7\}$; cQED and Zd |
| Fig 2 | $x$, $d$, $N$, model | $x \in [0.5, 25]$; $d \in \{3,5\}$; $N \in \{4,6\}$; cQED and Zd |
| Fig 3 | $T$, $d$, $N$ | $T \in [10, 100]$; $d \in \{3,5\}$; $N \in \{4,6\}$; cQED model |
| Fig 4 | $T$, $d$, $N$ | Same as Fig 3; Zd model |
| Fig 5 | $\lambda$, $d$ | $\lambda \in \{0, 10^{-4}, 5\times10^{-4}, 10^{-3}, 3\times10^{-3}\}$; cQED |
| Fig 6 | $\lambda$, $d$ | Same as Fig 5; Zd model |

---

## 5. Banned Libraries

| Banned | Reason |
|--------|--------|
| Pre-built quantum simulation libraries (QuTiP, Qiskit, etc.) | Must implement ED from scratch |
| Tensor network libraries (ITensor, TeNPy, etc.) | Use ED instead of MPS |
| GPU-accelerated libraries (cupy, pytorch for simulation) | Keep in pure NumPy/SciPy |

**Allowed**: `numpy`, `scipy` (sparse, linalg, special), `matplotlib`

---

## 6. Data Files to Reproduce

### 6.1 Figure 1: Energy Density vs Lattice Spacing
**File**: `data/fig1_energy_density.csv`

| Property | Value |
|----------|-------|
| Description | Ground state energy density for truncated cQED and Zd models |
| X-axis | $ga = 1/\sqrt{x}$, range: 0.2 - 1.0 |
| Y-axis | $\omega = E_0/N$ (energy density) |
| Data series | 6 curves: 2 models x 3 d-values (3, 5, 7) |
| Number of points | 5 per series |

**Columns**: `ga, omega_cqed_d3, omega_cqed_d5, omega_cqed_d7, omega_zd_d3, omega_zd_d5, omega_zd_d7`

---

### 6.2 Figure 2: Energy Gap vs x
**File**: `data/fig2_energy_gap.csv`

| Property | Value |
|----------|-------|
| Description | Gap between ground and first excited state for both models |
| X-axis | $x$, range: 0.5 - 25 |
| Y-axis | $|E_0 - E_1|$ |
| Data series | 8 curves: 2 models x 2 d-values x 2 N-values |
| Number of points | 9 per series |

**Columns**: `x, gap_cqed_d3_N4, gap_cqed_d3_N6, gap_cqed_d5_N4, gap_cqed_d5_N6, gap_zd_d3_N4, gap_zd_d3_N6, gap_zd_d5_N4, gap_zd_d5_N6`

---

### 6.3 Figure 3: Adiabatic Preparation (cQED)
**File**: `data/fig3_adiabatic_cqed.csv`

| Property | Value |
|----------|-------|
| Description | Overlap with ground state after adiabatic evolution (cQED model) |
| X-axis | $T$ (total evolution time), range: 10 - 100 |
| Y-axis | Overlap, range: 0.9 - 1.0 |
| Data series | 4 curves: 2 N-values x 2 d-values |
| Number of points | 10 per series |

**Columns**: `T, overlap_N4_d3, error_N4_d3, overlap_N4_d5, error_N4_d5, overlap_N6_d3, error_N6_d3, overlap_N6_d5, error_N6_d5`

---

### 6.4 Figure 4: Adiabatic Preparation (Zd)
**File**: `data/fig4_adiabatic_zd.csv`

| Property | Value |
|----------|-------|
| Description | Overlap with ground state after adiabatic evolution (Zd model) |
| X-axis | $T$ (total evolution time), range: 10 - 100 |
| Y-axis | Overlap, range: 0.9 - 1.0 |
| Data series | 4 curves: 2 N-values x 2 d-values |
| Number of points | 10 per series |

**Columns**: `T, overlap_N4_d3, error_N4_d3, overlap_N4_d5, error_N4_d5, overlap_N6_d3, error_N6_d3, overlap_N6_d5, error_N6_d5`

---

### 6.5 Figure 5: Noise Effects (cQED)
**File**: `data/fig5_noise_cqed.csv`

| Property | Value |
|----------|-------|
| Description | Penalty energy and overlap under noise (cQED model) |
| X-axis | $\lambda$ (noise strength) |
| Y-axis | P/N (penalty per site), Overlap |
| Data series | 2 curves (d=3, d=5) at N=4 |
| Number of points | 5 per series |

**Columns**: `lambda, overlap_N4_d3, error_N4_d3, penalty_N4_d3, overlap_N4_d5, error_N4_d5, penalty_N4_d5`

---

### 6.6 Figure 6: Noise Effects (Zd)
**File**: `data/fig6_noise_zd.csv`

| Property | Value |
|----------|-------|
| Description | Penalty energy and overlap under noise (Zd model) |
| X-axis | $\lambda$ (noise strength) |
| Y-axis | P/N (penalty per site), Overlap |
| Data series | 2 curves (d=3, d=5) at N=4 |
| Number of points | 5 per series |

**Columns**: `lambda, overlap_N4_d3, error_N4_d3, penalty_N4_d3, overlap_N4_d5, error_N4_d5, penalty_N4_d5`

---

## 7. Code Structure

### 7.1 Core Module
**Path**: `reproduction/schwinger_ed.py`

Contains:
- `get_physical_states(N, d, model)`: enumerate Gauss law satisfying states
- `build_hamiltonian_physical(N, d, x, model)`: build Hamiltonian in physical subspace
- `compute_ground_state(N, d, x, model, n_states)`: find lowest eigenvalues
- `energy_density(N, d, x, model)`: compute E0/N
- `energy_gap(N, d, x, model)`: compute gap E1-E0
- `adiabatic_evolution(N, d, xF, T, n_steps, model)`: simulate adiabatic protocol
- `noisy_adiabatic_evolution(N, d, xF, T, n_steps, lam, model)`: noisy version in full space

### 7.2 Figure Scripts
| Figure | Script Path |
|--------|-------------|
| Fig 1 | `reproduction/fig1_energy_density.py` |
| Fig 2 | `reproduction/fig2_energy_gap.py` |
| Fig 3 | `reproduction/fig3_adiabatic_cqed.py` |
| Fig 4 | `reproduction/fig4_adiabatic_zd.py` |
| Fig 5 | `reproduction/fig5_noise_cqed.py` |
| Fig 6 | `reproduction/fig6_noise_zd.py` |

### 7.3 Output Directories
- Data files: `data/`
- Plot images: `plots/`

---

## 8. Reproduction Requirements

### 8.1 Computational Resources
| Requirement | Value | Notes |
|-------------|-------|-------|
| **Total time limit** | 2 hours | All 6 figures |
| **Memory limit** | 2 GB | ED for small systems |
| **CPU** | Single-threaded | No parallelization required |

### 8.2 Numerical Precision
- Export data with 8 decimal places
- Use double precision (float64) throughout

### 8.3 Environment
```bash
pip install numpy scipy matplotlib
```

### 8.4 Execution
Run scripts from the `reproduction/` directory:
```bash
cd reproduction/
python3 fig1_energy_density.py
python3 fig2_energy_gap.py
python3 fig3_adiabatic_cqed.py
python3 fig4_adiabatic_zd.py
python3 fig5_noise_cqed.py
python3 fig6_noise_zd.py
```

### 8.5 Key Differences from Paper
Our reproduction uses exact diagonalization for small systems instead of MPS for large ones:
- **System sizes**: N=4-12 instead of N=50-200
- **Qualitative agreement**: Results show the same physics (convergence with d, gap behavior, noise scaling)
- **Quantitative differences**: Finite-size effects are larger for small N, but extrapolation trends match
