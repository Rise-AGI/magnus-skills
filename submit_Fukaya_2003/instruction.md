# Reproduction Instructions for AI Agent

This document provides complete instructions for an AI agent to reproduce all figures from the referenced physics paper.

---

## 1. Article Information

| Field | Value |
|-------|-------|
| **Title** | Lattice Study of the Massive Schwinger Model with theta Term under Luscher's "Admissibility" Condition |
| **Authors** | Hidenori Fukaya, Tetsuya Onogi |
| **Journal** | Physical Review D **68**, 074503 (2003) |
| **arXiv** | hep-lat/0305004 |
| **Source Markdown** | `fukaya2003.md` |

---

## 2. Main Formulas Used in Reproduction

### 2.1 Luscher's Gauge Action (Eq. 3)
$$
S_G = \beta \sum_{x,\mu>\nu} \frac{(1 - \text{Re}\,P_{\mu\nu}(x))}{1 - \|1 - P_{\mu\nu}(x)\|/\epsilon} \quad \text{if } \|1-P_{\mu\nu}(x)\| < \epsilon
$$
where $P_{\mu\nu}(x) = U_\mu(x) U_\nu(x+\hat\mu) U_\mu^\dagger(x+\hat\nu) U_\nu^\dagger(x)$ is the plaquette.

For U(1) gauge theory in 2D: $P(x) = e^{i\phi(x)}$ where $\phi(x)$ is the plaquette angle.
- $1 - \text{Re}\,P = 1 - \cos\phi = 2\sin^2(\phi/2)$
- $|1 - P| = 2|\sin(\phi/2)|$

The action is infinite if any plaquette violates the admissibility condition $\|1-P\| \geq \epsilon$.

### 2.2 Topological Charge (Eq. 17)
$$
N = -\frac{i}{2\pi} \sum_x \ln P_{12}(x)
$$
For U(1): $N = \frac{1}{2\pi} \sum_x \text{arg}(P(x))$, where the argument is taken in $(-\pi, \pi]$.

### 2.3 Classical Gauge Configuration (Eq. 19)
The classical config minimizing the action in sector $N$ with moduli $\nu_1, \nu_2$:
$$
U_1^{cl[N]}(x,y) = \exp\left(\frac{2\pi i \nu_1}{L}\right), \quad x \neq L-1
$$
$$
U_1^{cl[N]}(L-1,y) = \exp\left(\frac{2\pi i \nu_1}{L} - \frac{2\pi i N y}{L}\right)
$$
$$
U_2^{cl[N]}(x,y) = \exp\left(\frac{2\pi i \nu_2}{L} + \frac{2\pi i N x}{L^2}\right)
$$
All plaquettes have angle $\phi_0 = 2\pi N / L^2$, giving constant field strength.

### 2.4 Minimum Action per Sector
$$
S_{G,\min}^N = \beta L^2 \cdot \frac{1 - \cos(2\pi N/L^2)}{1 - 2|\sin(\pi N/L^2)|/\epsilon}
$$
For small $N/L^2$: $S_{G,\min}^N \approx \beta \cdot 2\pi^2 N^2 / L^2$ (quadratic in $N$).

### 2.5 Domain-Wall Fermion Operator (Eqs. 14-15)
The DWF operator $D_{DW}$ acts on vectors in $\mathbb{C}^{L \times L \times L_3 \times 2}$:
$$
D_{DW}(x,s; x',s') = \frac{1}{2}\sum_{\mu=1}^{2} [(1+\gamma_\mu) U_\mu(x) \delta_{x+\hat\mu,x'} + (1-\gamma_\mu) U_\mu^\dagger(x-\hat\mu) \delta_{x-\hat\mu,x'}] \delta_{s,s'}
+ (M-3)\delta_{x,x'}\delta_{s,s'} + P_+ \delta_{s+1,s'} + P_- \delta_{s-1,s'}
+ (m-1) P_+ \delta_{s,L_3}\delta_{s',1} + (m-1) P_- \delta_{s,1}\delta_{s',L_3}
$$
The Pauli-Villars operator $D_{AP}$ is the same but with boundary terms $-2 P_\pm$ instead of $(m-1) P_\pm$.

Gamma matrices: $\gamma_1 = \begin{pmatrix}0&1\\1&0\end{pmatrix}$, $\gamma_2 = \begin{pmatrix}0&-i\\i&0\end{pmatrix}$, $\gamma_3 = \begin{pmatrix}1&0\\0&-1\end{pmatrix}$.

### 2.6 Reweighting Factor (Eq. 29)
$$
R^N(\beta,m) = \exp(-\beta S_{G,\min}^N) \cdot Det^N \cdot \exp\left[\int_\beta^\infty d\beta' \, S_{\text{subtr}}^N(\beta',m)\right]
$$
where:
- $Det^N = \frac{\int d\nu_1 d\nu_2 \, \det(D_{DW}^N)^2 / \det(D_{AP}^N)^2}{\int d\nu_1 d\nu_2 \, \det(D_{DW}^0)^2 / \det(D_{AP}^0)^2}$
- $S_{\text{subtr}}^N(\beta',m) = \langle S_G - S_{G,\min}^N \rangle_{\beta',m}^N - \langle S_G \rangle_{\beta',m}^0$

### 2.7 Theta-Vacuum Correlator (Eq. 33)
$$
C_\pi(x;\theta) = \frac{\sum_{N=-N_{\max}}^{N_{\max}} e^{iN\theta} \, C_\pi^N(x) \, R^N(\beta,m)}{\sum_{N=-N_{\max}}^{N_{\max}} e^{iN\theta} \, R^N(\beta,m)}
$$

### 2.8 Pion Mass Dependence (Eq. 35)
$$
m_\pi(m) = a \cdot m^{2/3} + b
$$
From the continuum bosonization, and with theta dependence:
$$
m_\pi(\theta) \propto |\cos(\theta/2)|^{2/3}
$$

---

## 3. Main Methodology

### 3.1 Monte Carlo with Luscher's Action
1. Initialize gauge config from classical solution for target sector $N$
2. Use Metropolis algorithm: propose link angle change $\theta \to \theta + \delta\theta$
3. Accept/reject based on $\Delta S$ including admissibility check
4. The admissibility condition automatically prevents topology changes

### 3.2 Fermion Determinant on Classical Backgrounds
1. Build DWF operator $D_{DW}$ as sparse matrix (dimension $L \times L \times L_3 \times 2$)
2. Compute $\log\det$ using `slogdet` to avoid overflow
3. Integrate over moduli $(\nu_1, \nu_2)$ using weighted sum of $n_\nu \times n_\nu$ grid points

### 3.3 Reweighting Factor Computation
1. Compute $S_{G,\min}^N$ analytically from classical configuration
2. Compute $Det^N$ from fermion determinant on classical backgrounds
3. Measure $S_{\text{subtr}}^N(\beta',m)$ by Monte Carlo at $\beta' = 0.5, 1.0, 1.5, 2.0$
4. Integrate $S_{\text{subtr}}$ using trapezoidal rule

### 3.4 Meson Correlator Measurement
1. Generate gauge configs in each sector using Metropolis MC
2. Build DWF operator for each config, solve for fermion propagator
3. Contract propagators to form pion correlator $C_\pi(x)$
4. Fit $C_\pi(x)$ to $A \cosh(m_\pi(x - L/2))$ to extract mass

---

## 4. Input Parameters

### 4.1 Physical Parameters
| Parameter | Symbol | Default Value | Description |
|-----------|--------|---------------|-------------|
| Inverse coupling | $\beta$ | 0.5 | $\beta = 1/g^2$ |
| Admissibility | $\epsilon$ | 1.0 | Luscher condition parameter |
| Lattice size | $L$ | 16 | $L \times L$ 2D lattice |
| DWF extent | $L_3$ | 6 | Domain-wall extra dimension |
| DW mass | $M$ | 0.9 | Domain-wall mass ($0 < M < 1$) |
| Fermion masses | $m$ | 0.1, 0.2, 0.3, 0.4 | Physical fermion mass |
| HMC steps | $n_{MD}$ | 50 | Molecular dynamics steps per trajectory |
| Step size | $\Delta\tau$ | 0.02 | MD integration step size |

### 4.2 Derived Quantities
| Quantity | Value |
|----------|-------|
| Plaquette angle (N=1) | $2\pi/256 \approx 0.02454$ |
| $S_{G,\min}^1$ | $\approx 0.0395$ |
| $\|1-P\|$ for N=1 | $\approx 0.0245$ (well below $\epsilon=1.0$) |
| DWF matrix dimension | $16 \times 16 \times 6 \times 2 = 3072$ |
| Max sector | $|N| \leq 4$ (higher sectors contribute $< 1.2\%$) |

### 4.3 Parameter Variations per Figure
| Figure | Varied Parameter | Values |
|--------|-----------------|--------|
| Fig 1 | Action type | Wilson ($\beta=2.0$) vs Luscher ($\beta=0.5$, $\epsilon=1.0$) |
| Fig 2 | $S_{G,\min}^N$ | $N = 0, 1, 2, \ldots, 10$ |
| Fig 3 | Fermion mass | $m = 0.1, 0.2$ for $Det^N$ |
| Fig 4 | $\beta'$ values | $0.5, 1.0, 1.5, 2.0$ for $S_{\text{subtr}}$ integration |
| Fig 5 | Fermion mass | $m = 0.1, 0.2, 0.3, 0.4$ for pion mass |
| Fig 6 | Theta | $\theta/(2\pi) = 0$ to $1.0$ in 21 steps |

---

## 5. Banned Libraries

| Banned | Reason |
|--------|--------|
| Pre-built lattice QCD libraries (`pylqcd`, `lattice-qcd`) | Must implement U(1) gauge theory from scratch |
| GPU-accelerated libraries (cupy, pytorch for simulation) | Keep in pure NumPy/SciPy |
| External fermion operator libraries | Implement DWF operator directly |
| External Monte Carlo frameworks | Implement Metropolis/HMC directly |

**Allowed**: `numpy`, `scipy` (sparse, linalg, optimize), `matplotlib`

---

## 6. Data Files to Reproduce

### 6.1 Figure 1: Topology Evolution
**File**: `data/fig1_topology_evolution.csv`

| Property | Value |
|----------|-------|
| Description | MC evolution of topological charge: Wilson vs Luscher action |
| X-axis | MC sweep (0 to 500) |
| Y-axis | Topological charge $Q$ (integer) |
| Data series | 3 columns: Wilson($\beta=2.0$), Luscher($N=0$), Luscher($N=2$) |
| Key result | Wilson: $Q$ fluctuates; Luscher: $Q$ stays fixed |

**Columns**: `sweep, Q_wilson, Q_luscher_N0, Q_luscher_N2`

---

### 6.2 Figure 2: Classical Action
**File**: `data/fig2_classical_action.csv`

| Property | Value |
|----------|-------|
| Description | Minimum Luscher action $S_{G,\min}^N$ vs $|N|^2$ |
| X-axis | $|N|^2$ (0 to 100) |
| Y-axis | $S_{G,\min}^N$ (0 to 5) |
| Number of points | 11 ($N = 0$ to $10$) |
| Key result | $S_{G,\min}^N$ is quadratic in $N$ |

**Columns**: `N, N_squared, S_Gmin, S_fit_quadratic`

---

### 6.3 Figure 3: Fermion Determinant Ratio
**File**: `data/fig3_fermion_determinant.csv`

| Property | Value |
|----------|-------|
| Description | $Det^N$ vs $|N|$ for $m=0.1$ and $m=0.2$ |
| X-axis | $|N|$ (0 to 4) |
| Y-axis | $Det^N$ (0 to 1) |
| Data series | 2 curves: $m=0.1$, $m=0.2$ |
| Key result | $Det^N$ decreases with $|N|$ due to near-zero modes |

**Columns**: `N, DetN_m0.1, DetN_m0.2`

---

### 6.4 Figure 4: Reweighting Factor
**File**: `data/fig4_reweighting_factor.csv`

| Property | Value |
|----------|-------|
| Description | Total reweighting factor $R^N(\beta=0.5, m=0.2)$ vs $N$ |
| X-axis | $N$ (0 to 4) |
| Y-axis | $R^N$ (normalized, 0 to 1) |
| Key result | Higher sectors are exponentially suppressed |

**Columns**: `N, S_Gmin_N, Det_N, integral_S_subtr, R_N`

---

### 6.5 Figure 5: Pion Mass
**File**: `data/fig5_pion_propagator.csv`, `data/fig5_pion_mass_vs_m.csv`

| Property | Value |
|----------|-------|
| Description | Pion propagator per sector and mass vs fermion mass |
| Propagator X-axis | Lattice distance $x$ (0 to 15) |
| Mass X-axis | Fermion mass $m$ (0.1, 0.2, 0.3, 0.4) |
| Key result | $m_\pi \propto m^{2/3}$ consistent with continuum theory |

**Propagator columns**: `x, C_pi_N0, C_pi_N1, C_pi_N2, C_pi_N3`
**Mass columns**: `m_fermion, m_pion, m_pion_err`

---

### 6.6 Figure 6: Theta Dependence
**File**: `data/fig6_theta_dependence.csv`

| Property | Value |
|----------|-------|
| Description | Pion mass vs $\theta$ at $m=0.2$ |
| X-axis | $\theta/(2\pi)$ (0 to 1) |
| Y-axis | $m_\pi(\theta)$ |
| Data series | Lattice data + continuum prediction $|\cos(\theta/2)|^{2/3}$ |
| Key result | Perfect agreement with continuum for $\theta/(2\pi) < 0.5$ |

**Columns**: `theta_over_2pi, theta_rad, m_pion, m_pion_err`

---

## 7. Code Structure

### 7.1 Core Module
**Path**: `reproduction/schwinger_luscher.py`

Contains:
- U(1) lattice operations: `make_classical_config()`, `plaquette_angle()`, `topological_charge()`
- Gauge actions: `luscher_action()`, `wilson_action()`, `luscher_action_min()`
- Domain-wall fermion: `build_dwf_operator()`, `fermion_determinant_ratio()`, `compute_DetN()`
- Monte Carlo: `metropolis_sweep_luscher()`, `metropolis_sweep_wilson()`
- HMC: `hmc_trajectory()` with pseudofermion method
- Meson correlators: `measure_pion_correlator()`, `measure_eta_correlator()`, `fit_correlator_mass()`

### 7.2 Figure Scripts
| Figure | Script Path |
|--------|-------------|
| Fig 1 | `reproduction/fig1_topology_evolution.py` |
| Fig 2 | `reproduction/fig2_classical_action.py` |
| Fig 3 | `reproduction/fig3_fermion_determinant.py` |
| Fig 4 | `reproduction/fig4_reweighting_factor.py` |
| Fig 5 | `reproduction/fig5_pion_mass.py` |
| Fig 6 | `reproduction/fig6_theta_dependence.py` |

### 7.3 Output Directories
- Data files: `data/`
- Plot images: `plots/`

---

## 8. Reproduction Requirements

### 8.1 Computational Resources

| Requirement | Value | Notes |
|-------------|-------|-------|
| **Total time limit** | 3 hours | Figs 4-6 dominate (MC + DWF inversions) |
| **Memory limit** | 4 GB | DWF matrix is 3072x3072 complex |
| **CPU** | Single-threaded | No parallelization required |

### 8.2 Numerical Precision
- Export data with 8 decimal places
- Use `np.linalg.slogdet` for determinant to avoid overflow
- CG solver tolerance: $10^{-10}$

### 8.3 Environment
```bash
pip install numpy scipy matplotlib
```

### 8.4 Execution
Run scripts from the `reproduction/` directory:
```bash
cd reproduction/
python3 fig1_topology_evolution.py
python3 fig2_classical_action.py
python3 fig3_fermion_determinant.py
python3 fig4_reweighting_factor.py
python3 fig5_pion_mass.py
python3 fig6_theta_dependence.py
```
