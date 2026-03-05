# Reproduction Instructions for AI Agent

This document provides complete instructions for an AI agent to reproduce all computational figures from the referenced physics paper.

---

## 1. Article Information

| Field | Value |
|-------|-------|
| **Title** | The Lattice Schwinger Model: Confinement, Anomalies, Chiral Fermions and All That |
| **Authors** | Kirill Melnikov, Marvin Weinstein |
| **Affiliation** | Stanford Linear Accelerator Center, Stanford University |
| **arXiv** | hep-lat/0006029 |
| **Preprint** | SLAC-PUB-8523 |
| **Year** | 2000 |
| **Source Markdown** | `melnikov2000.md` |

---

## 2. Main Formulas Used in Reproduction

### 2.1 General Lattice Fermion Hamiltonian (Eq. 63)
The free lattice Hamiltonian in momentum space:
$$
H_f = \int_{-\pi/a}^{\pi/a} \frac{dk}{2\pi} \psi_k^\dagger \{Z_k \sigma_z + X_k \sigma_x\} \psi_k
$$
where $Z_k$ and $X_k$ define the fermion derivative. Different choices lead to different lattice fermion formulations.

### 2.2 One-Particle Energy Spectrum (Eq. 68)
$$
E_k = \sqrt{Z_k^2 + X_k^2}
$$

### 2.3 Fermion Derivative Definitions

**Naive derivative:**
$$Z_k = \sin(ka)/a, \quad X_k = 0$$

**Wilson derivative:**
$$Z_k = \sin(ka)/a, \quad X_k = \frac{r}{a}(1 - \cos(ka))$$

**SLAC derivative:**
$$Z_k = k, \quad X_k = 0$$

**Modified SLAC derivative** (Eq. 86, parametrized by $\mu \in (0,1)$):
$$
Z_k = k\,\theta(\mu\pi/a - k) + \frac{\mu}{\mu-1}\left(\frac{\pi}{a} - k\right)\theta(k - \mu\pi/a), \quad X_k = 0
$$

**Perfect Wilson derivative** (Eq. 100):
$$
Z_k = k\,\theta(\pi/(2a) - k) + \frac{\pi}{2a}\sin(\pi - ka)\,\theta(k - \pi/(2a))
$$
$$
X_k = \frac{\pi}{2a}\cos(\pi - ka)\,\theta(k - \pi/(2a))
$$

### 2.4 Anomalous Commutator (Eq. 83)
In the continuum limit $a \to 0$:
$$
\lim_{ak\to 0} W = \frac{k}{\pi} \int_0^\pi d\xi \left(\frac{d^2 Z_\xi}{d\xi^2}\cos\theta_\xi + \frac{d^2 X_\xi}{d\xi^2}\sin\theta_\xi\right)
$$
where $\xi = ka$, $\cos\theta = Z/E$, $\sin\theta = X/E$.

**Known results:**
- Naive: $W/(k/\pi) = -2$ (Eq. 84, fermion doubling)
- Wilson: $W/(k/\pi) = -2$ (Eq. 85, $r$-independent)
- Modified SLAC: $W/(k/\pi) = -(1 + \mu/(1-\mu))$ (Eq. 88)

### 2.5 Coupled Dispersion Relation (Eq. 93)
For the modified SLAC derivative, with $c = \mu/(1-\mu)$:
$$
\mathcal{M} = \begin{pmatrix} k^2 + e^2/\pi & \sqrt{c}\, e^2/\pi \\ \sqrt{c}\, e^2/\pi & c^2 k^2 + c\, e^2/\pi \end{pmatrix}
$$
Eigenvalues $\omega^2$ of $\mathcal{M}$ give the dispersion relations.

**Small-$k$ limit** ($ck^2 \ll e^2/\pi$, Eq. 92):
$$E_1 = \sqrt{c}\, k \quad (\text{massless Goldstone}), \quad E_2 = \sqrt{(1+c)\, e^2/\pi} \quad (\text{massive})$$

**Large-$k$ limit** ($ck^2 \gg e^2/\pi$, Eq. 94):
$$E_1 = \sqrt{k^2 + e^2/\pi}, \quad E_2 = \sqrt{c^2 k^2 + c\, e^2/\pi}$$

---

## 3. Main Methodology

### 3.1 Energy Spectrum Computation
1. Define $Z_k$ and $X_k$ for a chosen fermion derivative over $\xi = ka \in [0, \pi]$
2. Compute $E_k = \sqrt{Z_k^2 + X_k^2}$
3. Plot $E_k$ in units of $1/a$ as a function of $\xi/\pi$

### 3.2 Anomalous Commutator
1. For smooth derivatives (naive, Wilson): numerically evaluate $\int_0^\pi d\xi\, [d^2Z/d\xi^2 \cos\theta + d^2X/d\xi^2 \sin\theta]$ using finite differences
2. For piecewise-linear derivatives (modified SLAC): use the analytical formula from Eq. (88)
3. The result $W/(k/\pi)$ should equal $-1$ for a correct continuum limit

### 3.3 Coupled Dispersion Relation
1. Construct the $2\times 2$ matrix $\mathcal{M}$ for given $k$, $c$, $e^2/\pi$
2. Compute eigenvalues and take square roots to get dispersion branches $E_1(k)$, $E_2(k)$
3. Compare with small-$k$ and large-$k$ asymptotic limits

### 3.4 Key Numerical Details
- Use dimensionless momentum $\xi = ka \in [0, \pi]$
- Set $a = 1$ (lattice units) for spectrum plots
- For the coupled dispersion, set $e^2/\pi = 1$ (natural mass units)
- Number of grid points: 300-500 per figure is sufficient
- Second derivatives computed via centered finite differences

---

## 4. Input Parameters

### 4.1 Default Parameters
| Parameter | Value | Description |
|-----------|-------|-------------|
| $a$ | 1 | Lattice spacing (lattice units) |
| $e^2/\pi$ | 1.0 | Schwinger mass squared (natural units) |
| Grid points | 500 | Points in $[0, \pi]$ for spectra |

### 4.2 Parameter Variations per Figure
| Figure | Varied Parameter | Values |
|--------|-----------------|--------|
| Fig 1 | Derivative type | Naive, Wilson ($r=1$), SLAC, Modified SLAC ($\mu=0.5$) |
| Fig 2 | $\mu$ | 0.3, 0.7 |
| Fig 3 | $c = \mu/(1-\mu)$ | 2, 5, 20 |
| Fig 4 | (none) | Perfect Wilson derivative |
| Fig 5 | $\mu$ (scan) | 0.05 to 0.95 (100 points); discrete: naive, Wilson $r=0.5,1,2$ |

---

## 5. Banned Libraries

| Banned | Reason |
|--------|--------|
| Lattice QCD libraries (e.g., `pyquda`, `lattice`) | Must implement from scratch |
| GPU libraries (`cupy`, `pytorch`) | Keep in pure NumPy |

**Allowed**: `numpy`, `scipy`, `matplotlib`

---

## 6. Data Files to Reproduce

### 6.1 Figure 1: Derivative Comparison
**File**: `data/fig1_derivative_comparison.csv`

| Property | Value |
|----------|-------|
| Description | One-particle energy $E_k$ for 4 fermion derivatives |
| X-axis | $\xi = ka$ (range: $0$ to $\pi$) |
| Y-axis | $E_k$ (units of $1/a$) |
| Data series | 4 curves: Naive, Wilson ($r=1$), SLAC, Modified SLAC ($\mu=0.5$) |
| Number of points | 500 per series |

**Columns**: `xi, E_naive, E_wilson_r1, E_slac, E_mod_slac_0.5`

---

### 6.2 Figure 2: Modified SLAC Energy Spectrum
**File**: `data/fig2_modified_slac_spectrum.csv`

| Property | Value |
|----------|-------|
| Description | One-particle energy $E_k$ for modified SLAC at different $\mu$ |
| X-axis | $\xi = ka$ (range: $0$ to $\pi$) |
| Y-axis | $E_k$ (units of $1/a$) |
| Data series | 2 curves: $\mu = 0.3$, $\mu = 0.7$ |
| Number of points | 500 per series |

**Columns**: `xi, E_mu_0.3, E_mu_0.7`

---

### 6.3 Figure 3: Coupled Dispersion Relation
**File**: `data/fig3_coupled_dispersion.csv`

| Property | Value |
|----------|-------|
| Description | Eigenvalues of coupled system $\mathcal{M}$ (Eq. 93) |
| X-axis | $k$ (range: 0.001 to 3.0) |
| Y-axis | $E = \sqrt{\text{eigenvalue}}$ |
| Data series | 2 branches per $c$ value, 3 $c$ values |
| Number of points | 300 per series |

**Columns**: `k, E1_c_2.0, E2_c_2.0, E1_c_5.0, E2_c_5.0, E1_c_20.0, E2_c_20.0`

---

### 6.4 Figure 4: Perfect Wilson Spectrum
**File**: `data/fig4_perfect_wilson_spectrum.csv`

| Property | Value |
|----------|-------|
| Description | One-particle energy and components for perfect Wilson derivative |
| X-axis | $\xi = ka$ (range: $0$ to $\pi$) |
| Y-axis | $E_k$, $Z_k$, $X_k$ (units of $1/a$) |
| Number of points | 500 |

**Columns**: `xi, Z, X, E`

---

### 6.5 Figure 5: Anomalous Commutator
**File**: `data/fig5_anomalous_commutator.csv`

| Property | Value |
|----------|-------|
| Description | Analytical anomalous commutator $W/(k/\pi)$ vs $\mu$ for modified SLAC |
| X-axis | $\mu$ (range: 0.05 to 0.95) |
| Y-axis | $W/(k/\pi) = -(1 + \mu/(1-\mu))$ |
| Number of points | 100 |

**Columns**: `mu, W_analytical`

**Additional file**: `data/fig5_discrete_results.csv` — numerical results for smooth derivatives (naive, Wilson)

---

## 7. Code Structure

### 7.1 Core Module
**Path**: `reproduction/lattice_schwinger.py`

Contains:
- Fermion derivative functions: `naive_derivative()`, `wilson_derivative()`, `slac_derivative()`, `modified_slac_derivative()`, `perfect_wilson_derivative()`
- Energy spectrum: `energy_spectrum()`
- Anomalous commutator: `anomalous_commutator_continuum_limit()`, `anomalous_commutator_analytical()`
- Coupled dispersion: `coupled_dispersion_eigenvalues()`, `uncoupled_dispersion()`, `small_k_dispersion()`

### 7.2 Figure Scripts
| Figure | Script Path |
|--------|-------------|
| Fig 1 | `reproduction/fig1_derivative_comparison.py` |
| Fig 2 | `reproduction/fig2_modified_slac_spectrum.py` |
| Fig 3 | `reproduction/fig3_coupled_dispersion.py` |
| Fig 4 | `reproduction/fig4_perfect_wilson_spectrum.py` |
| Fig 5 | `reproduction/fig5_anomalous_commutator.py` |

### 7.3 Output Directories
- Data files: `data/`
- Plot images: `plots/`

---

## 8. Reproduction Requirements

### 8.1 Computational Resources

| Requirement | Value | Notes |
|-------------|-------|-------|
| **Total time limit** | 5 minutes | All 5 figures (purely analytical, no Monte Carlo) |
| **Memory limit** | 500 MB | Modest memory requirements |
| **CPU** | Single-threaded | No parallelization needed |

### 8.2 Numerical Precision
- Export data with 8 decimal places
- Use at least 500 grid points for smooth spectra
- Use at least 50000 points for numerical anomalous commutator integration

### 8.3 Environment
```bash
pip install numpy scipy matplotlib
```

### 8.4 Execution
Run scripts from the `reproduction/` directory:
```bash
cd reproduction/
python3 fig1_derivative_comparison.py
python3 fig2_modified_slac_spectrum.py
python3 fig3_coupled_dispersion.py
python3 fig4_perfect_wilson_spectrum.py
python3 fig5_anomalous_commutator.py
```
