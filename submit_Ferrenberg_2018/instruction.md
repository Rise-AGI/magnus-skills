# Reproduction Instructions for AI Agent

This document provides complete instructions for an AI agent to reproduce all figures from the referenced physics paper.

---

## 1. Article Information

| Field | Value |
|-------|-------|
| **Title** | Pushing the Limits of Monte Carlo Simulations for the 3d Ising Model |
| **Authors** | Alan M. Ferrenberg, Jiahao Xu, David P. Landau |
| **Journal** | Physical Review E **97**, 043301 (2018) |
| **DOI** | 10.1103/PhysRevE.97.043301 |
| **arXiv** | 1806.03558 |
| **Source Markdown** | `ferrenberg2018.md` |

---

## 2. Main Formulas Used in Reproduction

### 2.1 Ising Model Hamiltonian
$$
\mathcal{H} = -J \sum_{\langle i,j \rangle} \sigma_i \sigma_j
$$
where $\sigma_i = \pm 1$ are spin variables on a simple cubic lattice with periodic boundary conditions, and $\langle i,j \rangle$ denotes nearest-neighbor pairs.

The dimensionless coupling constant is $K = J / k_B T$.

### 2.2 Wolff Cluster Algorithm
Bonds are drawn to nearest neighbors of a growing cluster with probability:
$$
p = 1 - e^{-2K \delta_{\sigma_i \sigma_j}}
$$
A single seed spin is chosen randomly, a cluster is grown by adding aligned neighbors with probability $p$, then all spins in the cluster are flipped.

### 2.3 Histogram Reweighting
Given a histogram $H_0(E, M)$ collected at $K_0$, the probability distribution at arbitrary $K$ is:
$$
P_K(E, M) = \frac{H_0(E, M) \, e^{\Delta K \cdot E}}{\sum_{E,M} H_0(E, M) \, e^{\Delta K \cdot E}}
$$
where $\Delta K = K_0 - K$. Any function $f(E,M)$ can be evaluated as:
$$
\langle f(E,M) \rangle_K = \sum_{E,M} f(E,M) \, P_K(E,M)
$$

### 2.4 Thermodynamic Quantities
- Logarithmic derivative of magnetization:
$$
\frac{\partial \ln \langle |m|^i \rangle}{\partial K} = \frac{\langle |m|^i E \rangle}{\langle |m|^i \rangle} - \langle E \rangle
$$
- Magnetization cumulants ($i = 1, 2, 3$):
$$
U_{2i} = 1 - \frac{\langle |m|^{2i} \rangle}{3 \langle |m|^i \rangle^2}
$$
- Specific heat:
$$
C = K^2 L^{-d} (\langle E^2 \rangle - \langle E \rangle^2)
$$
- Finite-lattice susceptibility:
$$
\chi' = K L^d (\langle |m|^2 \rangle - \langle |m| \rangle^2)
$$

### 2.5 Finite-Size Scaling (FSS) with Corrections
The peak height of a thermodynamic quantity $X$ scales as:
$$
X_{\max} = X_0 L^{1/\nu} (1 + a_1 L^{-\omega_1} + a_2 L^{-\omega_2} + a_3 L^{-\omega_\nu})
$$
with correction exponents $\omega_1 = 0.83$, $\omega_2 = 4$, $\omega_\nu = 1.6$.

The effective critical coupling scales as:
$$
K_c(L) = K_c + A_0 L^{-1/\nu} (1 + A_1 L^{-\omega_1} + A_2 L^{-\omega_2} + A_3 L^{-\omega_\nu})
$$

### 2.6 Binder Cumulant Crossing Technique
The 4th-order magnetization cumulant $U_4$ is plotted vs $K$ for different $L$. As $L \to \infty$, $U_4 \to 0$ for $K < K_c$ and $U_4 \to 2/3$ for $K > K_c$. The crossing locations are fitted to:
$$
K_{\text{cross}}(L, b) = K_c + a_1 L^{-1/\nu - \omega_1} \left(\frac{b^{-\omega_1} - 1}{b^{1/\nu} - 1}\right) + \cdots
$$
where $b = L'/L$ is the ratio of two lattice sizes.

### 2.7 Cumulant Ansatz (Self-Consistency Check)
$$
U_4(L) = U^* (1 + c \, L^{-\omega_1})
$$
where $U^* = 0.46548(5)$ is the fixed-point value.

---

## 3. Main Methodology

### 3.1 Monte Carlo Simulation
1. Initialize the lattice (all spins up or random)
2. Thermalize using $2 \times 10^5$ Wolff cluster steps (discard)
3. Collect measurements: after each Wolff step, record $E$ (dimensionless energy) and $M$ (total magnetization)
4. All simulations performed at $K_0 = 0.221654$ (near $K_c$)
5. Use Mersenne Twister RNG (32-bit and 53-bit variants tested in the paper)

### 3.2 Histogram Reweighting Analysis
1. Construct the joint histogram $H(E, M)$ from simulation data at $K_0$
2. Use the reweighting formula (Sec. 2.3) to compute any observable at nearby $K$ values
3. Locate peaks in thermodynamic derivatives using golden-section search
4. All analysis done in quadruple (128-bit) precision in the original paper

### 3.3 FSS Analysis for $\nu$
1. For each $L_{\min}$, fit Eq. (2.5) to the peak heights of $\partial U_{2i}/\partial K$ and $\partial \ln |m|^i / \partial K$
2. Use cross-correlation jackknife method to combine estimates from different quantities
3. Best fit uses three correction exponents: $\omega_1 = 0.83$, $\omega_2 = 4$, $\omega_\nu = 1.6$
4. Final estimate: $\nu = 0.629912(86)$ from $L_{\min} = 80$--144

### 3.4 FSS Analysis for $K_c$
1. Fix $\nu = 0.629912$ in Eq. (2.5)
2. Fit $K_c(L)$ from peak locations of 11 thermodynamic quantities
3. Use cross-correlation jackknife analysis
4. Final estimate: $K_c = 0.2216546262(23)$ from FSS, $K_c = 0.221654628(2)$ from cumulant crossings

---

## 4. Input Parameters

### 4.1 Physical Parameters

| Parameter | Symbol | Value | Description |
|-----------|--------|-------|-------------|
| Simulation coupling | $K_0$ | 0.221654 | Inverse temperature for MC simulation |
| Critical coupling | $K_c$ | 0.221654626(5) | Best estimate from paper |
| Lattice sizes | $L$ | 16--1024 | Linear dimension of $L^3$ cubic lattice |
| Correction exponent 1 | $\omega_1$ | 0.83 | Leading confluent correction |
| Correction exponent 2 | $\omega_2$ | 4 | Sub-leading confluent correction |
| Correction exponent 3 | $\omega_\nu$ | 1.6 | Non-linear scaling field correction |
| Thermalization | — | $2 \times 10^5$ Wolff steps | Discarded before measurement |

### 4.2 Derived Quantities

| Quantity | Symbol | Value |
|----------|--------|-------|
| Critical exponent (correlation length) | $\nu$ | 0.629912(86) |
| Critical exponent (susceptibility) | $\gamma$ | 1.23708(33) |
| Critical exponent (magnetization) | $\beta$ | 0.32630(22) |
| Cumulant fixed point | $U^*$ | 0.46548(5) |
| Ratio $\gamma/\nu$ | — | 1.96390(45) |
| Ratio $\beta/\nu$ | — | 0.51801(35) |

### 4.3 Parameter Variations per Figure

| Figure | Key varying parameter | Range |
|--------|-----------------------|-------|
| Fig 2 | $L_{\min}$ (minimum lattice size in fit) | 16--160 |
| Fig 3 | $L_{\min}$ | 16--160 |
| Fig 4 | $L_{\min}$ and number of correction exponents | 16--160, 1--3 corrections |
| Fig 5 | $K$ (inverse temperature) and $L$ | $K \in [0.2210, 0.2225]$, $L \in \{8, 12, 16, 24, 32\}$ |
| Fig 6 | $L_{\min}$ and number of correction terms | 16--192, 1--2 corrections |
| Fig 7 | $L$ at three fixed $K$ values | $L \in \{8, 12, ..., 48\}$, $K = K_c \pm 6 \times 10^{-9}$ |

---

## 5. Banned Libraries

| Library type | Examples | Reason |
|-------------|----------|--------|
| Pre-built MC frameworks | `emcee`, `pymc`, `vegas` | Must implement Wolff algorithm from scratch |
| Lattice QCD / spin model libraries | `pySpin`, `lattice_mc` | Must implement Ising model directly |
| GPU computing libraries | `cupy`, `numba.cuda`, `pycuda` | Reproducibility on CPU only |
| Parallel computing | `mpi4py`, `dask` | Single-threaded for deterministic results |
| Pre-computed lookup tables | Any cached critical exponent databases | Results must come from computation |

---

## 6. Data Files to Reproduce

### 6.1 Figure 2: Critical Exponent $\nu$ vs $L_{\min}$
**File**: `data/fig2_nu_vs_lmin.csv`

| Property | Value |
|----------|-------|
| Description | $\nu$ estimated from FSS fits with 1, 2, and 3 correction exponents |
| X-axis | $L_{\min}$ (16 to 160) |
| Y-axis | $\nu$ |
| Parameters | $\omega_1 = 0.83$, $\omega_2 = 4$, $\omega_\nu = 1.6$ |
| Data series | 3 (one per number of correction terms) |
| Number of points | 11 per series (10 for 3-correction) |

**Columns**: `L_min, nu_1corr, nu_1corr_err, nu_2corr, nu_2corr_err, nu_3corr, nu_3corr_err`

### 6.2 Figure 3: $K_c$ vs $L_{\min}$ (One Correction, Unfixed $\omega$)
**File**: `data/fig3_kc_one_correction.csv`

| Property | Value |
|----------|-------|
| Description | $K_c$ from FSS with one correction term, $\omega$ unfixed |
| X-axis | $L_{\min}$ (16 to 160) |
| Y-axis | $K_c$ |
| Data series | 1 |
| Number of points | 11 |

**Columns**: `L_min, Kc, Kc_err`

### 6.3 Figure 4: $K_c$ vs $L_{\min}$ (Fixed Correction Exponents)
**File**: `data/fig4_kc_vs_lmin.csv`

| Property | Value |
|----------|-------|
| Description | $K_c$ with 1, 2, 3 fixed correction exponents |
| X-axis | $L_{\min}$ (16 to 160) |
| Y-axis | $K_c$ |
| Parameters | $\omega_1 = 0.83$, $\omega_2 = 4$, $\omega_\nu = 1.6$ |
| Data series | 3 |
| Number of points | 11 per series (10 for 3-correction) |

**Columns**: `L_min, Kc_1corr, Kc_1corr_err, Kc_2corr, Kc_2corr_err, Kc_3corr, Kc_3corr_err`

### 6.4 Figure 5: Binder Cumulant $U_4$ vs $K$
**File**: `data/fig5_cumulant_crossing.csv`

| Property | Value |
|----------|-------|
| Description | $U_4(K)$ from Wolff MC with histogram reweighting |
| X-axis | $K \in [0.2210, 0.2225]$ (200 points) |
| Y-axis | $U_4$ |
| Parameters | $K_0 = 0.221654$, $n_{\text{therm}} = 5000$, $n_{\text{meas}} = 50000$ |
| Data series | 5 (one per lattice size: $L = 8, 12, 16, 24, 32$) |
| Number of points | 200 per series |

**Columns**: `K, U4_L8, U4_L12, U4_L16, U4_L24, U4_L32`

### 6.5 Figure 6: $K_c$ from Cumulant Crossings
**File**: `data/fig6_kc_crossing.csv`, `data/fig6_kc_crossing_2corr.csv`

| Property | Value |
|----------|-------|
| Description | $K_c$ from Binder cumulant crossing technique |
| X-axis | $L_{\min}$ (16 to 192) |
| Y-axis | $K_c$ |
| Data series | 2 (1-correction and 2-correction fits) |
| Number of points | 12 (1-corr), 5 (2-corr) |

**Columns**: `L_min, Kc_1corr, Kc_1corr_err` / `L_min, Kc_2corr, Kc_2corr_err`

### 6.6 Figure 7: Self-Consistency Check $U_4$ vs $L$
**File**: `data/fig7_u4_vs_L.csv`

| Property | Value |
|----------|-------|
| Description | $U_4$ vs lattice size at three fixed $K$ values near $K_c$ |
| X-axis | $L \in \{8, 12, 16, 24, 32, 48\}$ |
| Y-axis | $U_4$ |
| Parameters | $K = 0.221654622, 0.221654628, 0.221654634$ |
| Data series | 3 (one per $K$ value) |
| Number of points | 6 per series |

**Columns**: `L, U4_K_0.221654622, U4_K_0.221654628, U4_K_0.221654634`

---

## 7. Code Structure

| File | Description |
|------|-------------|
| `reproduction/ising3d.py` | Core module: Ising3D class (Wolff algorithm), histogram reweighting, published table data, FSS fitting functions, physical constants |
| `reproduction/fig2_nu_scaling.py` | Fig 2: $\nu$ vs $L_{\min}$ from Table II |
| `reproduction/fig3_kc_one_correction.py` | Fig 3: $K_c$ vs $L_{\min}$ from Table IV |
| `reproduction/fig4_kc_scaling.py` | Fig 4: $K_c$ vs $L_{\min}$ from Table V |
| `reproduction/fig5_cumulant_crossing.py` | Fig 5: $U_4$ vs $K$ from MC simulation |
| `reproduction/fig6_kc_crossing.py` | Fig 6: $K_c$ from crossings from Tables VII/VIII |
| `reproduction/fig7_u4_consistency.py` | Fig 7: $U_4$ vs $L$ from MC simulation |

---

## 8. Reproduction Requirements

| Requirement | Value |
|-------------|-------|
| Total compute time | < 3 hours |
| Memory | < 4 GB |
| Floating-point precision | 8 decimal places in CSV output |
| Random seed | Deterministic (seeded RNG) |
| Dependencies | `numpy`, `scipy`, `matplotlib` only |

---

## 9. Environment

```bash
pip install numpy scipy matplotlib
```

---

## 10. Execution

```bash
cd reproduction
export MPLBACKEND=Agg
python3 fig2_nu_scaling.py
python3 fig3_kc_one_correction.py
python3 fig4_kc_scaling.py
python3 fig5_cumulant_crossing.py
python3 fig6_kc_crossing.py
python3 fig7_u4_consistency.py
```
