# Reproduction Instructions for AI Agent

This document provides complete instructions for an AI agent to reproduce all computational results from the referenced physics/computational science paper.

---

## 1. Article Information

| Field | Value |
|-------|-------|
| **Title** | Parallel Implementations of the Split-Step Fourier Method for Solving Nonlinear Schrodinger Systems |
| **Authors** | S.M. Zoldi, V. Ruban, A. Zenchuk, S. Burtsev |
| **Year** | 1997 |
| **arXiv** | physics/9711012v1 |
| **Source Markdown** | `zoldi1997.md` |

---

## 2. Main Formulas Used in Reproduction

### 2.1 Nonlinear Schrodinger Equation (NLSE)
$$
i A_t + \sigma \frac{d}{2} A_{xx} + |A|^2 A = G
$$
where $\sigma = \pm 1$ (anomalous/normal dispersion), $d$ is the normalized dispersion, and $G$ is a perturbation (set to 0 for integrable case).

### 2.2 Split-Step Fourier Method (Section 2)
The solution over a short time interval $\tau$:
$$
A(t+\tau, x) = \exp(\tau L) \exp(\tau N) A(t, x)
$$
where:
- **Nonlinear operator**: $\exp(\tau N) A = A \exp(i |A|^2 \tau)$
- **Linear operator** (in Fourier space): $\exp(\tau L) B = \mathcal{F}^{-1}[\exp(-i \sigma d k^2 \tau / 2) \mathcal{F}[B]]$

The four steps per time step:
1. Nonlinear step: $A_1 = \exp(\tau N) A(t, x)$
2. Forward FFT: $A_2 = \mathcal{F}[A_1]$
3. Linear step: $A_3 = \exp(-i \sigma d k^2 \tau / 2) A_2$
4. Backward FFT: $A(t+\tau) = \mathcal{F}^{-1}[A_3]$

### 2.3 Parallel FFT via 2D Matrix Decomposition (Eq. 2-4)
For $N = M_0 \times M_1$, the 1D array $A$ is written as a 2D matrix $a_{jk}$ of size $M_0 \times M_1$.

The FFT is decomposed as:
$$
F(M_1 k_1 + k_0) = \sum_{n_0=0}^{M_0-1} f(k_0, n_0) \exp\left(\frac{-2\pi i}{N} n_0 k_0\right) \exp\left(\frac{-2\pi i}{M_0} n_0 k_1\right)
$$

where:
$$
f(k_0, n_0) = \sum_{n_1=0}^{M_1-1} A(M_0 n_1 + n_0) \exp\left(\frac{-2\pi i}{M_1} n_1 k_0\right)
$$

This decomposes into three independent steps:
1. $M_0$ independent $M_1$-size FFTs on rows
2. Multiply by twiddle factors $E_{jk} = \exp(-2\pi i \cdot j \cdot k / N)$
3. $M_1$ independent $M_0$-size FFTs on columns

### 2.4 Complete Parallel SSF Algorithm (8 steps)
1. Nonlinear step
2. Row-FFT
3. Multiply by $E$
4. Column-FFT
5. Linear step (transposed linear operator)
6. Column-BFT (backward)
7. Multiply by $E^*$
8. Row-BFT

**Key insight**: The transposition stage normally required after the 2D FFT decomposition is eliminated in SSF because the linear operator can be applied in the transposed order.

### 2.5 Speedup Formula (Eq. 5-6)
$$
SU = \frac{\tau_M N \log(N)}{\tau_M N \log(\sqrt{N})/P + \tau_c N/P + f(N/P)^2}
$$

Simplified for $N = 2^{2K}$:
$$
SU = \frac{2P}{1 + \xi/K + f \cdot 2^K/(PK)}
$$

where $\xi = \tau_c / \tau_M$ (communication-to-computation ratio) and $f$ models memory contention.

### 2.6 Bright Soliton Solution
For $\sigma = 1$, $d = 1$, $G = 0$, the exact solution:
$$
A(x, t) = a \cdot \text{sech}(a(x - vt)) \cdot \exp\left(i\left(vx + \frac{a^2 - v^2/2}{2}t\right)\right)
$$

This is used to validate the SSF implementation.

---

## 3. Main Methodology

### 3.1 Serial SSF Implementation
1. Discretize spatial domain: $x_l = -L + l \cdot h$, $l = 0, \ldots, N-1$, $h = 2L/N$
2. Compute spatial frequencies: $k_j = 2\pi j / (Nh)$ for $j = 0, \ldots, N-1$
3. For each time step $\tau = T/n_{\text{steps}}$:
   - Apply nonlinear operator pointwise: $A_1(l) = A(l) \exp(i|A(l)|^2 \tau)$
   - Forward FFT: $\hat{A}_2 = \text{FFT}(A_1)$
   - Apply linear operator in frequency space: $\hat{A}_3(k) = \hat{A}_2(k) \exp(-i \sigma d k^2 \tau / 2)$
   - Inverse FFT: $A = \text{IFFT}(\hat{A}_3)$

### 3.2 Parallel FFT Implementation
1. Choose $M_0 = M_1 = \sqrt{N}$ (requires $N$ be a perfect square)
2. Reshape 1D array into $M_0 \times M_1$ matrix: $a[n_0, n_1] = A[M_0 n_1 + n_0]$
3. Perform $M_0$ independent $M_1$-size FFTs on rows
4. Multiply by twiddle factors $E[j,k] = \exp(-2\pi i j k / N)$
5. Perform $M_1$ independent $M_0$-size FFTs on columns
6. Result stored as: $F[M_1 k_1 + k_0]$ in position $[k_1, k_0]$

### 3.3 Validation Strategy
- **Soliton preservation**: Bright soliton should maintain shape over propagation time
- **FFT accuracy**: 2D decomposition should match direct FFT to machine precision (~$10^{-15}$)
- **Serial vs parallel**: Both SSF implementations should produce identical results

---

## 4. Input Parameters

### 4.1 SSF Simulation Parameters
| Parameter | Symbol | Default Value | Description |
|-----------|--------|---------------|-------------|
| Grid points | $N$ | 1024 | Number of spatial discretization points |
| Domain half-width | $L$ | 20.0 | Spatial domain is $[-L, L]$ |
| Total time | $T$ | 10.0 | Propagation time |
| Time steps | $n_{\text{steps}}$ | 5000 | Number of SSF steps |
| Dispersion sign | $\sigma$ | 1 | Anomalous dispersion |
| Dispersion coefficient | $d$ | 1.0 | Normalized dispersion |

### 4.2 Soliton Parameters
| Parameter | Symbol | Default Value | Description |
|-----------|--------|---------------|-------------|
| Amplitude | $a$ | 1.0 | Soliton amplitude |
| Velocity | $v$ | 0.0 | Soliton velocity |
| Center | $x_0$ | 0.0 | Initial center position |

### 4.3 Paper's Timing Parameters
| Array Size $N$ | $K$ ($N=2^{2K}$) | Time Steps (Shared) | Time Steps (MPI) |
|----------------|-------------------|---------------------|------------------|
| $2^{12} = 4096$ | 6 | 8000 | 8000 |
| $2^{14} = 16384$ | 7 | 2000 | 2000 |
| $2^{16} = 65536$ | 8 | 500 | 500 |
| $2^{18} = 262144$ | 9 | 125 | 125 |

### 4.4 Parameter Variations per Figure
| Figure | Varied Parameter | Values |
|--------|-----------------|--------|
| Fig 1 | (none) | Default soliton, N=1024, T=10 |
| Fig 2 | Array size N | 16, 64, 256, 1024, 4096 |
| Fig 3 | (none) | Empirical data from paper's Tables 1 & 2 |
| Fig 4 | (none) | N=256, serial vs parallel comparison |

---

## 5. Banned Libraries

| Banned | Reason |
|--------|--------|
| Pre-built NLSE solvers (`dedalus`, `py-pde`, etc.) | Must implement SSF from scratch |
| Parallel FFT libraries (`mpi4py`, `pyfftw` with threads) | Must implement 2D decomposition manually |
| GPU-accelerated libraries (cupy, pytorch for simulation) | Keep in pure NumPy |

**Allowed**: `numpy`, `scipy`, `matplotlib`

---

## 6. Data Files to Reproduce

### 6.1 Figure 1: Soliton Propagation
**File**: `data/fig1_soliton_propagation.csv`

| Property | Value |
|----------|-------|
| Description | Amplitude of bright soliton $|A(x,t)|$ at selected spatial points over time |
| X-axis | Time $t$, range: 0 to 10 |
| Y-axis | $|A(x,t)|$, range: 0 to 1.0 |
| Parameters | $N=1024$, $L=20$, $T=10$, $n_{\text{steps}}=5000$, $a=1$, $\sigma=d=1$ |
| Data series | ~200 spatial sample columns |
| Number of time snapshots | ~100 |

**Columns**: `t, x=-10.0000, x=-9.8047, ...` (selected spatial positions)

**Additional file**: `data/fig1_soliton_2d.csv` — Full 2D snapshot data (all x points at all times) for contour plotting.

**Validation**: Initial peak = 1.0, final peak should be within 0.01% of 1.0 (soliton preservation).

---

### 6.2 Figure 2: FFT 2D Decomposition Accuracy
**File**: `data/fig2_fft_accuracy.csv`

| Property | Value |
|----------|-------|
| Description | Maximum relative error between 2D decomposed FFT and direct FFT |
| X-axis | Array size $N$ |
| Y-axis | Max relative error |
| Data series | 5 data points: N = 16, 64, 256, 1024, 4096 |

**Columns**: `N, M, max_relative_error`

**Additional file**: `data/fig2_fft_detail_N256.csv` — Detailed per-frequency comparison for N=256.

**Expected results**: All errors should be at machine precision level ($\sim 10^{-16}$ to $10^{-15}$).

---

### 6.3 Figure 3: Speedup Analysis
**File**: `data/fig3_empirical_shared.csv`

| Property | Value |
|----------|-------|
| Description | Shared memory empirical timing data from Table 1 of the paper |
| Data points | 4 array sizes: $N = 2^{12}, 2^{14}, 2^{16}, 2^{18}$ |

**Columns**: `K, N, T1_sec, T2_sec, T4_sec, SU_P2, SU_P4`

**File**: `data/fig3_empirical_mpi.csv`

| Property | Value |
|----------|-------|
| Description | Distributed memory (MPI) empirical timing data from Table 2 of the paper |

**Columns**: `K, N, T1_sec, T2_sec, T4_sec, SU_P2, SU_P4`

**File**: `data/fig3_theoretical_speedup.csv` — Theoretical curves from Eq. (6).

**Key observations from data**:
- Shared memory: Maximum speedup at $N = 2^{16}$ (SU$_4 = 3.36$), decreases for larger $N$ due to contention
- MPI: Speedup for $P=4$ increases monotonically with $N$ (SU$_4 = 3.45$ at $N = 2^{18}$)
- Both paradigms achieve near-linear speedup for small $P$ at optimal problem sizes

---

### 6.4 Figure 4: Serial vs Parallel SSF Comparison
**File**: `data/fig4_ssf_comparison.csv`

| Property | Value |
|----------|-------|
| Description | Difference between serial and parallel SSF implementations over time |
| X-axis | Time $t$, range: 0 to 5 |
| Y-axis | Max difference $|A_{\text{serial}} - A_{\text{parallel}}|$ |
| Parameters | $N=256$, $L=20$, $T=5$, 2000 steps |

**Columns**: `t, rms_difference, max_difference`

**File**: `data/fig4_fft_timing.csv` — FFT timing comparison.

**Columns**: `N, serial_time_sec, parallel_2d_time_sec, parallel_serial_ratio`

**Expected**: Serial-parallel difference should be at round-off level ($\sim 10^{-13}$).

---

## 7. Code Structure

### 7.1 Core Module
**Path**: `reproduction/ssf_core.py`

Contains:
- `serial_ssf_step()`: One step of serial SSF method
- `parallel_fft_2d()`: 2D matrix decomposed FFT (Eq. 2-4)
- `parallel_bft_2d()`: Inverse of 2D decomposed FFT
- `parallel_ssf_step()`: One step of parallel SSF (8-step algorithm)
- `soliton_initial()`: Bright soliton initial condition
- `run_ssf_simulation()`: Full SSF simulation driver
- `time_fft_methods()`: FFT timing comparison utility
- `theoretical_speedup()`: Eq. (6) speedup formula

### 7.2 Figure Scripts
| Figure | Script Path |
|--------|-------------|
| Fig 1 | `reproduction/fig1_soliton_propagation.py` |
| Fig 2 | `reproduction/fig2_fft_accuracy.py` |
| Fig 3 | `reproduction/fig3_theoretical_speedup.py` |
| Fig 4 | `reproduction/fig4_ssf_comparison.py` |

### 7.3 Output Directories
- Data files: `data/`
- Plot images: `plots/`

---

## 8. Reproduction Requirements

### 8.1 Computational Resources

| Requirement | Value | Notes |
|-------------|-------|-------|
| **Total time limit** | 30 minutes | All 4 figures |
| **Memory limit** | 1 GB | Modest requirements |
| **CPU** | Single-threaded | No parallelization needed |

### 8.2 Numerical Precision
- Export data with 8 decimal places
- FFT accuracy at machine precision (~1e-15)
- Soliton preservation within 0.01%

### 8.3 Environment
```bash
pip install numpy scipy matplotlib
```

### 8.4 Execution
Run scripts from the `reproduction/` directory:
```bash
cd reproduction/
MPLBACKEND=Agg python3 fig1_soliton_propagation.py
MPLBACKEND=Agg python3 fig2_fft_accuracy.py
MPLBACKEND=Agg python3 fig3_theoretical_speedup.py
MPLBACKEND=Agg python3 fig4_ssf_comparison.py
```
