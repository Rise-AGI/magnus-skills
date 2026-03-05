# Reproduction Instructions for AI Agent

This document provides complete instructions for an AI agent to reproduce all figures from the referenced physics paper.

---

## 1. Article Information

| Field | Value |
|-------|-------|
| **Title** | A comparison between the Split Step Fourier and Finite-Difference method in analysing the soliton collision of a type of Nonlinear Schrodinger equation found in the context of optical pulses |
| **Author** | Luke Taylor |
| **Affiliation** | University of Cape Town, Department of Mathematics and Applied Mathematics |
| **Year** | 2017 |
| **DOI** | Not available (university report) |
| **Source Markdown** | `taylor2017.md` |

---

## 2. Main Formulas Used in Reproduction

### 2.1 Nonlinear Schrodinger Equation with Saturation Nonlinearity
$$
i \frac{\partial \psi}{\partial t} + \frac{1}{2} \frac{\partial^2 \psi}{\partial x^2} + \frac{|\psi|^2 \psi}{1 + S |\psi|^2} = 0
$$
where $S$ is the saturation parameter controlling the strength and sign of the nonlinearity.

### 2.2 Soliton Solution
$$
\psi(x, t) = \frac{2\sqrt{2} e^{\sqrt{2} x}}{1 + \left(\frac{3}{2} - 2S\right) e^{2\sqrt{2} x}} e^{i t + i v x}
$$
where $v$ is the velocity parameter (phase gradient).

### 2.3 Split-Step Method — Nonlinear Part
Solving $i \psi_t + \frac{|\psi|^2}{1+S|\psi|^2} \psi = 0$ analytically:
$$
\psi(x, t+\tau) = \psi(x, t) \exp\left(i \frac{|\psi|^2}{1 + S|\psi|^2} \tau\right)
$$

### 2.4 Split-Step Method — Linear Part (Fourier Space)
Writing $\psi = \sum_n \hat{\psi}_n e^{2\pi i n x / L}$, the linear PDE $i\psi_t + \frac{1}{2}\psi_{xx} = 0$ gives:
$$
\hat{\psi}_n(t+\tau) = \hat{\psi}_n(t) \exp\left(-i \frac{1}{2} k_n^2 \tau\right), \quad k_n = \frac{2\pi n}{L}
$$

### 2.5 Forward Difference (First Step of FD Method)
$$
\psi_{j,k+1} = i\tau \left(\frac{1}{2} \psi_{xx} + A_{j,k} \psi_{j,k}\right) + \psi_{j,k}
$$

### 2.6 Central Difference (Subsequent Steps of FD Method)
$$
\psi_{j,k+1} = 2i\tau \left(\frac{1}{2} \psi_{xx} + A_{j,k} \psi_{j,k}\right) + \psi_{j,k-1}
$$

where:
$$
\psi_{xx} = \frac{\psi_{j-1,k} - 2\psi_{j,k} + \psi_{j+1,k}}{h^2}, \quad A_{j,k} = \frac{|\psi_{j,k}|^2}{1 + S|\psi_{j,k}|^2}
$$

### 2.7 Stability Condition (Von Neumann Analysis)
For the Central Difference scheme applied to the linearized PDE:
$$
\tau < \frac{h^2}{2}
$$
where $h = L/N$ is the spatial step size.

### 2.8 Conservation Quantity
$$
\mathcal{N} = \int_{-\infty}^{\infty} |\psi|^2 \, dx
$$
Computed using the composite trapezoidal rule.

---

## 3. Main Methodology

### 3.1 Split-Step Fourier Method
1. Initialize: $\psi(x, 0) = \psi_1(x) + \psi_2(x)$ (superposition of two solitons)
2. For each time step $\tau$:
   a. **Nonlinear step**: Multiply $\psi$ by $\exp(i \frac{|\psi|^2}{1+S|\psi|^2} \tau)$
   b. **Fourier transform**: $\hat{\psi} = \text{FFT}(\psi)$, then shift to center
   c. **Linear step**: Multiply each Fourier mode by $\exp(-i k_n^2 \tau / 2)$
   d. **Inverse transform**: $\psi = \text{IFFT}(\text{shift}(\hat{\psi}))$
3. Store $|\psi(x, t)|$ at each time step for visualization

### 3.2 Central Finite Difference Method
1. Initialize: $\psi_{j,0} = \psi_1(x_j) + \psi_2(x_j)$
2. First step: Use Forward Difference to compute $\psi_{j,1}$
3. For subsequent steps: Use Central Difference to compute $\psi_{j,k+1}$ from $\psi_{j,k}$ and $\psi_{j,k-1}$
4. Periodic boundary conditions: $\psi_{0,k} = \psi_{N,k}$

### 3.3 Two-Soliton Initialization
Each soliton is initialized with:
- Center position offset (x1_off, x2_off)
- Velocity parameter (v1, v2 — opposite signs for collision)
- The initial wavefunction is the sum: $\psi_0 = \psi_1 + \psi_2$

### 3.4 Visualization
All figures are 3D surface plots of $|\psi(x,t)|$ vs $(x, t)$ using matplotlib's `plot_surface` with the jet colormap.

---

## 4. Input Parameters

### 4.1 Split-Step Parameters
| Parameter | Symbol | Value | Description |
|-----------|--------|-------|-------------|
| Grid points | $N$ | 512 | Number of spatial points |
| Domain length | $L$ | 64 | Spatial domain $[0, L]$ |
| Simulation time | $T$ | 1.0 | Total time |
| Time step | $\tau$ | 0.01 | Time step size |
| Soliton 1 offset | $x_{1,\text{off}}$ | 8.0 | Center of soliton 1 |
| Soliton 2 offset | $x_{2,\text{off}}$ | 18.0 | Center of soliton 2 |

### 4.2 Finite Difference Parameters
| Parameter | Symbol | Value | Description |
|-----------|--------|-------|-------------|
| Grid points | $N$ | 512 | Number of spatial points |
| Domain length | $L$ | 30 | Spatial domain $[0, L]$ (smaller for stability) |
| Simulation time | $T$ | 1.0 | Total time |
| Time step | $\tau$ | 0.001 | Time step size (10x smaller than SS) |
| Soliton 1 offset | $x_{1,\text{off}}$ | 10.0 | Center of soliton 1 |
| Soliton 2 offset | $x_{2,\text{off}}$ | 20.0 | Center of soliton 2 |

### 4.3 Parameter Variations per Figure
| Figure | Method | S | v1 | v2 | Paper Figure |
|--------|--------|---|----|----|-------------|
| Fig 1 | Split-Step | -0.1 | 20 | -20 | Fig 2 |
| Fig 2 | Finite Diff | -0.1 | 20 | -20 | Fig 1 |
| Fig 3 | Split-Step | -0.1 | 10 | -10 | Fig 3 |
| Fig 4 | Split-Step | 0.4 | 20 | -20 | Fig 6 |
| Fig 5 | Split-Step | -10 | 10 | -10 | Fig 8 |
| Fig 6 | Both | -0.1, 0, 0.4, -10, 2 | N/A | N/A | Section 2.6 |

---

## 5. Banned Libraries

| Banned | Reason |
|--------|--------|
| Pre-built PDE solvers (e.g., `py-pde`, `dedalus`, `fenics`) | Must implement Split-Step and Finite Difference from scratch |
| GPU-accelerated libraries (cupy, pytorch for simulation) | Keep implementation in pure NumPy for clarity |
| External spectral method libraries | Implement FFT-based Split-Step directly using numpy.fft |

**Allowed**: `numpy` (including `numpy.fft`), `scipy`, `matplotlib`

---

## 6. Data Files to Reproduce

### 6.1 Figure 1: Split-Step, S=-0.1, v=+/-20
**File**: `data/fig1_ss_small_neg_v20.csv`

| Property | Value |
|----------|-------|
| Description | 3D surface of soliton collision using Split-Step with small negative S |
| X-axis | x position, range: 0-64 |
| Y-axis | time t, range: 0-1 |
| Z-axis | |psi(x,t)|, range: 0-2 |
| Parameters | S=-0.1, v1=20, v2=-20, L=64, N=512, tau=0.01 |
| Grid | 51 time points x 128 spatial points (subsampled) |

**Columns**: `time, x_0.0000, x_0.5000, x_1.0000, ...` (128 x-positions)

---

### 6.2 Figure 2: Finite Difference, S=-0.1, v=+/-20
**File**: `data/fig2_fd_small_neg_v20.csv`

| Property | Value |
|----------|-------|
| Description | 3D surface of soliton collision using Finite Difference with small negative S |
| X-axis | x position, range: 0-30 |
| Y-axis | time t, range: 0-1 |
| Z-axis | |psi(x,t)|, range: 0-2 |
| Parameters | S=-0.1, v1=20, v2=-20, L=30, N=512, tau=0.001 |
| Grid | 51 time points x 128 spatial points (subsampled) |

**Columns**: `time, x_0.0000, x_0.2344, ...` (128 x-positions)

---

### 6.3 Figure 3: Split-Step, S=-0.1, v=+/-10
**File**: `data/fig3_ss_small_neg_v10.csv`

| Property | Value |
|----------|-------|
| Description | 3D surface with slower velocity collision |
| X-axis | x position, range: 0-64 |
| Y-axis | time t, range: 0-1 |
| Z-axis | |psi(x,t)|, range: 0-2 |
| Parameters | S=-0.1, v1=10, v2=-10, L=64, N=512, tau=0.01 |
| Grid | 51 time points x 128 spatial points |

**Columns**: `time, x_0.0000, x_0.5000, ...`

---

### 6.4 Figure 4: Split-Step, S=0.4, v=+/-20
**File**: `data/fig4_ss_small_pos_v20.csv`

| Property | Value |
|----------|-------|
| Description | 3D surface showing soliton height increase with positive S |
| X-axis | x position, range: 0-64 |
| Y-axis | time t, range: 0-1 |
| Z-axis | |psi(x,t)|, range: 0-3 |
| Parameters | S=0.4, v1=20, v2=-20, L=64, N=512, tau=0.01 |
| Grid | 51 time points x 128 spatial points |

**Columns**: `time, x_0.0000, x_0.5000, ...`

---

### 6.5 Figure 5: Split-Step, S=-10, v=+/-10
**File**: `data/fig5_ss_large_neg_v10.csv`

| Property | Value |
|----------|-------|
| Description | 3D surface showing chaotic dynamics with large negative S |
| X-axis | x position, range: 0-64 |
| Y-axis | time t, range: 0-1 |
| Z-axis | |psi(x,t)| |
| Parameters | S=-10, v1=10, v2=-10, L=64, N=512, tau=0.01 |
| Grid | 51 time points x 128 spatial points |

**Columns**: `time, x_0.0000, x_0.5000, ...`

---

### 6.6 Figure 6: Conservation Comparison
**File**: `data/fig6_conservation.csv`

| Property | Value |
|----------|-------|
| Description | Conservation error |Delta N| for both methods across different S values |
| X-axis | S parameter values: -10, -0.1, 0, 0.4, 2 |
| Y-axis | |Delta N| (log scale) |
| Parameters | SS: L=64, tau=0.01; FD: L=30, tau=0.001; 8 steps each |
| Data series | 2 columns: Split-Step and Finite Difference errors |

**Columns**: `S, delta_N_splitstep, delta_N_finitediff`

---

## 7. Code Structure

### 7.1 Core Module
**Path**: `reproduction/nls_solver.py`

Contains:
- `soliton()`: Exact soliton solution with velocity and time parameters
- `init_two_soliton()`: Two-soliton superposition initialization
- `init_one_soliton()`: Single soliton initialization
- `compute_norm()`: Conservation quantity via trapezoidal rule
- `splitstep_advance()`: One Split-Step Fourier time step
- `run_splitstep()`: Full Split-Step simulation with storage
- `fd_forward()`: Forward Difference step (first time step)
- `fd_central()`: Central Difference step (subsequent steps)
- `run_finitediff()`: Full Finite Difference simulation with storage

### 7.2 Figure Scripts
| Figure | Script Path |
|--------|-------------|
| Fig 1 | `reproduction/fig1_ss_small_neg_v20.py` |
| Fig 2 | `reproduction/fig2_fd_small_neg_v20.py` |
| Fig 3 | `reproduction/fig3_ss_small_neg_v10.py` |
| Fig 4 | `reproduction/fig4_ss_small_pos_v20.py` |
| Fig 5 | `reproduction/fig5_ss_large_neg_v10.py` |
| Fig 6 | `reproduction/fig6_conservation.py` |

### 7.3 Output Directories
- Data files: `data/`
- Plot images: `plots/`

---

## 8. Reproduction Requirements

### 8.1 Computational Resources

| Requirement | Value | Notes |
|-------------|-------|-------|
| **Total time limit** | 1 hour | Figs 1,3,4,5 are fast (~5s each); Fig 2 (FD) is slowest (~5-10 min) |
| **Memory limit** | 2 GB | Modest requirements |
| **CPU** | Single-threaded | No parallelization required |

### 8.2 Numerical Precision
- Export data with 8 decimal places

### 8.3 Environment
```bash
pip install numpy scipy matplotlib
```

### 8.4 Execution
Run scripts from the `reproduction/` directory:
```bash
cd reproduction/
python3 fig1_ss_small_neg_v20.py
python3 fig2_fd_small_neg_v20.py
python3 fig3_ss_small_neg_v10.py
python3 fig4_ss_small_pos_v20.py
python3 fig5_ss_large_neg_v10.py
python3 fig6_conservation.py
```
