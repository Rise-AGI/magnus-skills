# Reproduction Instructions for AI Agent

This document provides complete instructions for an AI agent to reproduce all computational figures from the referenced paper.

---

## 1. Article Information

| Field | Value |
|-------|-------|
| **Title** | The Energy Conserving Particle-in-Cell Method |
| **Authors** | Stefano Markidis, Giovanni Lapenta |
| **Journal** | Journal of Computational Physics 230 (2011) 7037-7052 |
| **DOI** | 10.1016/j.jcp.2011.05.033 |
| **Source Markdown** | `markidis2011.md` |

---

## 2. Main Formulas Used in Reproduction

### 2.1 Particle Equations of Motion (Midpoint Rule)
The particles obey (CGS units):
$$
\frac{d\mathbf{x}_p}{dt} = \mathbf{v}_p, \quad
\frac{d\mathbf{v}_p}{dt} = \frac{q_s}{m_s}\left(\mathbf{E}_p + \frac{\mathbf{v}_p}{c} \times \mathbf{B}_p\right)
$$

Discretized with the implicit midpoint rule:
$$
\mathbf{v}_p^{n+1} = \mathbf{v}_p^n + \frac{q_s}{m_s}\Delta t\left(\bar{\mathbf{E}}_p + \frac{\bar{\mathbf{v}}_p}{c}\times\bar{\mathbf{B}}_p\right)
$$
$$
\mathbf{x}_p^{n+1} = \mathbf{x}_p^n + \bar{\mathbf{v}}_p \Delta t
$$
where bars denote time-averages: $\bar{q} = (q^{n+1} + q^n)/2$.

### 2.2 Cloud-in-Cell Interpolation
$$
W(\mathbf{x}_g - \mathbf{x}_p) = \begin{cases} 1 - |\mathbf{x}_g - \mathbf{x}_p|/\Delta x & \text{if } |\mathbf{x}_g - \mathbf{x}_p| < \Delta x \\ 0 & \text{otherwise} \end{cases}
$$

### 2.3 Electrostatic Field Equations
For the electrostatic limit (unmagnetized plasma):
$$
\bar{\mathbf{v}}_p = \mathbf{v}_p^n + \frac{q_s}{4m_s}(\mathbf{E}_p^{n+1}(\bar{\mathbf{x}}_p) + \mathbf{E}_p^n(\bar{\mathbf{x}}_p))\Delta t
$$
$$
\mathbf{E}_g^{n+1} - \mathbf{E}_g^n = -4\pi \bar{\mathbf{J}}_g \Delta t
$$

### 2.4 Average Current Density
$$
\bar{\mathbf{J}}_g = \sum_s \sum_p q_s \bar{\mathbf{v}}_p W(\mathbf{x}_g - \bar{\mathbf{x}}_p) / V_g
$$

### 2.5 Cold Plasma Numerical Dispersion (EC-PIC)
$$
\tan\left(\frac{\omega\Delta t}{2}\right)\sin\left(\frac{\omega\Delta t}{2}\right) = \left(\frac{\omega_{pe}\Delta t}{2}\right)^2
$$
This always has **real roots**, so the scheme is linearly unconditionally stable.

### 2.6 Cold Plasma Numerical Dispersion (Explicit PIC)
$$
4\sin^2\left(\frac{\omega\Delta t}{2}\right) = (\omega_{pe}\Delta t)^2
$$
This has exponential growth for $\omega_{pe}\Delta t > 2$.

### 2.7 Two-Stream Instability Growth Rate
For two cold counter-streaming beams at $\pm v_0$ with $k = 1 \omega_{pe}/c$:
$$
\gamma = \sqrt{3}/2^{5/3} \cdot \omega_{pe} \approx 0.35355\,\omega_{pe}
$$

### 2.8 Weibel Instability Growth Rate
For bi-Maxwellian with anisotropy $a = v_{thy}^2/v_{thx}^2 - 1 = 15$, the $k=1$ mode:
$$
\gamma \approx 0.22\,\omega_{pe}
$$

### 2.9 Stability Condition for Momentum Non-Conservation
$$
v_c \Delta t / \Delta x < 1.5
$$
where $v_c$ is the characteristic velocity (thermal or drift).

---

## 3. Main Methodology

### 3.1 Energy-Conserving PIC Algorithm (Electrostatic)
1. Initialize particles (positions, velocities) and solve Poisson equation for initial E field
2. At each time step, solve the coupled system of particle velocity equations and Ampere's law simultaneously using a Jacobian-Free Newton-Krylov (JFNK) solver
3. The unknowns are: average velocities $\bar{\mathbf{v}}_p$ for all particles and the new electric field $\mathbf{E}_g^{n+1}$ on all grid points
4. The residual function evaluates: (a) velocity residual from midpoint equation, (b) field residual from Ampere's law with current computed from average quantities
5. Update particle positions and velocities from the converged average velocity

### 3.2 Explicit Momentum-Conserving PIC (for comparison)
1. Standard leapfrog/Boris push scheme
2. Deposit charge -> Solve Poisson -> Push particles
3. Does NOT conserve energy; exhibits finite grid instability when $\Delta x > 2\lambda_D$

### 3.3 Electromagnetic Extension (for Weibel test)
1. Full Maxwell equations with Yee lattice discretization
2. Boris velocity push with half B-step / push / half B-step
3. Faraday + Ampere laws with staggered grid

### 3.4 Key Numerical Details
- **JFNK solver**: scipy.optimize.newton_krylov with relative tolerance matching paper's tolerance settings
- **Interpolation**: Cloud-in-Cell (linear, 1st order) on periodic grid
- **Poisson solver**: Spectral (FFT-based) for periodic boundaries
- **Boundary conditions**: Periodic in all cases

---

## 4. Input Parameters

### 4.1 Common Parameters
| Parameter | Symbol | Value | Description |
|-----------|--------|-------|-------------|
| Plasma frequency | $\omega_{pe}$ | 1.0 | Normalized |
| Speed of light | $c$ | 1.0 | Normalized |
| Electron charge-to-mass | $q_e/m_e$ | -1.0 | Normalized |
| Ion mass ratio | $m_i/m_e$ | 1836 | Proton mass |

### 4.2 Two-Stream Instability (Figures 5, 8, 9)
| Parameter | Value | Description |
|-----------|-------|-------------|
| Beam drift | $\pm 0.2c$ | Counter-streaming |
| Box length $L$ | $2.053\,c/\omega_{pe}$ | One wavelength of k=1 mode |
| Grid cells $N_g$ | 64 | |
| Particles $N$ | 20000 | Reduced from 200000 |
| Time step $\Delta t$ | $0.1\,\omega_{pe}^{-1}$ | |
| Cycles | 300 | |
| Perturbation | $x_0 \to x_0 + \Delta x\sin(2\pi x_0/L)$ | Seed k=1 mode |

### 4.3 Maxwellian Plasma / Finite Grid Test (Figure 7)
| Parameter | Value | Description |
|-----------|-------|-------------|
| Thermal velocity | $0.2c$ | |
| Box length | $50\pi\,c/\omega_{pe}$ | Large domain |
| Grid cells | 64 | $\Delta x \approx 12\lambda_D$ |
| Particles | 10000 | Reduced from 50000 |
| Time step | $0.5\,\omega_{pe}^{-1}$ | |
| Cycles | 200 | |

### 4.4 Weibel Instability (Figures 10, 11)
| Parameter | Value | Description |
|-----------|-------|-------------|
| $v_{thy}$ | $0.4c$ | Hot direction |
| $v_{thx}$ | $0.1c$ | Cold direction ($a=15$) |
| Box length | $2\pi\,c/\omega_{pe}$ | |
| Grid cells | 64 | |
| Particles | 50000 | Reduced from 100000 |
| Time step | $0.25\,\omega_{pe}^{-1}$ | |
| Cycles | 400 | |

---

## 5. Banned Libraries

| Banned | Reason |
|--------|--------|
| Pre-built PIC simulation codes (e.g., `EPOCH`, `SMILEI`, `PIConGPU`) | Must implement PIC from scratch |
| GPU libraries (`cupy`, `torch` for simulation) | Pure NumPy/SciPy |
| Pre-built plasma dispersion solvers | Implement directly |

**Allowed**: `numpy`, `scipy` (optimize, sparse, fft, special), `matplotlib`

---

## 6. Data Files to Reproduce

### 6.1 Cold Plasma Dispersion Relation
**File**: `data/dispersion_cold_plasma.csv`

| Property | Value |
|----------|-------|
| Description | Numerical dispersion relation for EC-PIC vs explicit PIC in cold plasma |
| X-axis | $\omega_{pe}\Delta t$, range: 0.1-8.0 |
| Y-axis | $\omega/\omega_{pe}$ |
| Series | EC-PIC (always stable), Explicit PIC (unstable for $\omega_{pe}\Delta t > 2$) |
| Points | 200 |

**Columns**: `wpe_dt, omega_over_wpe_EC, omega_over_wpe_explicit`

---

### 6.2 Figure 5: Energy Conservation vs Solver Tolerance
**File**: `data/fig5_tolerance_energy.csv`

| Property | Value |
|----------|-------|
| Description | Energy variation in two-stream instability for different JFNK tolerances |
| X-axis | Time [$\omega_{pe}^{-1}$] |
| Y-axis | $\Delta E/E_0$ [%] |
| Series | Tolerances: $10^{-5}$, $10^{-7}$, $10^{-9}$ |
| Points | 201 per series |

**Columns**: `time, dE_tol_1e-05_percent, dE_tol_1e-07_percent, dE_tol_1e-09_percent`

---

### 6.3 Figure 7: Finite Grid Instability Energy History
**File**: `data/fig7_finite_grid_energy.csv`

| Property | Value |
|----------|-------|
| Description | Energy history for Maxwellian plasma: EC-PIC (stable) vs explicit PIC (heating) |
| X-axis | Time [$\omega_{pe}^{-1}$] |
| Y-axis | Total energy, $\Delta E/E_0$ [%] |
| Series | EC-PIC, Explicit PIC |
| Points | 201 |

**Columns**: `time, energy_ec, energy_explicit, dE_ec_percent, dE_explicit_percent`

---

### 6.4 Figure 8: Two-Stream Instability Growth
**File**: `data/fig8_twostream_growth.csv`

| Property | Value |
|----------|-------|
| Description | Growth of $k=1$ spectral component vs linear theory ($\gamma=0.35355\omega_{pe}$) |
| X-axis | Time [$\omega_{pe}^{-1}$] |
| Y-axis | $|E_{k=1}|$ (log scale) |
| Series | EC-PIC simulation, linear theory |
| Points | 300 |

**Columns**: `time, E_k1_simulation, E_k1_linear_theory`

---

### 6.5 Figure 9: Two-Stream Energy Conservation
**File**: `data/fig9_twostream_energy.csv`

| Property | Value |
|----------|-------|
| Description | Total energy: EC-PIC (~$10^{-4}$% variation) vs explicit (~5% variation) |
| X-axis | Time [$\omega_{pe}^{-1}$] |
| Y-axis | Total energy, $\Delta E/E_0$ [%] |
| Series | EC-PIC, Explicit PIC |
| Points | 301 |

**Columns**: `time, energy_ec, energy_explicit, dE_ec_percent, dE_explicit_percent`

---

### 6.6 Figure 10: Weibel Instability Growth
**File**: `data/fig10_weibel_growth.csv`

| Property | Value |
|----------|-------|
| Description | Growth of $B_z$ $k=1$ component vs linear theory ($\gamma=0.22\omega_{pe}$) |
| X-axis | Time [$\omega_{pe}^{-1}$] |
| Y-axis | $|B_{z,k=1}|$ (log scale) |
| Series | PIC simulation, linear theory |
| Points | 400 |

**Columns**: `time, Bz_k1_simulation, Bz_k1_linear_theory`

---

### 6.7 Figure 11: Weibel Energy and Momentum
**File**: `data/fig11_weibel_energy_momentum.csv`

| Property | Value |
|----------|-------|
| Description | Energy conservation and momentum oscillation in Weibel instability |
| X-axis | Time [$\omega_{pe}^{-1}$] |
| Y-axis | Energy, momentum, $\Delta E/E_0$ [%] |
| Points | 400 |

**Columns**: `time, total_energy, total_momentum, dE_percent`

---

## 7. Code Structure

### 7.1 Core Modules
| Module | Description |
|--------|-------------|
| `reproduction/ecpic_es.py` | Electrostatic EC-PIC with JFNK solver (scipy.optimize.newton_krylov) |
| `reproduction/ecpic_em.py` | Electromagnetic PIC with Boris push (for Weibel test) |

### 7.2 Figure Scripts
| Figure | Script |
|--------|--------|
| Dispersion | `reproduction/dispersion_cold_plasma.py` |
| Fig 5 | `reproduction/fig5_tolerance_energy.py` |
| Fig 7 | `reproduction/fig7_finite_grid_energy.py` |
| Fig 8 | `reproduction/fig8_twostream_growth.py` |
| Fig 9 | `reproduction/fig9_twostream_energy.py` |
| Figs 10,11 | `reproduction/fig10_11_weibel.py` |

### 7.3 Output Directories
- Data: `data/`
- Plots: `plots/`

---

## 8. Reproduction Requirements

### 8.1 Computational Resources

| Requirement | Value | Notes |
|-------------|-------|-------|
| **Total time limit** | 3 hours | All figures |
| **Memory** | 2 GB | Moderate |
| **CPU** | Single-threaded | No parallelization |

### 8.2 Numerical Precision
- Export data with 8 decimal places
- JFNK solver tolerance: 1e-7 (default), varies for Figure 5

### 8.3 Environment
```bash
pip install numpy scipy matplotlib
```

### 8.4 Execution
```bash
export MPLBACKEND=Agg
cd reproduction/
python3 dispersion_cold_plasma.py
python3 fig5_tolerance_energy.py
python3 fig7_finite_grid_energy.py
python3 fig8_twostream_growth.py
python3 fig9_twostream_energy.py
python3 fig10_11_weibel.py
```

### 8.5 Key Implementation Notes

1. **JFNK Solver**: The heart of the EC-PIC method is solving the coupled particle-field system simultaneously. Use `scipy.optimize.newton_krylov` with the residual function that evaluates both the velocity equation (midpoint rule) and the field equation (Ampere's law with current from average quantities).

2. **Energy Conservation**: The midpoint rule + consistent current density definition = exact energy conservation (up to solver tolerance). The key is that the current is deposited using **average** positions and velocities, not the old or new values.

3. **Interpolation Consistency**: The same CIC interpolation must be used for both gathering fields to particles AND scattering current to the grid. Any inconsistency breaks energy conservation.

4. **Poisson Solver**: Use FFT-based solver for periodic boundaries. This gives the initial E field satisfying Gauss' law.

5. **Particle Count**: The paper uses 200000 particles for two-stream and 100000 for Weibel. We reduce these for computational feasibility. The qualitative behavior (growth rates, energy conservation) is preserved; statistical noise increases.

6. **Normalization**: All quantities are in CGS-like normalized units with $\omega_{pe} = 1$, $c = 1$, $m_e = 1$, $|e| = 1$.
