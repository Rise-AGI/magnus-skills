# Reproduction Instructions for AI Agent

This document provides complete instructions for an AI agent to reproduce all figures from the referenced physics paper.

---

## 1. Article Information

| Field | Value |
|-------|-------|
| **Title** | Destabilization of magnetosonic-whistler waves by a relativistic runaway beam |
| **Authors** | T. Fulop, G. Pokol, P. Helander, M. Lisak |
| **Journal** | Physics of Plasmas **13**, 062506 (2006) |
| **DOI** | 10.1063/1.2208327 |
| **Source Markdown** | `fulop2006.md` |

---

## 2. Main Formulas Used in Reproduction

### 2.1 Plasma Frequencies and Cyclotron Frequencies
$$
\omega_{pe} = \sqrt{\frac{n_e e^2}{m_e \varepsilon_0}}, \quad
\omega_{pi} = \sqrt{\frac{n_e e^2}{m_i \varepsilon_0}}, \quad
\omega_{ce} = \frac{eB}{m_e}, \quad
\omega_{ci} = \frac{eB}{m_i}
$$

### 2.2 Background Dispersion Relation (Magnetosonic-Whistler)
$$
\omega_0 = k v_A \sqrt{1 + \frac{k_\parallel^2 c^2}{\omega_{pi}^2}}
$$
where $v_A = B / \sqrt{\mu_0 n_e m_i}$ is the Alfven speed.

**Limiting cases:**
- Magnetosonic: $\omega = k v_A$ (low frequency, $k_\parallel^2 c^2 / \omega_{pi}^2 \ll 1$)
- Whistler: $\omega = k k_\parallel v_A^2 / \omega_{ci}$ (high frequency, $k_\parallel^2 c^2 / \omega_{pi}^2 \gg 1$)

### 2.3 Cold-Plasma Dispersion Equation (Eq. 10)
$$
(\varepsilon_{11} - k_\parallel^2 c^2 / \omega^2)(\varepsilon_{22} - k^2 c^2 / \omega^2) + \varepsilon_{12}^2 = 0
$$
where $\varepsilon_{ij}$ are Stix cold-plasma dielectric tensor elements (S, D, P parameters).

### 2.4 Growth Rate -- Simplified Analytical (Eq. 22)
$$
\gamma_i(\omega_0, k, k_\parallel) = \frac{\pi}{4 c_Z} \frac{\omega_{pr}^2}{\omega_{pi}^2} \frac{k^2 v_A^2}{\omega_0} \exp\left[-\frac{\omega_{ce}}{(k_\parallel c - \omega_0) c_Z}\right]
$$

### 2.5 Growth Rate -- Full Analytical with Bessel Functions (Eq. 20)
$$
\frac{\gamma_i}{\omega_0} = \hat{C} \int_0^\infty p_\perp \, dp_\perp \, n^2 J_n^2(K_\perp p_\perp) \exp(a_n p_\perp^2) \left[b_n + a_n(1-y) K_\parallel p_\perp^2 / n\right]
$$
For $n = -1$ (anomalous Doppler resonance), evaluated using modified Bessel functions $I_0(\lambda)$, $I_1(\lambda)$.

### 2.6 Collisional Damping Rate
$$
\gamma_d = -\frac{1.5}{\tau_{ei}}, \quad \tau_{ei} = \frac{3\pi^{3/2} m_e^2 v_{Te}^3 \varepsilon_0^2}{n_e Z^2 e^4 \ln\Lambda}
$$

**IMPORTANT**: The thermal velocity convention is $v_{Te} = \sqrt{2 T_e / m_e}$ (rms thermal speed), NOT $\sqrt{T_e / m_e}$. This factor of $2\sqrt{2}$ in $\tau_{ei}$ is critical for correct threshold values.

### 2.7 Stability Threshold (Eq. 25)
$$
\frac{n_r}{n_e} > \frac{Z^2 B_T}{20 \, T_{eV}^{3/2}}
$$
Solving for critical temperature:
$$
T_{eV} = \left(\frac{Z^2 B_T}{20 \cdot n_r/n_e}\right)^{2/3}
$$

### 2.8 Maximum Growth Rate (Eq. 24)
$$
\gamma_{i,\max} \approx 1.3 \times 10^{-9} \frac{n_r}{B_T} \quad [\text{rad/s, SI units}]
$$

---

## 3. Main Methodology

### 3.1 Analytical Approach
1. Compute background dispersion $\omega_0(k, k_\parallel)$ from Section 2.2
2. Use Eq.(22) or Eq.(20) to compute growth rate $\gamma_i$ as a function of $k$ and $\theta_k$
3. Compare $\gamma_i$ to collisional damping $|\gamma_d|$; instability when $\gamma_i / |\gamma_d| > 1$

### 3.2 Numerical Approach
1. Solve the full Stix cold-plasma dispersion (Eq. 10) using Newton's method with continuation
2. Compute growth rate using Eq.(20) with the Stix-solved $\omega$ (instead of analytical $\omega_0$)
3. Key differences from analytical: full Stix susceptibility deviates from $\omega_0$ at high $k$

### 3.3 Key Numerical Details
- **Resonance condition**: Only $n = -1$ anomalous Doppler resonance is retained; requires $k_\parallel c > \omega_0$ (i.e., $p_{\text{res}} > 0$)
- **Stix solver**: Newton iteration with 10% step-size limit; use continuation (sweep from low $k$) for robustness
- **Numerical validity range**: Only near-perpendicular propagation ($\theta_k \gtrsim 75$) satisfies $|k| \gg |k_\parallel|$ assumption
- **Overflow protection**: Guard exponent $< -500$ to prevent underflow in exponential terms

### 3.4 Parameter Scans (Figs 4-6)
For each propagation angle $\theta_k$:
1. Scan wave number $k$ from 10 to 5000 m$^{-1}$ (400 points for analytical, 100 for numerical)
2. Find the $k$ giving maximum growth rate
3. Compute normalized ratio $\gamma_i / |\gamma_d| - 1$; threshold at 0

---

## 4. Input Parameters

### 4.1 Default Physical Parameters (Pure Deuterium)
| Parameter | Symbol | Default Value | Description |
|-----------|--------|---------------|-------------|
| Electron density | $n_e$ | $5 \times 10^{19}$ m$^{-3}$ | Background plasma density |
| Magnetic field | $B_T$ | 2.0 T | Toroidal magnetic field |
| Temperature | $T_{eV}$ | 10.0 eV | Electron temperature |
| Ion mass | $m_i$ | $2 m_p$ | Deuterium |
| Charge number | $Z$ | 1 | Pure deuterium |
| Runaway fraction | $n_r / n_e$ | $5 \times 10^{-3}$ | Runaway electron fraction |
| Normalized E-field | $E$ | 50 | $E \gg 1$ (typical) |
| Coulomb logarithm | $\ln\Lambda$ | 15 | Standard tokamak value |

### 4.2 Derived Quantities (at default parameters)
| Quantity | Symbol | Value |
|----------|--------|-------|
| Ion plasma frequency | $\omega_{pi}$ | $\sim 6.58 \times 10^9$ rad/s |
| Electron cyclotron frequency | $\omega_{ce}$ | $\sim 3.52 \times 10^{11}$ rad/s |
| Alfven speed | $v_A$ | $\sim 4.36 \times 10^6$ m/s |
| Collisional damping | $|\gamma_d|$ | $\sim 1.03 \times 10^8$ rad/s |
| Thermal velocity | $v_{Te}$ | $\sim 1.88 \times 10^6$ m/s |

### 4.3 Parameter Variations per Figure
| Figure | Varied Parameter | Values |
|--------|-----------------|--------|
| Fig 1 | $n_r/n_e$ | $5 \times 10^{-4}$, $10^{-3}$, $5 \times 10^{-3}$ |
| Fig 2 | (none) | Fixed: $\theta_k = 85$deg |
| Fig 3 | (none) | Fixed: $\theta_k = 85$deg |
| Fig 4 | $B_T$ | 1.5, 2.0, 3.0 T |
| Fig 5 | $n_e$ | $2 \times 10^{19}$, $5 \times 10^{19}$, $10^{20}$ m$^{-3}$ |
| Fig 6 | $T_{eV}$ | 5, 10, 15 eV |

---

## 5. Banned Libraries

The following libraries/approaches are **NOT allowed** for the core physics computation:

| Banned | Reason |
|--------|--------|
| Pre-built plasma physics libraries (e.g., `plasmapy`, `pyplasma`) | Must implement dispersion and growth rate from scratch |
| GPU-accelerated libraries (cupy, pytorch for simulation) | Keep implementation in pure NumPy/SciPy for clarity |
| External dispersion solvers or MHD codes | Implement Stix cold-plasma solver directly |

**Allowed**: `numpy`, `scipy` (constants, special functions, integrate, optimize), `matplotlib`

---

## 6. Data Files to Reproduce

### 6.1 Figure 1: Stability Threshold
**File**: `data/fig1_threshold.csv`

| Property | Value |
|----------|-------|
| Description | Critical temperature vs magnetic field for different runaway fractions |
| X-axis | $B_T$ [T], range: 1-4 |
| Y-axis | $T_{eV}$, range: 0-50 |
| Data series | 3 curves: $n_r/n_e = 5 \times 10^{-4}$, $10^{-3}$, $5 \times 10^{-3}$ |
| Number of points | 200 per series |

**Columns**: `B_T, T_crit_nr_5e-4, T_crit_nr_1e-3, T_crit_nr_5e-3`

---

### 6.2 Figure 2: Dispersion Relation Comparison
**File**: `data/fig2_dispersion.csv`

| Property | Value |
|----------|-------|
| Description | Dispersion relations at $\theta_k = 85$deg: analytical, magnetosonic, whistler, numerical |
| X-axis | $k$ [m$^{-1}$], range: 0-900 |
| Y-axis | $\omega_0 / \omega_{ci}$, range: 0-150 |
| Parameters | Default plasma parameters |
| Data series | 3 analytical curves (300 pts) + 1 numerical series (60 pts) |

**Columns**: `k, omega_magnetosonic, omega_whistler, omega_analytic` (+ numerical appended as comments)

---

### 6.3 Figure 3: Growth Rate vs Wave Number
**File**: `data/fig3_growth_rate.csv`

| Property | Value |
|----------|-------|
| Description | Growth rate comparison at $\theta_k = 85$deg: Eq.(22), Eq.(20), numerical |
| X-axis | $k$ [m$^{-1}$], range: 0-900 |
| Y-axis | $\gamma_i / \omega_{ci}$, range: 0-1.8 |
| Parameters | Default plasma parameters |
| Data series | 2 analytical curves (200 pts) + 1 numerical series (40 pts) |

**Columns**: `k, gamma_eq22, gamma_eq20` (+ numerical appended as comments)

---

### 6.4 Figure 4: Growth Rate vs Angle (Different B)
**File**: `data/fig4_angle_B.csv`

| Property | Value |
|----------|-------|
| Description | Normalized growth rate vs propagation angle for $B = 1.5, 2, 3$ T |
| X-axis | $\theta_k$ [deg], range: 0-90 |
| Y-axis | $\gamma_i / |\gamma_d| - 1$, range: $-1.1$ to $1.0$ |
| Parameters | $n_e = 5 \times 10^{19}$, $T = 10$ eV, $n_r/n_e = 5 \times 10^{-3}$ |
| Data series | 3 analytical curves (80 pts each) + 3 numerical series (12 pts each) |

**Columns**: `theta_deg, ratio_analytic_B1.5, ratio_analytic_B2.0, ratio_analytic_B3.0` (+ numerical appended as comments)

---

### 6.5 Figure 5: Growth Rate vs Angle (Different Density)
**File**: `data/fig5_angle_density.csv`

| Property | Value |
|----------|-------|
| Description | Normalized growth rate vs angle for $n_e = 2 \times 10^{19}$, $5 \times 10^{19}$, $10^{20}$ m$^{-3}$ |
| X-axis | $\theta_k$ [deg], range: 0-90 |
| Y-axis | $\gamma_i / |\gamma_d| - 1$, range: $-1.1$ to $1.0$ |
| Parameters | $B = 2$ T, $T = 10$ eV, $n_r/n_e = 5 \times 10^{-3}$ |
| Data series | 3 analytical curves (80 pts each) + 3 numerical series (12 pts each) |

**Columns**: `theta_deg, ratio_analytic_ne_2e19, ratio_analytic_ne_5e19, ratio_analytic_ne_1e20` (+ numerical appended as comments)

---

### 6.6 Figure 6: Growth Rate vs Angle (Different Temperature)
**File**: `data/fig6_angle_temperature.csv`

| Property | Value |
|----------|-------|
| Description | Normalized growth rate vs angle for $T = 5, 10, 15$ eV |
| X-axis | $\theta_k$ [deg], range: 0-90 |
| Y-axis | $\gamma_i / |\gamma_d| - 1$ |
| Parameters | $B = 2$ T, $n_e = 5 \times 10^{19}$, $n_r/n_e = 5 \times 10^{-3}$ |
| Data series | 3 analytical curves (80 pts each) + 3 numerical series (12 pts each) |

**Columns**: `theta_deg, ratio_analytic_T5, ratio_analytic_T10, ratio_analytic_T15` (+ numerical appended as comments)

---

## 7. Code Structure

### 7.1 Core Module
**Path**: `reproduction/plasma_params.py`

Contains:
- `PlasmaParams` class: all plasma parameters and derived quantities
- Background dispersion: `omega0_dispersion()`, `omega_magnetosonic()`, `omega_whistler()`
- Stix cold-plasma solver: `stix_parameters()`, `cold_dispersion_det()`, `solve_cold_dispersion()`, `solve_cold_dispersion_sweep()`
- Analytical growth rates: `growth_rate_eq22()`, `growth_rate_eq20()`
- Numerical growth rate: `growth_rate_eq19_numerical()`, `growth_rate_numerical()`
- Utilities: `find_max_growth_rate()`, `find_max_growth_rate_numerical()`

### 7.2 Figure Scripts
| Figure | Script Path |
|--------|-------------|
| Fig 1 | `reproduction/fig1_threshold.py` |
| Fig 2 | `reproduction/fig2_dispersion.py` |
| Fig 3 | `reproduction/fig3_growth_rate.py` |
| Fig 4 | `reproduction/fig4_angle_B.py` |
| Fig 5 | `reproduction/fig5_angle_density.py` |
| Fig 6 | `reproduction/fig6_angle_temperature.py` |

### 7.3 Output Directories
- Data files: `data/`
- Plot images: `plots/`

---

## 8. Reproduction Requirements

### 8.1 Computational Resources

| Requirement | Value | Notes |
|-------------|-------|-------|
| **Total time limit** | 1 hour | All 6 figures (Figs 4-6 are slowest due to parameter scans) |
| **Memory limit** | 2 GB | Modest memory requirements |
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
python3 fig1_threshold.py
python3 fig2_dispersion.py
python3 fig3_growth_rate.py
python3 fig4_angle_B.py
python3 fig5_angle_density.py
python3 fig6_angle_temperature.py
```
