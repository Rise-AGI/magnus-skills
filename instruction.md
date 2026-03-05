# Reproduction Instructions for AI Agent

This document provides complete instructions for an AI agent to reproduce all computational results from the referenced physics paper.

---

## 1. Article Information

| Field | Value |
|-------|-------|
| **Title** | Surface plasmon polariton propagation around bends at a metal-dielectric interface |
| **Authors** | K. Hasegawa, J. U. Nockel, M. Deutsch |
| **Journal** | Applied Physics Letters **84**, 1835 (2004) |
| **DOI** | 10.1063/1.1675942 |
| **Source Markdown** | `hasegawa2004.md` |

---

## 2. Main Formulas Used in Reproduction

### 2.1 Silver Dielectric Function (Drude Model)
$$
\varepsilon_i(\omega) = \varepsilon_\infty - \frac{\omega_p^2}{\omega(\omega + i\gamma)}
$$
Parameters: $\varepsilon_\infty = 3.7$, $\omega_p = 9.01$ eV, $\gamma = 0.048$ eV.

In terms of wavelength: $E = 1240/\lambda_{\text{nm}}$ eV, then $\varepsilon_i = 3.7 - 9.01^2 / (E(E + 0.048i))$.

### 2.2 SPP Wavevector on Flat Interface
$$
k = \frac{\omega}{c} \sqrt{\frac{\varepsilon_i \varepsilon_o}{\varepsilon_i + \varepsilon_o}}
$$
where $\varepsilon_o = 1$ for air.

### 2.3 Decay Constants
$$
\gamma_i = -\frac{\omega \varepsilon_i}{c} \sqrt{\frac{-1}{\varepsilon_i + \varepsilon_o}}, \quad
\gamma_o = \frac{\omega \varepsilon_o}{c} \sqrt{\frac{-1}{\varepsilon_i + \varepsilon_o}}
$$

### 2.4 Boundary Matching Equation (Eq. 2)
$$
\frac{1}{k_i} \frac{J'_n(k_i R)}{J_n(k_i R)} = \frac{1}{k_o} \frac{H_n^{(1)'}(k_o R)}{H_n^{(1)}(k_o R)}
$$
where $k_i = \omega\sqrt{\varepsilon_i}/c$, $k_o = \omega\sqrt{\varepsilon_o}/c$.

**CRITICAL**: Both $n$ (the mode index) and the arguments $k_i R$, $k_o R$ are complex. Standard `scipy.special.jv` does NOT support complex order with complex argument. You MUST use `mpmath.besselj` and `mpmath.bessely` instead.

### 2.5 Solving for the Fundamental Mode $m$
The fundamental mode $m$ minimizes the field profile mismatch at the boundary. In the short-wavelength limit ($\omega R/c \gg 1$), $m \approx kR$ serves as an excellent initial guess.

Use Newton's method on the boundary equation:
1. Start with $m_0 = kR$
2. Compute residual $f(m) = \text{LHS} - \text{RHS}$ of Eq. 2
3. Compute derivative numerically: $f'(m) \approx (f(m+h) - f(m))/h$ with $h = m \times 10^{-8}$
4. Update: $m_{n+1} = m_n - f(m_n)/f'(m_n)$
5. Converge when $|f(m)| < 10^{-12}$

For the recurrence relations:
- $J'_n(z) = J_{n-1}(z) - (n/z) J_n(z)$
- $H^{(1)}_n(z) = J_n(z) + i Y_n(z)$

### 2.6 Transmittance (Eq. 3)
$$
T = \left| \frac{4mkR}{-e^{im\theta}(m - kR)^2 + e^{-im\theta}(m + kR)^2} \right|^2
$$

### 2.7 Reflectance
$$
R = \left| \frac{-(m^2 - (kR)^2) \cdot 2i \sin(m\theta)}{-e^{im\theta}(m - kR)^2 + e^{-im\theta}(m + kR)^2} \right|^2
$$

### 2.8 Upper Bound Transmittance (Eq. 6)
$$
T_u = \exp(-2 \operatorname{Im}[m] \theta)
$$

### 2.9 Coupling Efficiency (Eq. 5)
$$
\Delta^2 = \frac{\int_R^{R + \eta/\gamma_o} \left| \frac{H_m^{(1)}(k_o r)}{H_m^{(1)}(k_o R)} - e^{-\gamma_o(r-R)} \right|^2 dr}{\int_R^{R + \eta/\gamma_o} |e^{-\gamma_o(r-R)}|^2 dr}
$$
with $\eta = 3$.

---

## 3. Main Methodology

### 3.1 Key Physical Insight
The problem is analogous to 1D scattering from a finite potential well. The SPP field in the bend region (Region II) is described by cylindrical Bessel functions. Matching at the two boundary lines gives T, R, and radiation loss.

### 3.2 Numerical Approach
1. Compute $\varepsilon_i$ from the Drude model at the given wavelength
2. Compute $k$, $k_i$, $k_o$ from $\omega$ and dielectric constants
3. Solve Eq. 2 for complex mode index $m$ using Newton iteration with `mpmath`
4. Evaluate T, R, Tu from the formulas above

### 3.3 Key Numerical Details
- **mpmath precision**: Use `mpmath.mp.dps = 15` (15 decimal places). Higher precision for small $\omega R/c$.
- **Continuation**: When scanning over R, use the previous solution as the initial guess for the next point. Scale: $m_{\text{guess}} = m_{\text{prev}} \times R / R_{\text{prev}}$.
- **Newton step**: Use $h = m \times 10^{-8}$ for the numerical derivative. If $|h| < 10^{-12}$, use $h = 10^{-8} + 10^{-8}i$.
- **Convergence**: Typically 3-8 iterations per point.

### 3.4 Important Physics Note
$\operatorname{Im}[m]$ has two contributions:
1. **Absorption loss** from $\operatorname{Im}[kR]$ (dominant at large R)
2. **Radiation loss** from tunneling through the angular momentum barrier (dominant at small R)

Crucially, $\operatorname{Im}[m] < \operatorname{Im}[kR]$ because the curved geometry reduces field overlap with the metal. This means propagation on curved surfaces can have LOWER losses than flat surfaces — a key result of the paper.

---

## 4. Input Parameters

### 4.1 Physical Constants
| Parameter | Value | Description |
|-----------|-------|-------------|
| $c$ | $2.998 \times 10^8$ m/s | Speed of light |
| $\varepsilon_\infty$ | 3.7 | Silver high-frequency dielectric |
| $\omega_p$ | 9.01 eV | Silver plasma frequency |
| $\gamma$ | 0.048 eV | Silver Drude damping |
| $\varepsilon_o$ | 1.0 | Air dielectric constant |

### 4.2 Paper's Specific Test Case
| Parameter | Value |
|-----------|-------|
| $\varepsilon_i$ | $-15 + 0.5i$ |
| $\omega R / c$ | 800 |
| $\theta$ | $90°$ ($\pi/2$ radians) |

Expected results for this case:
- **Lossless** ($\varepsilon_i = -15$): $T = 0.997$, $R = 1.19 \times 10^{-8}$, $P \approx 0.003$
- **Lossy** ($\varepsilon_i = -15 + 0.5i$): $T = 0.0516$, $R = 1.18 \times 10^{-6}$, $P \approx 0.00282$

### 4.3 Parameter Variations per Figure
| Figure | Computation | Parameters |
|--------|-------------|------------|
| Fig 1 | Field profiles | $\omega R/c = 800$, $\varepsilon_i = -15 + 0.5i$, $\theta = 90°$ |
| Fig 2 (main) | $T_u$ vs $R$ | $\lambda = 500, 600, 700$ nm; $R = 1$--$200$ $\mu$m; $\theta = 90°$ |
| Fig 2 (inset) | $T_u$ map | $\lambda = 400$--$800$ nm; $R = 1$--$200$ $\mu$m |

---

## 5. Banned Libraries

| Banned | Reason |
|--------|--------|
| Pre-built plasmonics libraries | Must implement boundary equation from scratch |
| GPU-accelerated libraries | Pure NumPy/SciPy/mpmath for clarity |

**Allowed**: `numpy`, `scipy` (constants, optimization), `matplotlib`, `mpmath` (for complex-order Bessel functions)

---

## 6. Data Files to Reproduce

### 6.1 Figure 2 (Main): Tu vs Bend Radius
**File**: `data/fig2_transmittance_bound.csv`

| Property | Value |
|----------|-------|
| Description | Upper bound transmittance $T_u$ vs bend radius for three wavelengths |
| X-axis | $R$ [$\mu$m], range: 1--200 (log scale) |
| Y-axis | $T_u$, range: 0--1 |
| Data series | 3 curves: $\lambda = 500, 600, 700$ nm |
| Number of points | 120 per series |

**Columns**: `R_um, Tu_500nm, Tu_600nm, Tu_700nm`

---

### 6.2 Figure 2 (Inset): Tu Grayscale Map
**File**: `data/fig2_inset_Tu_map.csv`

| Property | Value |
|----------|-------|
| Description | $T_u$ as function of $R$ and $\lambda$ |
| X-axis | $R$ [$\mu$m], range: 1--200 (log scale) |
| Y-axis | $\lambda$ [nm], range: 400--800 |
| Grid | 25 wavelengths $\times$ 40 radii |

---

### 6.3 Scattering Coefficients
**File**: `data/scattering_coefficients.csv`

| Property | Value |
|----------|-------|
| Description | T, R, A, Tu for lossless and lossy cases at $\omega R/c = 800$ |
| Cases | 2: lossless ($\varepsilon_i = -15$) and lossy ($\varepsilon_i = -15 + 0.5i$) |

**Columns**: `case, eps_i_real, eps_i_imag, T, R_coeff, A, Tu, flat_surface, m_real, m_imag`

---

### 6.4 Field Profiles
**File**: `data/fig1_field_profiles.csv`

| Property | Value |
|----------|-------|
| Description | SPP field profiles for flat and cylindrical geometries |
| X-axis | $(r - R)/R$, range: $-0.08$ to $0.08$ |
| Y-axis | Normalized field amplitude |

**Columns**: `r_minus_R_over_R, flat_spp_field, cylindrical_mode_field`

---

### 6.5 Coupling Efficiency
**File**: `data/coupling_efficiency.csv`

| Property | Value |
|----------|-------|
| Description | $\Delta^2$ vs $\omega R/c$ for lossless metal |
| Points | 9 values: $\omega R/c = 50, 100, 150, 200, 300, 400, 600, 800, 1000$ |

**Columns**: `omega_R_over_c, Delta2`

---

## 7. Code Structure

### 7.1 Core Module
**Path**: `reproduction/spp_core.py`

Contains:
- `silver_epsilon()`: Drude model for silver
- `spp_wavevector()`: SPP k on flat surface
- `solve_mode_index()`: Newton solver for boundary equation using mpmath
- `transmittance()`, `reflectance()`: Eq. 3 and reflectance formula
- `upper_bound_transmittance()`: Eq. 6
- `compute_Tu_vs_R()`: Scan over R with continuation

### 7.2 Figure Scripts
| Figure | Script Path |
|--------|-------------|
| Fig 1 + coefficients | `reproduction/fig1_coefficients.py` |
| Fig 2 | `reproduction/fig2_transmittance_bound.py` |

### 7.3 Output Directories
- Data files: `data/`
- Plot images: `plots/`

---

## 8. Reproduction Requirements

### 8.1 Computational Resources
| Requirement | Value | Notes |
|-------------|-------|-------|
| **Total time limit** | 3 hours | Fig 2 with inset is the slowest (~30 min) |
| **Memory limit** | 2 GB | Modest requirements |
| **CPU** | Single-threaded | mpmath is inherently serial |

### 8.2 Numerical Precision
- Export data with 8 decimal places
- mpmath precision: 15 decimal places

### 8.3 Environment
```bash
pip install numpy scipy matplotlib mpmath
```

### 8.4 Execution
Run scripts from the `reproduction/` directory:
```bash
cd reproduction/
python3 fig1_coefficients.py
python3 fig2_transmittance_bound.py
```
