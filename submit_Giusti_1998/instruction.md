# Reproduction Instructions for AI Agent

This document provides complete instructions for an AI agent to reproduce all computational results from the referenced lattice QCD paper.

---

## 1. Article Information

| Field | Value |
|-------|-------|
| **Title** | The QCD Chiral Condensate from the Lattice |
| **Authors** | L. Giusti, F. Rapuano, M. Talevi, A. Vladikas |
| **Journal** | Nucl. Phys. B538 (1999) 249-277 |
| **arXiv** | hep-lat/9807014 |
| **DOI** | 10.1016/S0550-3213(98)00659-2 |
| **Source Markdown** | `giusti1998.md` |

---

## 2. Main Formulas Used in Reproduction

### 2.1 Beta-Function Coefficients (Eq. B.4)
The QCD beta-function at NNLO in the MSbar scheme:
$$
\frac{\beta(\alpha_s)}{4\pi} = -\sum_{i=0}^{\infty} \beta_i \left(\frac{\alpha_s}{4\pi}\right)^{i+2}
$$

For $N_c = 3$, $N_f = 0$ (quenched):
$$
\beta_0 = \frac{11}{3} N_c = 11, \quad
\beta_1 = \frac{34}{3} N_c^2 = 102, \quad
\beta_2^{\overline{MS}} = \frac{2857}{54} N_c^3 = 1428.5
$$

### 2.2 Quark Mass Anomalous Dimension (Eq. B.5)
$$
\gamma_m(\alpha_s) = \sum_{i=0}^{\infty} \gamma_m^{(i)} \left(\frac{\alpha_s}{4\pi}\right)^{i+1}
$$

Coefficients for $N_c = 3$, $N_f = 0$:
$$
\gamma_m^{(0)} = 3 \frac{N_c^2 - 1}{N_c} = 8
$$
$$
\gamma_m^{(1)} = \frac{N_c^2-1}{N_c^2}\left(-\frac{3}{4} + \frac{203}{12}N_c^2\right) = 134.667
$$
$$
\gamma_m^{(2)} = \frac{N_c^2-1}{N_c^3}\left[\frac{129}{8} - \frac{129}{8}N_c^2 + \frac{11413}{108}N_c^4\right] = 2498.0
$$

### 2.3 Running Coupling at NNLO (Eq. B.7)
$$
\frac{\alpha_s^{\overline{MS}}}{4\pi}(q^2) = \frac{1}{\beta_0 L} - \frac{\beta_1}{\beta_0^3} \frac{\ln L}{L^2} + \frac{1}{\beta_0^5 L^3}\left(\beta_1^2 \ln^2 L - \beta_1^2 \ln L + \beta_2 \beta_0 - \beta_1^2\right)
$$
where $L = \ln(q^2)$, $q^2 = (\mu/\Lambda_{QCD}^{\overline{MS}})^2$.

**Key parameter**: $\Lambda_{QCD}^{\overline{MS}} = 0.251 \pm 0.021$ GeV (quenched, from ref. [37]).

### 2.4 RI/MSbar Matching Coefficient (Eq. B.1-B.3)
$$
\Delta Z^{RI/\overline{MS}} = 1 + \frac{\alpha_s}{4\pi} C^{(1)} + \left(\frac{\alpha_s}{4\pi}\right)^2 C^{(2)}
$$

For $N_c = 3$, $N_f = 0$:
$$
C^{(1)} = \frac{8(N_c^2-1)}{4 N_c} = \frac{16}{3} \approx 5.333
$$
$$
C^{(2)} = \frac{N_c^2-1}{96 N_c^2}\left(-309 + 3029 N_c^2 - 288\zeta_3 - 576 N_c^2 \zeta_3\right) \approx 188.651
$$

where $\zeta_3 = 1.202056903...$

### 2.5 RG Evolution Coefficient (Eq. B.6)
$$
c_S^{\overline{MS}}(\mu) = \alpha_s(\mu)^{\bar{\gamma}_S^{(0)}} \left\{ 1 + \frac{\alpha_s}{4\pi}(\bar{\gamma}_S^{(1)} - \bar{\beta}_1 \bar{\gamma}_S^{(0)}) + \frac{1}{2}\left(\frac{\alpha_s}{4\pi}\right)^2 [\ldots] \right\}
$$

where $\bar{\beta}_i = \beta_i/\beta_0$ and $\bar{\gamma}_S^{(i)} = -\gamma_m^{(i)}/(2\beta_0)$.

The chiral condensate at different scales is related by:
$$
\langle\bar{\psi}\psi\rangle^{\overline{MS}}(\mu') = \frac{c_S^{\overline{MS}}(\mu')}{c_S^{\overline{MS}}(\mu)} \langle\bar{\psi}\psi\rangle^{\overline{MS}}(\mu)
$$

### 2.6 Chiral Condensate from GMOR Relation (Eqs. 41-42)

**Method 1** (Wilson action only, uses $Z_S$):
$$
\frac{1}{N_f}\langle\bar{\psi}\psi\rangle_1 = -\frac{1}{2} a^{-1} f_\chi^2 Z_S C^{HS}
$$

**Method 2** (both actions, uses $Z_P/Z_A$):
$$
\frac{1}{N_f}\langle\bar{\psi}\psi\rangle_2 = -\frac{1}{2} a^{-1} f_\chi^2 \frac{Z_P}{Z_A} C^{AWI}
$$

where $f_\chi = 0.1282$ GeV is the pion decay constant in the chiral limit.

---

## 3. Main Methodology

### 3.1 Perturbative QCD Calculations
1. Compute beta-function coefficients $\beta_0, \beta_1, \beta_2$ for quenched QCD ($N_c=3$, $N_f=0$)
2. Compute anomalous dimension coefficients $\gamma_m^{(0)}, \gamma_m^{(1)}, \gamma_m^{(2)}$
3. Compute the running coupling $\alpha_s^{\overline{MS}}(\mu)$ at NNLO using Eq. B.7
4. Compute the RI/MSbar matching coefficient $\Delta Z^{RI/\overline{MS}}$ at NNLO using Eq. B.1
5. Compute the evolution coefficient $c_S^{\overline{MS}}(\mu)$ at NNLO using Eq. B.6

### 3.2 Chiral Condensate Determination
1. Use lattice input data from Tables 1-4 (inverse lattice spacing $a^{-1}$, renormalization constants $Z_S$, $Z_P$, $Z_A$, slopes $C^{HS}$, $C^{AWI}$)
2. Compute the condensate in the RI scheme at the lattice scale $\mu \simeq a^{-1}$
3. Convert to MSbar scheme using RI/MSbar matching
4. Run to $\mu = 2$ GeV and $\mu = 1$ GeV using NNLO RG evolution
5. Compare across different $\beta$ values and actions (Wilson vs Clover)

### 3.3 Key Numerical Details
- **Quenched approximation**: $N_f = 0$ throughout all perturbative formulas
- **$\Lambda_{QCD}$**: Use $\Lambda_{QCD}^{\overline{MS}} = 0.251 \pm 0.021$ GeV
- **Lattice scales**: $a^{-1}$ from $K^*$ meson mass; values in Table 1
- **Non-perturbative renormalization**: $Z_S$, $Z_P$, $Z_A$ values from Table 2
- **Overflow/underflow protection**: None needed (all calculations are perturbative)

---

## 4. Input Parameters

### 4.1 Physical Constants
| Parameter | Symbol | Value | Description |
|-----------|--------|-------|-------------|
| Colors | $N_c$ | 3 | SU(3) gauge group |
| Flavors | $N_f$ | 0 | Quenched approximation |
| $\Lambda_{QCD}$ | $\Lambda_{QCD}^{\overline{MS}}$ | 0.251 GeV | Quenched QCD scale |
| $\Lambda_{QCD}$ error | | 0.021 GeV | Uncertainty |
| Pion decay constant | $f_\chi$ | 0.1282 GeV | Chiral limit value |
| Riemann zeta | $\zeta_3$ | 1.2020569... | Used in C^(2) |

### 4.2 Lattice Input (Tables 1-2)

| Run | Action | $\beta$ | $a^{-1}$ [GeV] | $Z_S$ | $Z_P$ | $Z_A$ |
|-----|--------|---------|-----------------|--------|--------|--------|
| C60 | Clover | 6.0 | 2.12(6) | 0.83(2) | 0.41(6) | 1.05(3) |
| W60 | Wilson | 6.0 | 2.26(5) | 0.68(1) | 0.45(6) | 0.81(1) |
| C62 | Clover | 6.2 | 2.7(1) | 0.85(2) | 0.47(5) | 1.02(2) |
| W62 | Wilson | 6.2 | 3.00(9) | 0.72(1) | 0.50(5) | 0.81(1) |
| C64 | Clover | 6.4 | 4.0(2) | 0.85(2) | 0.55(3) | 1.01(1) |
| W64 | Wilson | 6.4 | 4.1(2) | 0.74(1) | 0.57(4) | 0.82(1) |

### 4.3 Lattice Results (Table 4)

| Run | $C^{HS}$ | $C^{AWI}$ |
|-----|----------|-----------|
| C60a | 2.98(8) | 3.9(1) |
| C60b | 3.04(7) | 4.1(1) |
| W60 | 2.40(5) | 3.01(7) |
| C62 | 2.9(1) | 3.7(2) |
| W62 | 2.52(8) | 2.98(9) |
| C64 | 3.5(2) | 4.2(2) |
| W64 | 2.9(1) | 3.2(1) |

### 4.4 Parameter Variations per Computation

| Computation | Varied Parameter | Values |
|-------------|-----------------|--------|
| Fig 1 | $\mu$ | 0.8 - 10.0 GeV (running coupling) |
| Fig 2 | $\mu$ | 1.0 - 10.0 GeV (matching & evolution) |
| Fig 3 | Run/Action/$\beta$ | All runs from Table 4 (condensate comparison) |
| Fig 4 | $\mu$ | 1.0 - 10.0 GeV (scale dependence, LO/NLO/NNLO) |

---

## 5. Banned Libraries

| Banned | Reason |
|--------|--------|
| Pre-built QCD libraries (e.g., `pyqcd`, `qcdlib`) | Must implement perturbative formulas from scratch |
| GPU-accelerated libraries (cupy, pytorch) | Keep implementation in pure NumPy for clarity |
| External lattice QCD analysis packages | Implement renormalization and running directly |

**Allowed**: `numpy`, `scipy` (constants), `matplotlib`

---

## 6. Data Files to Reproduce

### 6.1 Computation 1: Running Coupling
**File**: `data/fig1_running_coupling.csv`

| Property | Value |
|----------|-------|
| Description | Running coupling $\alpha_s^{\overline{MS}}(\mu)$ at NNLO |
| X-axis | $\mu$ [GeV], range: 0.8 - 10.0 |
| Y-axis | $\alpha_s$, range: 0 - 0.8 |
| Data series | Central value + $\Lambda_{QCD}$ error band |
| Number of points | 500 |

**Columns**: `mu_GeV, alpha_s_central, alpha_s_low, alpha_s_high`

**Verification**: $\alpha_s(2\text{ GeV}) \approx 0.208$, $\alpha_s(3\text{ GeV}) \approx 0.177$

---

### 6.2 Computation 2: RI/MSbar Matching and RG Evolution
**File**: `data/fig2_matching_evolution.csv`

| Property | Value |
|----------|-------|
| Description | RI/MSbar matching coefficient and RG evolution |
| X-axis | $\mu$ [GeV], range: 1.0 - 10.0 |
| Y-axis | $\Delta Z$, $c_S$, and RG ratio |
| Data series | 3 columns |
| Number of points | 400 |

**Columns**: `mu_GeV, delta_z_ri_msbar, c_s_msbar, rg_ratio_to_2GeV`

**Verification**: $\Delta Z(2\text{ GeV}) \approx 1.140$

---

### 6.3 Computation 3: Chiral Condensate Comparison
**File**: `data/fig3_condensate_comparison.csv`

| Property | Value |
|----------|-------|
| Description | Chiral condensate from lattice data, all runs |
| Data series | Both methods, all runs |

**Columns**: `run, method, action, beta, computed, paper_value, paper_stat_err, paper_sys_err`

**Verification**: Best estimate $-\langle\bar{\psi}\psi\rangle^{\overline{MS}}(2\text{ GeV}) = 0.0147(8)(16)(12)$ GeV$^3$

---

### 6.4 Computation 4: Scale Dependence
**File**: `data/fig4_scale_dependence.csv`

| Property | Value |
|----------|-------|
| Description | RG running of chiral condensate at LO, NLO, NNLO |
| X-axis | $\mu$ [GeV], range: 1.0 - 10.0 |
| Y-axis | $-\langle\bar{\psi}\psi\rangle^{\overline{MS}}(\mu)$ [GeV$^3$] |
| Data series | 3 curves (LO, NLO, NNLO) |
| Number of points | 400 |

**Columns**: `mu_GeV, condensate_lo_GeV3, condensate_nlo_GeV3, condensate_nnlo_GeV3`

**Verification**:
- At $\mu = 1$ GeV: $0.0124$ GeV$^3$ (paper Eq. 46)
- RGI: $0.0088$ GeV$^3$ (paper Eq. 45)

---

## 7. Code Structure

### 7.1 Core Module
**Path**: `reproduction/qcd_perturbative.py`

Contains:
- `beta_coefficients()`: $\beta_0, \beta_1, \beta_2$ for arbitrary $N_c, N_f$
- `gamma_m_coefficients()`: $\gamma_m^{(0)}, \gamma_m^{(1)}, \gamma_m^{(2)}$
- `alpha_s_nnlo()`, `alpha_s_value()`: Running coupling at NNLO
- `ri_msbar_matching_coefficients()`: $C^{(1)}, C^{(2)}$
- `delta_z_ri_msbar()`: RI/MSbar matching coefficient
- `evolution_coefficient()`: $c_S^{\overline{MS}}(\mu)$ at NNLO
- `run_condensate()`: RG running of the chiral condensate
- `compute_condensate_method1()`, `compute_condensate_method2()`: Lattice to physical
- All lattice data from Tables 1-4 as dictionaries

### 7.2 Computation Scripts
| Computation | Script Path |
|-------------|-------------|
| Fig 1 | `reproduction/fig1_running_coupling.py` |
| Fig 2 | `reproduction/fig2_matching_evolution.py` |
| Fig 3 | `reproduction/fig3_condensate_comparison.py` |
| Fig 4 | `reproduction/fig4_scale_dependence.py` |

### 7.3 Output Directories
- Data files: `data/`
- Plot images: `plots/`

---

## 8. Reproduction Requirements

### 8.1 Computational Resources

| Requirement | Value | Notes |
|-------------|-------|-------|
| **Total time limit** | 1 minute | All computations are analytical/perturbative |
| **Memory limit** | 512 MB | Minimal memory requirements |
| **CPU** | Single-threaded | No parallelization needed |

### 8.2 Numerical Precision
- Export data with 8 decimal places
- All perturbative series truncated at NNLO ($\mathcal{O}(\alpha_s^3)$ corrections neglected)

### 8.3 Environment
```bash
pip install numpy scipy matplotlib
```

### 8.4 Execution
Run scripts from the `reproduction/` directory:
```bash
cd reproduction/
python3 fig1_running_coupling.py
python3 fig2_matching_evolution.py
python3 fig3_condensate_comparison.py
python3 fig4_scale_dependence.py
```
