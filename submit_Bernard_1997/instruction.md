# Reproduction Instructions for AI Agent

This document provides complete instructions for an AI agent to reproduce figures from the referenced lattice QCD paper.

---

## 1. Article Information

| Field | Value |
|-------|-------|
| **Title** | QCD Thermodynamics with an Improved Lattice Action |
| **Authors** | C. Bernard, J. Hetrick, T. DeGrand, M. Wingate, C. DeTar, S. Gottlieb, U. Heller, K. Rummukainen, D. Toussaint, R. Sugar (MILC Collaboration) |
| **Journal** | Phys. Rev. D **56**, 5584 (1997) |
| **Preprint** | hep-lat/9612025 |
| **Source Markdown** | `bernard1997.md` |

---

## 2. Main Formulas Used in Reproduction

### 2.1 Symanzik-Improved Gauge Action
$$
S_g = \beta \sum_{\mathrm{plaq}} \frac{1}{3}\mathrm{Re\,Tr}(1 - U_{\mathrm{plaq}})
+ \beta_1 \sum_{\mathrm{rect}} \frac{1}{3}\mathrm{Re\,Tr}(1 - U_{\mathrm{rect}})
+ \beta_2 \sum_{\mathrm{twist}} \frac{1}{3}\mathrm{Re\,Tr}(1 - U_{\mathrm{twist}})
$$

### 2.2 Coupling Constants
$$
\beta = \frac{6}{g^2 u_0^4} \frac{5}{3}(1 - 0.1020 g^2 + \mathcal{O}(g^4))
$$
$$
\beta_1 = -\frac{\beta}{20 u_0^2}(1 - 0.6264 \ln(u_0)), \quad
\beta_2 = \frac{\beta}{u_0^2} 0.04335 \ln(u_0)
$$

### 2.3 Mean Link and Strong Coupling
$$
u_0 \equiv \left(\frac{1}{3}\mathrm{Re\,Tr}\langle U_{\mathrm{plaq}}\rangle\right)^{1/4}, \quad
\frac{g^2}{4\pi} \equiv -\frac{\ln(\frac{1}{3}\mathrm{Re\,Tr}\langle U_{\mathrm{plaq}}\rangle)}{3.06839}
$$

### 2.4 Clover-Improved Fermion Action (Sheikholeslami-Wohlert)
$$
S_f = S_W - \frac{\kappa}{u_0^3}\sum_x \sum_{\mu<\nu} \bar{\psi}(x)\, i\sigma_{\mu\nu} F_{\mu\nu}\, \psi(x)
$$

### 2.5 Heavy Quark Potential
$$
V(r) = V_0 + \sigma r - \frac{e}{r} - f\left(G_L(r) - \frac{1}{r}\right)
$$
where $G_L$ is the lattice Coulomb potential.

### 2.6 Sommer Parameter
$$
r_0^2 F(r_0) = 1.65, \quad F(r) = \frac{\partial V}{\partial r}
$$
corresponding to $r_0 = 0.49$ fm from potential models.

### 2.7 Critical Temperature
On an $N_t$ lattice: $T_c = 1/(N_t a)$. All scaling quantities are:
- $T_c/M_V = 1/(N_t \cdot aM_V)$
- $T_c/\sqrt{\sigma} = 1/(N_t \sqrt{a^2\sigma})$
- $r_0 T_c = (r_0/a)/N_t$

---

## 3. Main Methodology

### 3.1 Simulation Approach (Context)
The paper performed large-scale Monte Carlo simulations using the Hybrid Monte Carlo (HMC) algorithm on $8^3 \times 4$ (finite temperature) and $8^3 \times 16$ (zero temperature) lattices with SU(3) gauge fields and two flavors of dynamical fermions.

### 3.2 Reproduction Approach (This Package)
The figures reproduced here are the **scaling analysis plots** (Figs 9, 11-16) which derive physical predictions from the simulation data in the paper's Tables 1-3. These involve:
1. Extracting lattice observables from the paper's tabulated results
2. Computing derived dimensionless ratios with proper error propagation
3. Comparing improved Wilson, standard Wilson, and Kogut-Susskind formulations
4. Constructing the phase diagram from spectroscopy data

### 3.3 Key Computational Steps
- **Error propagation**: Gaussian error propagation for all derived quantities
- **Interpolation**: Linear and polynomial interpolation for phase boundary lines
- **Phase diagram**: Extrapolating $\kappa_c(\beta)$ from spectroscopy data using $(aM_{PS})^2$ vs $1/\kappa$ fits

---

## 4. Input Parameters

### 4.1 Improved Wilson Fermion Data (Table 1)
Hadron masses on $8^3 \times 16$ lattice. Points on the $N_t=4$ crossover marked with *.

| $\beta$ | $\kappa$ | $u_0$ | $aM_{PS}$ | $aM_V$ | $M_{PS}/M_V$ | On crossover |
|---------|----------|-------|-----------|---------|--------------|-------------|
| 6.40 | 0.1450 | 0.826 | 0.931(4) | 1.351(18) | 0.689(10) | No |
| *6.40 | 0.1475 | 0.828 | 0.664(8) | 1.26(6) | 0.527(26) | Yes |
| 6.60 | 0.1400 | 0.834 | 1.173(4) | 1.481(9) | 0.792(6) | No |
| *6.60 | 0.1430 | 0.841 | 0.927(4) | 1.280(8) | 0.724(6) | Yes |
| 6.60 | 0.1460 | 0.855 | 0.468(15) | 1.04(13) | 0.45(6) | No |
| 6.80 | 0.1325 | 0.842 | 1.494(3) | 1.700(7) | 0.879(4) | No |
| *6.80 | 0.1370 | 0.849 | 1.187(3) | 1.421(6) | 0.835(4) | Yes |
| 6.80 | 0.1400 | 0.857 | 0.885(8) | 1.182(16) | 0.749(12) | No |
| *7.20 | 0.1180 | 0.864 | 1.915(3) | 1.994(3) | 0.960(2) | Yes |
| *7.30 | 0.1140 | 0.8695 | 2.043(4) | 2.106(5) | 0.970(3) | Yes |

### 4.2 Heavy Quark Potential (Table 2 - Improved Wilson)

| $\beta$ | $\kappa$ | $a^2\sigma$ | $r_0/a$ |
|---------|----------|------------|---------|
| 6.40 | 0.1475 | 0.41(8) | 1.52(4) |
| 6.60 | 0.1430 | 0.42(3) | 1.77(2) |
| 6.80 | 0.1370 | 0.346(15) | 1.913(13) |
| 7.20 | 0.1180 | 0.253(6) | 2.287(13) |

### 4.3 Heavy Quark Potential (Table 3 - Kogut-Susskind)

| $\beta$ | $am_q$ | $a^2\sigma$ | $r_0/a$ | $N_t$ |
|---------|--------|------------|---------|------|
| 5.2875 | 0.025 | 0.30(3) | 1.99(4) | 4 |
| 5.3200 | 0.050 | 0.29(4) | 2.17(11) | 4 |
| 5.3750 | 0.100 | 0.288(23) | 2.20(4) | 4 |
| 5.4150 | 0.0125 | 0.130(5) | 3.14(5) | 6 |

---

## 5. Banned Libraries

| Banned | Reason |
|--------|--------|
| Pre-built lattice QCD libraries | Must compute derived quantities from scratch |
| GPU-accelerated libraries | Keep in pure NumPy/SciPy |
| External QCD analysis frameworks | Implement error propagation and fits directly |

**Allowed**: `numpy`, `scipy`, `matplotlib`

---

## 6. Data Files to Reproduce

### 6.1 Figure 9: Phase Diagram
**File**: `data/fig9_phase_diagram.csv`

| Property | Value |
|----------|-------|
| Description | Phase diagram showing $\kappa_T(\beta)$ crossover and $\kappa_c(\beta)$ critical lines |
| X-axis | $\beta$, range: 6.2-7.5 |
| Y-axis | $\kappa$, range: 0.11-0.155 |
| Data series | Crossover line (5 pts), critical line (5 pts), zero-T points (10 pts) |
| Number of points | 200 interpolated + discrete points |

**Columns**: `beta_T_interp, kappa_T_interp, beta_C_interp, kappa_C_interp`

---

### 6.2 Figure 11: Polyakov Loop vs $(aM_{PS})^2$
**File**: `data/fig11_polyakov_vs_amps2.csv`

| Property | Value |
|----------|-------|
| Description | Polyakov loop vs pseudoscalar mass squared at fixed $\beta$ |
| X-axis | $(aM_{PS})^2$ |
| Y-axis | $\langle\mathrm{Re}\,P\rangle$, range: 0-0.4 |
| Data series | Improved Wilson ($\beta=6.8$), unimproved Wilson ($\beta=4.9$) |
| Number of points | 100 per series |

**Columns**: `aMPS2_improved, Re_P_improved, aMPS2_unimproved, Re_P_unimproved`

---

### 6.3 Figure 12: $T_c/M_V$ vs $M_{PS}/M_V$
**File**: `data/fig12_tc_over_mv.csv`

| Property | Value |
|----------|-------|
| Description | Critical temperature / vector meson mass vs mass ratio |
| X-axis | $M_{PS}/M_V$, range: 0-1 |
| Y-axis | $T_c/M_V$, range: 0-0.35 |
| Data series | Improved Wilson Nt=4 (5 pts), KS Nt=4 (3 pts), KS Nt=6 (1 pt) |

**Columns**: `MPS_MV, MPS_MV_err, Tc_MV, Tc_MV_err, beta, action`

---

### 6.4 Figure 13: $T_c/\sqrt{\sigma}$ and $r_0 T_c$ vs $M_{PS}/M_V$
**File**: `data/fig13_tc_scaling.csv`

| Property | Value |
|----------|-------|
| Description | Critical temperature scaled by string tension and Sommer parameter |
| X-axis | $M_{PS}/M_V$, range: 0-1 |
| Y-axis (a) | $T_c/\sqrt{\sigma}$, range: 0.3-0.7 |
| Y-axis (b) | $r_0 T_c$, range: 0.3-0.65 |
| Data series | Improved Wilson (4 pts), KS Nt=4 (3 pts), KS Nt=6 (1 pt) |

**Columns**: `MPS_MV, MPS_MV_err, Tc_sqrtsigma, Tc_sqrtsigma_err, r0_Tc, r0_Tc_err, action`

---

### 6.5 Figure 14: $T_c/\sqrt{\sigma}$ vs $a\sqrt{\sigma}$
**File**: `data/fig14_tc_vs_a.csv`

| Property | Value |
|----------|-------|
| Description | Scaling of $T_c/\sqrt{\sigma}$ with lattice spacing |
| X-axis | $a\sqrt{\sigma}$, range: 0-0.8 |
| Y-axis | $T_c/\sqrt{\sigma}$, range: 0.3-0.6 |

**Columns**: `a_sqrtsigma, a_sqrtsigma_err, Tc_sqrtsigma, Tc_sqrtsigma_err, action`

---

### 6.6 Figure 15: $r_0\sqrt{\sigma}$ vs $a/r_0$
**File**: `data/fig15_r0_sigma.csv`

| Property | Value |
|----------|-------|
| Description | Dimensionless ratio $r_0\sqrt{\sigma}$ vs lattice spacing |
| X-axis | $a/r_0$, range: 0-0.8 |
| Y-axis | $r_0\sqrt{\sigma}$, range: 0.6-1.6 |

**Columns**: `a_over_r0, a_over_r0_err, r0_sqrtsigma, r0_sqrtsigma_err, action`

---

### 6.7 Figure 16: $M_V r_0$ vs $M_{PS}/M_V$
**File**: `data/fig16_mv_r0.csv`

| Property | Value |
|----------|-------|
| Description | Vector meson mass times Sommer parameter |
| X-axis | $M_{PS}/M_V$, range: 0-1 |
| Y-axis | $M_V r_0$, range: 0-6 |

**Columns**: `MPS_MV, MPS_MV_err, MV_r0, MV_r0_err, action`

---

## 7. Code Structure

### 7.1 Core Module
**Path**: `reproduction/paper_data.py`

Contains:
- `TABLE1`, `TABLE2`, `TABLE3`: Tabulated data from the paper
- `KS_MESON_DATA`: Kogut-Susskind meson mass ratios
- `KAPPA_T`, `KAPPA_C`: Phase boundary data
- Derived quantity functions: `compute_tc_over_mv()`, `compute_tc_over_sqrtsigma()`, `compute_r0_tc()`, `compute_a_sqrtsigma()`, `compute_r0_sqrtsigma()`, `compute_a_over_r0()`, `compute_mv_r0()`
- Phase diagram utilities: `get_crossover_data()`, `estimate_kappa_c()`, `interpolate_line()`

### 7.2 Figure Scripts
| Figure | Script Path |
|--------|-------------|
| Fig 9 | `reproduction/fig9_phase_diagram.py` |
| Fig 11 | `reproduction/fig11_polyakov_vs_amps2.py` |
| Fig 12 | `reproduction/fig12_tc_over_mv.py` |
| Fig 13 | `reproduction/fig13_tc_scaling.py` |
| Fig 14 | `reproduction/fig14_tc_vs_a.py` |
| Fig 15 | `reproduction/fig15_r0_sigma.py` |
| Fig 16 | `reproduction/fig16_mv_r0.py` |

### 7.3 Output Directories
- Data files: `data/`
- Plot images: `plots/`

---

## 8. Reproduction Requirements

### 8.1 Computational Resources

| Requirement | Value | Notes |
|-------------|-------|-------|
| **Total time limit** | 5 minutes | All 7 figures (data analysis, no MC simulation) |
| **Memory limit** | 256 MB | Very modest requirements |
| **CPU** | Single-threaded | No parallelization required |

### 8.2 Numerical Precision
- Export data with 8 decimal places
- Gaussian error propagation for all derived quantities

### 8.3 Environment
```bash
pip install numpy scipy matplotlib
```

### 8.4 Execution
Run scripts from the `reproduction/` directory:
```bash
cd reproduction/
python3 fig9_phase_diagram.py
python3 fig11_polyakov_vs_amps2.py
python3 fig12_tc_over_mv.py
python3 fig13_tc_scaling.py
python3 fig14_tc_vs_a.py
python3 fig15_r0_sigma.py
python3 fig16_mv_r0.py
```
