# Reproduction Instructions for AI Agent

This document provides complete instructions for an AI agent to reproduce computational analyses from the referenced lattice QCD paper.

---

## 1. Article Information

| Field | Value |
|-------|-------|
| **Title** | Ab-initio Determination of Light Hadron Masses |
| **Authors** | S. Durr, Z. Fodor, J. Frison, C. Hoelbling, R. Hoffmann, S. D. Katz, S. Krieg, T. Kurth, L. Lellouch, T. Lippert, K.K. Szabo, G. Vulvert |
| **Journal** | Science **322**, 1224 (2008) |
| **DOI** | 10.1126/science.1163233 |
| **arXiv** | 0906.3599 |
| **Source Markdown** | `durr2008.md` |

---

## 2. Main Formulas Used in Reproduction

### 2.1 QCD Euclidean Lagrangian
$$
\mathcal{L} = -\frac{1}{2g^2} \text{Tr}\, F_{\mu\nu} F_{\mu\nu} + \bar{\psi}[\gamma_\mu(\partial_\mu + A_\mu) + m]\psi
$$
where $F_{\mu\nu} = \partial_\mu A_\nu - \partial_\nu A_\mu + [A_\mu, A_\nu]$ and $A_\mu$ is an SU(3) gauge field.

### 2.2 Chiral Extrapolation — Ratio Method (Chiral Fit)
$$
r_X = r_X^{(\text{ref})} + \alpha_X [r_\pi^2 - r_\pi^{2(\text{ref})}] + \beta_X [r_K^2 - r_K^{2(\text{ref})}] + \gamma_X r_\pi^3
$$
where $r_X = M_X / M_\Xi$ (or $M_X / M_\Omega$), $r_\pi = M_\pi / M_\Xi$, $r_K = M_K / M_\Xi$.

The $r_\pi^3$ term is the next-to-leading order chiral perturbation theory prediction (Ref. S11: Langacker & Pagels, Phys. Rev. D10, 2904, 1974).

### 2.3 Chiral Extrapolation — Taylor Fit
$$
r_X = r_X^{(\text{ref})} + \alpha_X [r_\pi^2 - r_\pi^{2(\text{ref})}] + \beta_X [r_K^2 - r_K^{2(\text{ref})}] + \delta_X [r_\pi^2 - r_\pi^{2(\text{ref})}]^2
$$

Taylor expansion around a non-singular reference point (midpoint of physical and max simulated $M_\pi$).

### 2.4 Continuum Extrapolation
$$
r_X(a) = r_X^{(\text{cont})} + c_a \cdot a^2
$$

Based on the $O(a^2)$-improved Symanzik gauge action (Ref. S1) with clover-improved Wilson fermions (Ref. S2). Leading discretization effects are $O(a^2)$.

### 2.5 Finite-Volume Correction (Type I, Luscher)
$$
M_X(L) = M_X(\infty) + c_X(M_\pi) \cdot \frac{\exp(-M_\pi L)}{(M_\pi L)^{3/2}}
$$

where $c_X(M_\pi) \propto M_\pi^2$ (Refs. S9, S10: Colangelo et al.).

### 2.6 Luscher Resonance Quantization Condition (Type II)
$$
n\pi - \delta_{11}(k) = \phi(q)
$$

where $q = kL/(2\pi)$, $\phi(q)$ is Luscher's kinematic function, and $\delta_{11}$ is the $\pi\pi$ scattering phase shift in the $I=1$, $J=1$ channel.

### 2.7 Effective Range Formula (Rho Channel)
$$
\frac{k^3}{W} \cot\delta_{11} = a + b k^2
$$
where $a = -b k_\rho^2 = 4k_\rho^5 / (M_\rho^2 \Gamma_\rho)$ and $k_\rho = \sqrt{M_\rho^2/4 - M_\pi^2}$.

### 2.8 Luscher Kinematic Function
$$
\phi(q) \approx \begin{cases} q^3 & \text{for small } q \\ \pi q^2 & \text{for } q \geq 0.1 \end{cases}
$$

---

## 3. Main Methodology

### 3.1 Hadron Spectrum (Figure 3)
1. Read lattice QCD results from Table 1 for both Xi-set and Omega-set
2. Compare to experimental (PDG) values with isospin averaging
3. Plot with combined statistical + systematic error bars (added in quadrature)
4. Resonances ($\rho$, $K^*$, $\Delta$, $\Sigma^*$, $\Xi^*$) shown with decay width bands

### 3.2 Chiral Extrapolation (Figure 2)
1. Generate mass ratio data $r_X = M_X / M_\Xi$ at each simulation point
2. For each lattice spacing ($\beta = 3.30, 3.57, 3.70$): fit $r_X$ vs $r_\pi^2$
3. Apply chiral fit ($r_\pi^3$ NLO term) and Taylor fit ($r_\pi^4$ quadratic term)
4. Extrapolate to physical $M_\pi$ and $a \to 0$ simultaneously
5. Cross mark shows continuum-extrapolated result at physical pion mass

### 3.3 Finite-Volume Analysis (Figure S4)
1. Use volume study data at $M_\pi \approx 320$ MeV, $a \approx 0.125$ fm
2. Fit pion and nucleon masses to Luscher exponential formula (Sec. 2.5)
3. Verify $M_\pi L \gtrsim 4$ rule-of-thumb for negligible finite-size effects
4. Verify $c_X(M_\pi) \propto M_\pi^2$ prediction

### 3.4 Systematic Error Analysis (Figure S5)
1. Combine $2 \times 2 \times 3 \times 2 \times 18 = 432$ fitting procedures
2. Weight each by fit quality ($\chi^2$/dof)
3. Build distribution; median gives central value, 68% CI gives systematic error
4. Bootstrap (2000 samples) for statistical error (SEM)
5. Add statistical and systematic errors in quadrature

### 3.5 Luscher Resonance Analysis (Supplementary)
1. Compute kinematic function $\phi(q)$ (approximation and exact)
2. Compute $\pi\pi$ scattering phase shift using effective range formula
3. Show free two-pion energy levels in finite volume
4. Demonstrate avoided level crossing for resonances ($\rho$, $\Delta$)

---

## 4. Input Parameters

### 4.1 Lattice Simulation Parameters
| Parameter | Symbol | Values | Description |
|-----------|--------|--------|-------------|
| Gauge coupling | $\beta = 6/g^2$ | 3.30, 3.57, 3.70 | Three lattice spacings |
| Lattice spacing | $a$ | 0.125, 0.085, 0.065 fm | Corresponding spacings |
| Pion mass range | $M_\pi$ | 190 — 650 MeV | Light quark mass variation |
| Strange quark mass | $m_s$ | $\approx$ physical | Fixed near physical value |
| Lattice sizes | $L^3 \cdot T$ | $16^3 \cdot 32$ to $48^3 \cdot 64$ | Spatial and temporal extents |
| $M_\pi L$ | | $\gtrsim 4$ | Finite-volume criterion |

### 4.2 Physical Constants
| Quantity | Symbol | Value | Description |
|----------|--------|-------|-------------|
| $\hbar c$ | | 0.19733 GeV fm | Conversion factor |
| Pion mass | $M_\pi$ | 0.135 GeV | Physical (isospin-averaged) |
| Kaon mass | $M_K$ | 0.495 GeV | Physical (isospin-averaged) |
| Xi mass | $M_\Xi$ | 1.318 GeV | Scale-setting (Xi-set) |
| Omega mass | $M_\Omega$ | 1.672 GeV | Scale-setting (Omega-set) |

### 4.3 Parameter Variations per Figure
| Figure | Varied Parameter | Values |
|--------|-----------------|--------|
| Fig 2 | Pion mass, lattice spacing | 14 simulation points at 3 lattice spacings |
| Fig 3 | (none) | Table 1 results vs experiment |
| Fig S4 | Spatial volume L | L/a = 16, 24, 32 at $M_\pi \approx 320$ MeV |
| Fig S5 | Fitting procedure | 432 combinations of analysis choices |
| Luscher | Box size L, pion momentum k | L = 1—8 fm, k = 0—500 MeV |

---

## 5. Banned Libraries

| Banned | Reason |
|--------|--------|
| Pre-built lattice QCD libraries (e.g., `qlua`, `Grid`, `chroma`) | Must implement analysis from scratch |
| GPU-accelerated libraries (cupy, pytorch) | Keep implementation in pure NumPy/SciPy |
| External QCD data analysis frameworks | Implement fitting procedures directly |

**Allowed**: `numpy`, `scipy` (optimize, constants), `matplotlib`

---

## 6. Data Files to Reproduce

### 6.1 Figure 3: Hadron Spectrum
**File**: `data/fig3_spectrum.csv`

| Property | Value |
|----------|-------|
| Description | Lattice QCD hadron masses vs experimental values |
| X-axis | Particle species (categorical) |
| Y-axis | Mass [GeV], range: 0—2.0 |
| Data series | 12 particles: $\pi$, $K$, $\rho$, $K^*$, $N$, $\Lambda$, $\Sigma$, $\Xi$, $\Delta$, $\Sigma^*$, $\Xi^*$, $\Omega$ |
| Number of points | 12 (one per particle) |

**Columns**: `particle, exp_mass_GeV, lat_mass_GeV, lat_stat_err, lat_sys_err, lat_total_err`

---

### 6.2 Figure 2: Chiral Extrapolation
**File**: `data/fig2_chiral.csv`

| Property | Value |
|----------|-------|
| Description | Mass ratios $M_X / M_\Xi$ vs $(M_\pi / M_\Xi)^2$ for N and Omega |
| X-axis | $(M_\pi / M_\Xi)^2$, range: 0—0.25 |
| Y-axis | $M_X / M_\Xi$, range: 0.6—1.4 |
| Data series | 2 particles $\times$ 3 lattice spacings = 6 series |
| Number of points | 14 per particle (one per simulation point) |

**Columns**: `particle, beta, a_fm, Mpi_sq_GeV2, r_pi_sq, r_X, r_X_err`

---

### 6.3 Figure S4: Volume Dependence
**File**: `data/figS4_volume.csv`

| Property | Value |
|----------|-------|
| Description | Pion and nucleon masses vs spatial volume at $M_\pi \approx 320$ MeV |
| X-axis | $L/a$, range: 12—40 |
| Y-axis | $aM_X$ (lattice units) |
| Data series | 2 channels (pion, nucleon) |
| Number of points | 3 per channel + fit curves |

**Columns**: `L_over_a, MpiL, aM_pi, aM_pi_err, aM_N, aM_N_err`

---

### 6.4 Figure S5: Systematic Error Distribution
**File**: `data/figS5_distribution.csv`

| Property | Value |
|----------|-------|
| Description | Weighted histogram of nucleon mass from 432 fitting procedures |
| X-axis | $M_N$ [GeV] |
| Y-axis | Weighted probability density |
| Data series | 1 histogram (50 bins) |
| Number of points | 50 |

**Columns**: `bin_center_GeV, weighted_density`

---

### 6.5 Luscher Phase Shift Analysis
**File**: `data/figS_luscher.csv` and `data/figS_luscher_phase.csv`

| Property | Value |
|----------|-------|
| Description | Luscher kinematic function $\phi(q)$ and $\pi\pi$ scattering phase shift |
| X-axis (phi) | $q = kL/(2\pi)$, range: 0—2 |
| Y-axis (phi) | $\phi(q)$ |
| X-axis (phase) | $k$ [MeV], range: 0—500 |
| Y-axis (phase) | $\delta_{11}$ [degrees] |
| Number of points | 500 (phi), 500 (phase) |

**Columns**: `q, phi_approx` (phi file); `k_MeV, delta_11_deg` (phase file)

---

## 7. Code Structure

### 7.1 Core Module
**Path**: `reproduction/lattice_qcd.py`

Contains:
- Hadron spectrum data: `EXPERIMENTAL_MASSES`, `LATTICE_XI_SET`, `LATTICE_OMEGA_SET`
- Simulation parameters: `SIMULATION_POINTS`, `LATTICE_SPACINGS`
- Chiral fits: `chiral_fit_function()`, `taylor_fit_function()`
- Continuum extrapolation: `continuum_extrapolation()`
- Finite-volume: `finite_volume_correction()`, `volume_dependence_model()`
- Luscher formulas: `luscher_phi()`, `luscher_zeta()`, `pipi_phase_shift_rho()`
- Systematic analysis: `systematic_error_analysis()`
- Data generation: `generate_chiral_data()`

### 7.2 Figure Scripts
| Figure | Script Path |
|--------|-------------|
| Fig 3 | `reproduction/fig3_spectrum.py` |
| Fig 2 | `reproduction/fig2_chiral.py` |
| Fig S4 | `reproduction/figS4_volume.py` |
| Fig S5 | `reproduction/figS5_distribution.py` |
| Luscher | `reproduction/figS_luscher.py` |

### 7.3 Output Directories
- Data files: `data/`
- Plot images: `plots/`
- Paper images: `images/`

---

## 8. Reproduction Requirements

### 8.1 Computational Resources

| Requirement | Value | Notes |
|-------------|-------|-------|
| **Total time limit** | 30 minutes | All 5 computations (mostly fast) |
| **Memory limit** | 1 GB | Modest memory requirements |
| **CPU** | Single-threaded | No parallelization required |

### 8.2 Numerical Precision
- Export data with 8 decimal places
- Masses in GeV, lattice spacing in fm

### 8.3 Environment
```bash
pip install numpy scipy matplotlib
```

### 8.4 Execution
Run scripts from the `reproduction/` directory:
```bash
cd reproduction/
python3 fig3_spectrum.py
python3 fig2_chiral.py
python3 figS4_volume.py
python3 figS5_distribution.py
python3 figS_luscher.py
```

### 8.5 Important Notes
- The original lattice QCD Monte Carlo simulations used Blue Gene supercomputers and are NOT reproducible on standard hardware.
- The scripts reproduce the **analysis procedures** applied to the simulation data: chiral extrapolation, continuum extrapolation, finite-volume corrections, and systematic error estimation.
- The hadron spectrum data in Table 1 is the primary result; the analysis scripts demonstrate how that data was obtained from raw simulation data.
- Synthetic data is generated for chiral extrapolation that is consistent with the paper's parameters and results.
