# Reproduction Instructions for AI Agent

This document provides complete instructions for an AI agent to reproduce all computational figures from the referenced lattice QCD paper.

---

## 1. Article Information

| Field | Value |
|-------|-------|
| **Title** | A high precision study of the QQ potential from Wilson loops in the regime of string breaking |
| **Authors** | B. Bolder, T. Struckmann, G.S. Bali, N. Eicker, T. Lippert, B. Orth, K. Schilling, P. Ueberholz (SESAM-TXL Collaboration) |
| **Year** | 2000 |
| **arXiv** | hep-lat/0005018 |
| **DOI** | 10.1103/PhysRevD.63.074504 |
| **Source Markdown** | `bolder2000.md` |

---

## 2. Main Formulas Used in Reproduction

### 2.1 R-Vector Constraint (Eq. 1)
All integer 3-vectors R = (C_min, C_mid, C_max) satisfying:
$$
R_{\min}^2 \le R^2 = C_{\min}^2 + C_{\mid}^2 + C_{\max}^2 \le R_{\max}^2
$$
with |C_min| <= |C_mid| <= |C_max| and |C_max| <= 12.

The paper uses R_min = 10 and R_max = 12*sqrt(3) ~ 20.78 in lattice units.

### 2.2 Bresenham Algorithm (Section II)
2D lattice path construction from (0,0) to (cmax, cmin):
```
cmax2 := 2 * cmax
cmin2 := 2 * cmin
chi := cmin2 - cmax
FOR i := 1 TO cmax DO
    step in max-direction
    IF chi >= 0 THEN
        chi := chi - cmax2
        step in min-direction
    ENDIF
    chi := chi + cmin2
ENDDO
```
3D generalization: combine two 2D algorithms for max-mid and max-min in one loop.

### 2.3 Cornell Potential
$$
V(R) = V_0 + \sigma a^2 \cdot R - \frac{e}{R}
$$
where:
- sigma*a^2 = K = 0.0372(8) (string tension in lattice units)
- R_0 = 5.89(3) (Sommer radius in lattice units, r_0 ~ 0.5 fm)
- e is determined by the Sommer condition: R_0^2 * dV/dR|_{R=R_0} = 1.65
  => e = 1.65 - R_0^2 * K

### 2.4 String Breaking Threshold
$$
2 m_{PS} a = 1.256(13)
$$
String breaking expected at r_c ~ 2.3 r_0 where V(r_c) = 2*m_PS*a.

### 2.5 APE Link Smearing (Eq. 2)
$$
\text{link} \to \alpha \times \text{link} + \text{staples}
$$
with alpha = 2.3, 26 iterations, followed by gauge group projection.

### 2.6 Inverse-Variance Filter (Eqs. 6-7)
$$
B_f(R_i) = N_i^{-1} \sum_{|R_j - R_i| \le R_f} \sigma_j^{-2} B(R_j)
$$
$$
N_i = \sum_{|R_j - R_i| \le R_f} \sigma_j^{-2}
$$
with filter radius R_f = 0.5.

### 2.7 Standard vs ARA Directions
Standard approach uses only lattice vectors that are multiples of:
(1,0,0), (1,1,0), (1,1,1), (2,1,0), (2,1,1), (2,2,1).
ARA uses ALL integer vectors satisfying the distance constraint.

---

## 3. Main Methodology

### 3.1 R-Vector Enumeration (Table I)
1. Generate all integer 3-vectors (C1, C2, C3) with |Ci| <= 12
2. Filter by distance constraint: 10 <= |R| <= 12*sqrt(3)
3. Group by R^2 value (vectors with the same R^2 contribute to the same R bin)
4. Count distinct R values and total vector count
5. Compare standard (6 direction families) vs ARA (all directions)

### 3.2 Bresenham Path Construction (Figure 1)
1. Implement the 2D Bresenham algorithm as described in Section II
2. Demonstrate path for C = (5, 3) as in the paper
3. Extend to 3D by combining two 2D algorithms

### 3.3 Static Potential (Figure 2)
1. Compute Cornell potential V(R) = V0 + K*R - e/R
2. Use K = 0.0372 (string tension), R0 = 5.89 (Sommer radius)
3. Determine e from Sommer condition
4. Set V0 so crossing with 2*m_PS*a occurs at r_c ~ 2.3*r0
5. Plot V(R) vs r/r0 with string breaking threshold band

### 3.4 Overlap Coefficients (Figure 3)
1. Model ground/excited state overlaps as linear functions of r/r0
2. Ground state: decreases from ~0.85 at r/r0=1.7 to ~0.31 at r/r0=3.5
3. Excited state: increases from ~0.05 to ~0.59 over same range
4. Sum approximately 0.9 (10% remainder)

### 3.5 Relative Errors on C21 (Figure 6 of paper)
1. Model exponential growth of relative errors with R
2. Compare three cases: no smearing, link smearing, source+link smearing
3. Demonstrate error reduction factors at T=1

### 3.6 Filtered Transition Correlator (Figure 7 of paper)
1. Generate synthetic ln(C21) data with R-dependent noise
2. Apply inverse-variance weighted filter (Eqs. 6-7) with R_f = 0.5
3. Demonstrate noise reduction in the string breaking region

---

## 4. Input Parameters

### 4.1 Lattice Parameters
| Parameter | Symbol | Value | Description |
|-----------|--------|-------|-------------|
| Lattice size | L_sigma^3 * L_tau | 24^3 x 40 | Spatial x temporal extent |
| Gauge coupling | beta | 5.6 | Inverse coupling |
| Sea quark hopping | kappa_sea | 0.1575 | Corresponds to m_pi/m_rho = 0.704(5) |
| Sea quark flavours | N_f | 2 | Full QCD with 2 light quarks |
| Configurations | | 184 | Separated by 1 autocorrelation length |

### 4.2 Derived Parameters
| Parameter | Symbol | Value | Description |
|-----------|--------|-------|-------------|
| Sommer radius | R_0 | 5.89(3) | In lattice units |
| Physical Sommer radius | r_0 | ~0.5 fm | |
| Lattice spacing | a | ~0.085 fm | From R_0*a ~ 0.5 fm |
| String tension | K = sigma*a^2 | 0.0372(8) | In lattice units |
| PS meson mass | 2*m_PS*a | 1.256(13) | String breaking threshold |
| Physical volume | L*a | ~2 fm | Spatial extent |

### 4.3 Smearing Parameters
| Parameter | Value | Description |
|-----------|-------|-------------|
| APE alpha | 2.3 | Link smearing weight |
| APE iterations | 26 | Number of smearing steps |
| Source alpha | 4.0 | Quark source smearing weight |
| Source iterations | 50 | Source smearing steps |
| Source cutoff R_q | 5 | |r_s - x| <= R_q |
| Filter radius R_f | 0.5 | For inverse-variance filter |

### 4.4 Wilson Loop Fit Parameters
| Parameter | Value | Description |
|-----------|-------|-------------|
| T_min (single exp) | 4 | Minimum temporal extent |
| T_max | 8 | Maximum temporal extent |
| T_min (double exp) | 1 | For two-exponential fits |
| C_max limit | 12 | Maximum component magnitude |

---

## 5. Banned Libraries

| Banned | Reason |
|--------|--------|
| Pre-built lattice QCD libraries | Must implement enumeration and algorithms from scratch |
| GPU-accelerated libraries (cupy, pytorch) | Keep in pure NumPy for clarity |
| External Monte Carlo frameworks | Not needed for these reproductions |

**Allowed**: `numpy`, `scipy`, `matplotlib`

---

## 6. Data Files to Reproduce

### 6.1 Figure 1: Bresenham Path
**File**: `data/fig1_bresenham_path.csv`

| Property | Value |
|----------|-------|
| Description | Bresenham lattice path for C = (5, 3) |
| Columns | step, x, y |
| Number of points | 6 (steps 0 to 5) |

**Columns**: `step, x, y`

---

### 6.2 Table I: R-Vector Enumeration
**File**: `data/table1_r_vectors.csv`

| Property | Value |
|----------|-------|
| Description | All distinct R^2 values and their vector counts in the ARA approach |
| Range | 100 <= R^2 <= 432 (i.e. 10 <= R <= 12*sqrt(3)) |
| Number of entries | 175 distinct R values |

**Columns**: `R_squared, R, n_vectors`

**File**: `data/table1_summary.csv`

| Property | Value |
|----------|-------|
| Description | Summary statistics: standard vs ARA |
| Expected values | Standard: 21 R-values, 302 vectors, 14.4 avg; ARA: 175 R-values, 11486 vectors, 65.6 avg |

**Columns**: `method, n_R_values, n_R_vectors, avg_vectors_per_R`

---

### 6.3 Figure 2: Static Potential
**File**: `data/fig2_potential.csv`

| Property | Value |
|----------|-------|
| Description | Cornell potential V(R) with string breaking threshold |
| X-axis | r/r_0 (range: 1.7 to 3.5) |
| Y-axis | V(R)*a (lattice units) |
| Number of points | 200 |
| Key feature | Crossing at r_c ~ 2.3 r_0 |

**Columns**: `r_over_r0, R_lattice, V_lattice`

---

### 6.4 Figure 3: Overlap Coefficients
**File**: `data/fig3_overlaps.csv`

| Property | Value |
|----------|-------|
| Description | Ground and excited state overlap coefficients vs r/r0 |
| X-axis | r/r_0 (range: 1.7 to 3.5) |
| Y-axis | Overlap coefficient (range: 0 to 1) |
| Data series | c0 (ground), c1 (excited), sum |
| Number of points | 200 per series |

**Columns**: `r_over_r0, c0_ground, c1_excited, c_sum`

---

### 6.5 Figure 4: Relative Errors
**File**: `data/fig4_relative_errors.csv`

| Property | Value |
|----------|-------|
| Description | Relative errors on C21(R,T=1) for different smearing methods |
| X-axis | R (lattice units, range: 8 to 20) |
| Y-axis | Relative error (log scale) |
| Data series | No smearing, link smearing, source+link smearing |
| Number of points | 24 per series |

**Columns**: `R, err_no_smearing, err_link_smearing, err_source_link_smearing`

---

### 6.6 Figure 5: Filtered Correlator
**File**: `data/fig5_filtering.csv`

| Property | Value |
|----------|-------|
| Description | ln(C21) before and after inverse-variance filtering at T=5 |
| X-axis | R (lattice units, range: 8 to 20) |
| Y-axis | ln(C21) |
| Filter radius | R_f = 0.5 |
| Number of points | 175 |

**Columns**: `R, B_true, B_noisy, sigma, B_filtered`

---

## 7. Code Structure

### 7.1 Core Module
**Path**: `reproduction/lattice_vectors.py`

Contains:
- `enumerate_r_vectors()`: enumerate all integer 3-vectors in distance range
- `classify_vectors()`: separate standard vs ARA vector sets
- `bresenham_2d()`, `bresenham_3d()`, `bresenham_general()`: lattice path construction
- `cornell_potential()`: V(R) = V0 + K*R - e/R
- `string_breaking_threshold()`: 2*m_PS*a computation

### 7.2 Figure Scripts
| Figure | Script Path |
|--------|-------------|
| Fig 1 (Bresenham path) | `reproduction/fig1_bresenham.py` |
| Table I (R-vector enumeration) | `reproduction/table1_enumeration.py` |
| Fig 2 (Static potential) | `reproduction/fig2_potential.py` |
| Fig 3 (Overlap coefficients) | `reproduction/fig3_overlaps.py` |
| Fig 4 (Relative errors) | `reproduction/fig4_relative_errors.py` |
| Fig 5 (Filtered correlator) | `reproduction/fig5_filtering.py` |

### 7.3 Output Directories
- Data files: `data/`
- Plot images: `plots/`

---

## 8. Reproduction Requirements

### 8.1 Computational Resources

| Requirement | Value | Notes |
|-------------|-------|-------|
| **Total time limit** | 5 minutes | All figures are fast (no Monte Carlo needed) |
| **Memory limit** | 1 GB | R-vector enumeration is the largest computation |
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
python3 fig1_bresenham.py
python3 table1_enumeration.py
python3 fig2_potential.py
python3 fig3_overlaps.py
python3 fig4_relative_errors.py
python3 fig5_filtering.py
```
