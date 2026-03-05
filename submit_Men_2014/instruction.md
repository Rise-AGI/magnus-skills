# Reproduction Instructions for AI Agent

This document provides complete instructions for an AI agent to reproduce all figures from the referenced physics paper.

---

## 1. Article Information

| Field | Value |
|-------|-------|
| **Title** | Robust topology optimization of three-dimensional photonic-crystal band-gap structures |
| **Authors** | H. Men, K. Y. K. Lee, R. M. Freund, J. Peraire, S. G. Johnson |
| **DOI** | 10.1364/OE.22.022632 |
| **Journal** | Optics Express 22(19), 22632-22657 (2014) |
| **Source Markdown** | `men2014.md` |

---

## 2. Main Formulas Used in Reproduction

### 2.1 Maxwell Eigenvalue Problem
The photonic band structure is determined by the master equation for the magnetic field H in a periodic dielectric structure:
$$
\nabla_\mathbf{k} \times \left(\frac{1}{\varepsilon(\mathbf{r})} \nabla_\mathbf{k} \times \mathbf{H}\right) = \frac{\omega^2}{c^2} \mathbf{H}
$$
where $\nabla_\mathbf{k} = \nabla + i\mathbf{k}$ for Bloch wavevector $\mathbf{k}$.

### 2.2 Plane-Wave Expansion
Expand fields in plane waves:
$$
\mathbf{H}(\mathbf{r}) = \sum_\mathbf{G} \mathbf{h}_\mathbf{G} e^{i(\mathbf{k}+\mathbf{G})\cdot\mathbf{r}}
$$
The eigenvalue equation becomes a matrix problem:
$$
\sum_{\mathbf{G}'} \eta(\mathbf{G}-\mathbf{G}') [(\mathbf{k}+\mathbf{G}) \times e_a] \cdot [(\mathbf{k}+\mathbf{G}') \times e_b] \, h_{\mathbf{G}',b} = \frac{\omega^2}{c^2} h_{\mathbf{G},a}
$$
where $\eta$ is the inverse of the epsilon Fourier matrix (Ho-Chan-Soukoulis inverse method).

### 2.3 Inverse Method (Ho-Chan-Soukoulis)
For high index contrast, use:
$$
\eta = [\varepsilon_{\mathbf{G}-\mathbf{G}'}]^{-1}
$$
rather than the direct Fourier transform of $1/\varepsilon$.

### 2.4 Fractional Band Gap
$$
\Delta\omega / \bar{\omega} = \frac{\min_\mathbf{k} \omega_{n+1}(\mathbf{k}) - \max_\mathbf{k} \omega_n(\mathbf{k})}{(\min_\mathbf{k} \omega_{n+1}(\mathbf{k}) + \max_\mathbf{k} \omega_n(\mathbf{k}))/2}
$$

### 2.5 Normalized Frequency
$$
\tilde{\omega} = \frac{\omega a}{2\pi c}
$$
where $a$ is the conventional lattice constant.

---

## 3. Main Methodology

### 3.1 Plane-Wave Expansion Method
1. Define crystal lattice vectors $\mathbf{a}_1, \mathbf{a}_2, \mathbf{a}_3$
2. Compute reciprocal vectors $\mathbf{b}_i$ satisfying $\mathbf{a}_i \cdot \mathbf{b}_j = 2\pi\delta_{ij}$
3. Generate G-vectors: $\mathbf{G} = m_1\mathbf{b}_1 + m_2\mathbf{b}_2 + m_3\mathbf{b}_3$ for $|m_i| \le N_{max}$
4. Sample $\varepsilon(\mathbf{r})$ on a real-space grid and compute FFT
5. Build $\eta$ matrix using inverse method
6. For each k-point:
   - Build Maxwell matrix M using cross-product formulation
   - Solve Hermitian eigenvalue problem (MUST use full complex matrix, NOT real part only)
   - Extract frequencies $\tilde{\omega} = \sqrt{\lambda}/(2\pi)$
7. Compute band gap from global band extrema

### 3.2 Crystal Structures
- **SC5**: Simple cubic lattice, hollow dielectric spheres + cylinders along cubic axes
- **Diamond2**: FCC primitive cell with 2-atom basis, tetrahedral cylinder bonds
- **FCC8**: FCC lattice, hollow spheres + rods to 12 nearest neighbors

### 3.3 Critical Implementation Notes
- The Maxwell matrix M is complex Hermitian; `scipy.linalg.eigh(M)` must be called on the FULL complex matrix, NOT `eigh(M.real)`. Using only the real part gives completely wrong results.
- Do NOT filter out zero eigenvalues, as this shifts band indices at Gamma and corrupts gap calculations.
- Use the inverse method for $\eta$ to achieve better convergence with high index contrast.

---

## 4. Input Parameters

### 4.1 Physical Parameters
| Parameter | Symbol | Value | Description |
|-----------|--------|-------|-------------|
| Dielectric constant | $\varepsilon_{hi}$ | 12.96 | Silicon (n=3.6) |
| Background | $\varepsilon_{lo}$ | 1.0 | Air |
| SC5 inner radius | $r_1$ | 0.14a | Air sphere radius |
| SC5 outer radius | $r_2$ | 0.36a | Dielectric sphere radius |
| SC5 cylinder radius | $r_3$ | 0.105a | Connecting cylinder radius |
| Diamond2 rod radius | $r$ | 0.1a | Cylinder radius |
| FCC8 air radius | $r_1$ | 0.12a | Air hole radius |
| FCC8 sphere radius | $r_2$ | 0.19a | Dielectric sphere radius |
| FCC8 rod radius | $r_3$ | 0.08a | Connecting rod radius |

### 4.2 Computational Parameters
| Parameter | Typical Value | Description |
|-----------|--------------|-------------|
| $N_{max}$ | 4-5 | Max Miller index (total PWs = $(2N_{max}+1)^3$) |
| Grid resolution | 32 | Real-space FFT grid per dimension |
| k-points per segment | 20 | Along high-symmetry path |

---

## 5. Banned Libraries

| Banned | Reason |
|--------|--------|
| `meep` or `mpb` | Must implement PWE solver from scratch |
| Pre-built photonic crystal libraries | Core computation must be original |
| GPU-accelerated eigensolvers | Keep in NumPy/SciPy for clarity |

Allowed: `numpy`, `scipy`, `matplotlib`

---

## 6. Data Files to Reproduce

### 6.1 Figure 5: SC5 Band Structure
**File**: `data/fig5_sc5_bands.csv`

| Property | Value |
|----------|-------|
| Description | Band structure of SC5 photonic crystal along Gamma-X-M-Gamma-R-X |
| X-axis | k-path distance |
| Y-axis | Normalized frequency $\omega a/(2\pi c)$ |
| Parameters | r1=0.14, r2=0.36, r3=0.105, eps=12.96 |
| Expected gap | ~17% between bands 5-6 |
| Number of bands | 10 |

**Columns**: `k_dist, band_1, band_2, ..., band_10`

**Gap file**: `data/fig5_sc5_gap.csv`
**Columns**: `structure, gap_bands, gap_percent, n_max, grid_res`

---

### 6.2 Figure 6: Diamond2 Band Structure
**File**: `data/fig6_diamond2_bands.csv`

| Property | Value |
|----------|-------|
| Description | Band structure of Diamond2 along FCC BZ high-symmetry path |
| Parameters | r=0.1, eps=12.96 |
| Expected gap | ~31.56% between bands 2-3 |
| Number of bands | 8 |

**Columns**: `k_dist, band_1, band_2, ..., band_8`

---

### 6.3 Figure 7: FCC8 Band Structure
**File**: `data/fig7_fcc8_bands.csv`

| Property | Value |
|----------|-------|
| Description | Band structure of FCC8 along FCC BZ high-symmetry path |
| Parameters | r1=0.12, r2=0.19, r3=0.08, eps=12.96 |
| Expected gap | ~18.3% between bands 8-9 |
| Number of bands | 14 |

**Columns**: `k_dist, band_1, band_2, ..., band_14`

---

### 6.4 Figure 8: Gap vs. Index Contrast
**File**: `data/fig8_gap_vs_contrast.csv`

| Property | Value |
|----------|-------|
| Description | Diamond2 fractional gap vs. dielectric contrast |
| X-axis | Index contrast $n_{hi}/n_{lo}$ |
| Y-axis | Fractional band gap (%) |
| Contrast values | eps_hi = 1.5, 2, 3, 4, 5, 6, 8, 10, 12.96, 16 |
| Expected threshold | gap opens around n~2 (eps~4) |

**Columns**: `eps_hi, n_contrast, gap_percent`

---

### 6.5 Figure 9: Robustness to Fabrication Error
**File**: `data/fig9_robustness.csv`

| Property | Value |
|----------|-------|
| Description | Diamond2 gap sensitivity to radius perturbation |
| X-axis | Fabrication error delta/a |
| Y-axis | Fractional band gap (%) |
| Delta range | -0.03 to 0.03 in steps of 0.005 |

**Columns**: `delta, gap_percent`

---

## 7. Code Structure

### 7.1 Core Modules
| File | Description |
|------|-------------|
| `reproduction/pwe3d.py` | 3D Plane-Wave Expansion solver |
| `reproduction/crystals.py` | Crystal structures, epsilon functions, k-paths |

### 7.2 Figure Scripts
| Figure | Script |
|--------|--------|
| Fig 5 | `reproduction/fig5_sc5_bands.py` |
| Fig 6 | `reproduction/fig6_diamond2_bands.py` |
| Fig 7 | `reproduction/fig7_fcc8_bands.py` |
| Fig 8 | `reproduction/fig8_gap_vs_contrast.py` |
| Fig 9 | `reproduction/fig9_robustness.py` |

### 7.3 Output Directories
- Data files: `data/`
- Plot images: `plots/`
- Paper images: `images/`

---

## 8. Reproduction Requirements

### 8.1 Computational Resources
| Requirement | Value | Notes |
|-------------|-------|-------|
| **Total time limit** | 3 hours | All 5 figures |
| **Memory limit** | 4 GB | n_max=5 uses ~2 GB |
| **CPU** | Single-threaded | SciPy eigh is internally threaded |

### 8.2 Numerical Precision
- Export data with 8 decimal places
- Band gaps accurate to ~1-2% relative to converged MPB values at n_max=4-5
