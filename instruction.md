# Reproduction Instructions for AI Agent

This document provides complete instructions for an AI agent to reproduce all computationally reproducible figures from the referenced physics paper.

---

## 1. Article Information

| Field | Value |
|-------|-------|
| **Title** | First principles phonon calculations in materials science |
| **Authors** | Atsushi Togo, Isao Tanaka |
| **Journal** | Scripta Materialia **108**, 1-5 (2015) |
| **DOI** | 10.1016/j.scriptamat.2015.07.021 |
| **Source Markdown** | `togo2015.md` |

---

## 2. Main Formulas Used in Reproduction

### 2.1 Dynamical Matrix (Eq. 3-4)
$$
D_{\kappa\kappa'}^{\alpha\beta}(\mathbf{q}) = \sum_{l'} \frac{\Phi_{\alpha\beta}(0\kappa, l'\kappa')}{\sqrt{m_\kappa m_{\kappa'}}} e^{i\mathbf{q} \cdot [\mathbf{r}(l'\kappa') - \mathbf{r}(0\kappa)]}
$$
Eigenvalue problem: $D(\mathbf{q}) \mathbf{e}_{\mathbf{q}j} = \omega_{\mathbf{q}j}^2 \mathbf{e}_{\mathbf{q}j}$

### 2.2 Born-von Karman Force Constants
For a neighbor at position $\mathbf{R}$ with unit direction $\hat{d} = \mathbf{R}/|\mathbf{R}|$:
$$
\Phi_{\alpha\beta}(\mathbf{R}) = \alpha_n \hat{d}_\alpha \hat{d}_\beta + \beta_n (\delta_{\alpha\beta} - \hat{d}_\alpha \hat{d}_\beta)
$$
where $\alpha_n$ (longitudinal) and $\beta_n$ (transverse) are the force constants for the $n$-th neighbor shell.

The dynamical matrix for a monatomic lattice:
$$
D_{\alpha\beta}(\mathbf{q}) = \frac{1}{M} \sum_{\mathbf{R}} \Phi_{\alpha\beta}(\mathbf{R}) [1 - e^{i\mathbf{q}\cdot\mathbf{R}}]
$$

### 2.3 Phonon Density of States (Eq. 6)
$$
g(\omega) = \frac{1}{N} \sum_{\mathbf{q}j} \delta(\omega - \omega_{\mathbf{q}j})
$$
Normalized so that $\int g(\omega) d\omega = 3n_a$ (3 modes per atom).

### 2.4 Thermal Properties

**Energy** (Eq. 8):
$$
E = \sum_{\mathbf{q}j} \hbar\omega_{\mathbf{q}j} \left[\frac{1}{2} + \frac{1}{e^{\hbar\omega_{\mathbf{q}j}/k_BT} - 1}\right]
$$

**Heat capacity** $C_V$ (Eq. 9):
$$
C_V = \sum_{\mathbf{q}j} k_B \left(\frac{\hbar\omega_{\mathbf{q}j}}{k_BT}\right)^2 \frac{e^{\hbar\omega/k_BT}}{(e^{\hbar\omega/k_BT}-1)^2}
$$

**Helmholtz free energy** (Eq. 10):
$$
F = \frac{1}{2}\sum_{\mathbf{q}j} \hbar\omega_{\mathbf{q}j} + k_BT \sum_{\mathbf{q}j} \ln\left[1 - e^{-\hbar\omega_{\mathbf{q}j}/k_BT}\right]
$$

**Entropy** (Eq. 11):
$$
S = \sum_{\mathbf{q}j} k_B \left[\frac{x}{e^x-1} - \ln(1-e^{-x})\right], \quad x = \hbar\omega_{\mathbf{q}j}/k_BT
$$

### 2.5 Gruneisen Parameter (Eq. 15)
$$
\gamma_{\mathbf{q}j}(V) = -\frac{V}{\omega_{\mathbf{q}j}(V)} \frac{\partial\omega_{\mathbf{q}j}(V)}{\partial V}
$$

### 2.6 Quasi-Harmonic Approximation (Eq. 17)
$$
G(T,p) = \min_V \left[F(T;V) + pV\right]
$$
where $F(T;V) \approx U_{el}(V) + F_{ph}(T;V)$.

### 2.7 Heat Capacity at Constant Pressure (Eq. 18)
$$
C_P(T,p) = C_V(T,V(T,p)) + T \frac{\partial V}{\partial T} \frac{\partial S}{\partial V}\bigg|_{V=V(T,p)}
$$

### 2.8 Vinet Equation of State
$$
E(V) = E_0 + \frac{4V_0 B_0}{\eta^2}\left[1 - (1+\eta(1-x))e^{-\eta(1-x)}\right]
$$
where $x = (V/V_0)^{1/3}$, $\eta = \frac{3}{2}(B_0'-1)$.

---

## 3. Main Methodology

### 3.1 Born-von Karman Force Constant Model
Since the original paper uses DFT (VASP) to compute force constants, and DFT is not available for reproduction, we use a Born-von Karman (BvK) force constant model for FCC-Al with parameters fitted to match experimental phonon dispersion data from neutron scattering (Stedman & Nilsson, Gilat & Nicklow).

The model uses 4 neighbor shells with 8 free parameters ($\alpha_n, \beta_n$ for $n=1,...,4$), fitted to minimize the deviation from experimental phonon frequencies at high-symmetry points (X, L, W, K and intermediate points).

### 3.2 Phonon Calculation Procedure
1. Generate FCC neighbor shells up to 4th nearest neighbors
2. Fit force constants to experimental Al phonon dispersion
3. Construct dynamical matrix $D(\mathbf{q})$ at each wave vector
4. Diagonalize to obtain phonon frequencies $\omega_{\mathbf{q}j}$
5. Compute band structure along high-symmetry path L-$\Gamma$-X-W
6. Compute DOS by Gaussian smearing on uniform q-grid

### 3.3 Thermal Properties
1. From phonon DOS, compute $C_V$, $S$, $F$ using Eqs. 9-11
2. For $C_P$, use Gruneisen relation: $C_P - C_V = \gamma^2 C_V^2 T / (V B)$

### 3.4 Quasi-Harmonic Approximation
1. At 10 volumes around equilibrium, scale force constants: $\Phi(a) = \Phi(a_0)(a_0/a)^p$
2. Compute phonon DOS at each volume
3. Add electronic energy from Vinet EOS: $U_{el}(V)$
4. At each temperature, minimize $F_{total}(V) = U_{el}(V) + F_{ph}(T,V)$
5. Extract equilibrium $V(T)$, thermal expansion $\beta(T)$

### 3.5 Key Numerical Details
- **q-grid for DOS**: 30x30x30 uniform grid in BZ (35x35x35 for higher accuracy)
- **Gaussian smearing width**: 0.12-0.15 THz
- **DOS normalization**: Integral over frequency = 3 (3 modes per atom)
- **Force constant fitting**: Nelder-Mead optimization with ~100,000 iterations
- **Force constant power law**: $p \approx 7$ for volume scaling
- **Overflow protection**: Clip $x = \hbar\omega/k_BT$ to max 500 in thermal property calculations

---

## 4. Input Parameters

### 4.1 Physical Parameters for FCC-Al
| Parameter | Symbol | Value | Description |
|-----------|--------|-------|-------------|
| Lattice constant | $a_0$ | 4.05 A | Equilibrium lattice constant |
| Atomic mass | $M$ | 26.9815385 amu | Aluminum mass |
| Bulk modulus | $B_0$ | 79 GPa | Equilibrium bulk modulus |
| $B_0'$ | $B_0'$ | 4.4 | Pressure derivative of bulk modulus |
| Gruneisen parameter | $\gamma$ | 2.2 | Macroscopic Gruneisen parameter |

### 4.2 FCC Structure
| Property | Value |
|----------|-------|
| Space group | Fm-3m (225) |
| Atoms per primitive cell | 1 |
| Phonon branches | 3 (all acoustic) |
| Conventional cell volume | $a^3 = 66.43$ A$^3$ |
| Volume per atom | $a^3/4 = 16.61$ A$^3$ |

### 4.3 Neighbor Shells in FCC
| Shell | Neighbors | Distance | Positions |
|-------|-----------|----------|-----------|
| 1st | 12 | $a/\sqrt{2}$ | $(a/2)(\pm1,\pm1,0)$ and permutations |
| 2nd | 6 | $a$ | $(a)(\pm1,0,0)$ and permutations |
| 3rd | 24 | $a\sqrt{3/2}$ | $(a/2)(\pm2,\pm1,\pm1)$ and permutations |
| 4th | 12 | $a\sqrt{2}$ | $(a/2)(\pm2,\pm2,0)$ and permutations |

### 4.4 Target Phonon Frequencies for Fitting (THz)
| Point | Coordinates (2pi/a) | $\omega_1$ | $\omega_2$ | $\omega_3$ |
|-------|---------------------|------------|------------|------------|
| X | (1,0,0) | 5.82 | 5.82 | 9.69 |
| L | (1/2,1/2,1/2) | 4.18 | 4.18 | 9.64 |
| W | (1,1/2,0) | 5.41 | 7.82 | 9.31 |
| K | (3/4,3/4,0) | 4.85 | 6.50 | 8.15 |

### 4.5 Parameter Variations per Figure
| Figure | Description | Parameters |
|--------|-------------|------------|
| Fig 1 | Phonon band structure and DOS | Fixed: equilibrium $a_0$ |
| Fig 2 | Thermal properties | Temperature: 0-800 K |
| Fig 3(a) | Frequencies vs volume | 10 volumes: 59-79 A$^3$ (conventional cell) |
| Fig 3(b) | Free energy vs volume | 9 temperatures: 0-800 K, 100 K step |
| Fig 3(c) | Thermal expansion | Temperature: 0-800 K |

---

## 5. Banned Libraries

| Banned | Reason |
|--------|--------|
| `phonopy`, `ase`, `pymatgen` | Must implement phonon model from scratch |
| `spglib` | Crystal symmetry not needed for this monatomic model |
| GPU libraries (cupy, pytorch) | Keep in pure NumPy/SciPy |

**Allowed**: `numpy`, `scipy` (constants, optimize, special), `matplotlib`

---

## 6. Data Files to Reproduce

### 6.1 Figure 1: Phonon Band Structure and DOS
**Files**: `data/fig1_bandstructure.csv`, `data/fig1_dos.csv`

| Property | Value |
|----------|-------|
| Description | Phonon band structure along L-Gamma-X-W and phonon DOS |
| X-axis (band) | Wave vector distance (arbitrary units) |
| Y-axis | Frequency (THz), range 0-10 |
| Band points | 80 per segment, 3 segments = 240 points |
| DOS bins | 300 bins, Gaussian smearing 0.12 THz |

**Band columns**: `distance, branch_1_THz, branch_2_THz, branch_3_THz`
**DOS columns**: `frequency_THz, dos_states_per_THz_per_atom`

### 6.2 Figure 2: Thermal Properties
**File**: `data/fig2_thermal.csv`

| Property | Value |
|----------|-------|
| Description | Cv, S, F, Cp vs temperature for FCC-Al |
| X-axis | Temperature (K), range 0-800 |
| Y-axis | Shared axis: J/K/mol for Cv,S,Cp; kJ/mol for F |
| Points | 801 (1 K step) |
| Experimental data | Cp from JANAF tables (Chase 1998) |

**Columns**: `T_K, Cv_J_per_K_per_mol, S_J_per_K_per_mol, F_kJ_per_mol, Cp_J_per_K_per_mol`

### 6.3 Figure 3: QHA Properties
**Files**: `data/fig3a_freq_vs_volume.csv`, `data/fig3c_thermal_expansion.csv`

| Property | Value |
|----------|-------|
| Description | (a) Phonon freq vs volume, (b) F(V,T), (c) thermal expansion |
| Panel (a) X-axis | Volume (A^3), range 59-79 |
| Panel (a) Y-axis | Frequency (THz) at X and L points |
| Panel (c) X-axis | Temperature (K), 0-800 |
| Panel (c) Y-axis | Thermal expansion coefficient (10^-6 K^-1) |
| Volume points | 10 |
| Temperature points | 200 (panel c) |

**Panel (a) columns**: `volume_A3, freq_X_1_THz, freq_X_2_THz, freq_X_3_THz, freq_L_1_THz, freq_L_2_THz, freq_L_3_THz`
**Panel (c) columns**: `T_K, V_eq_A3, beta_K_inv`

---

## 7. Code Structure

### 7.1 Core Module
**Path**: `reproduction/phonon_model.py`

Contains:
- `generate_fcc_shells()`: Generate neighbor shells for FCC lattice
- `dynamical_matrix()`: Compute 3x3 dynamical matrix at wave vector q
- `phonon_frequencies_Hz()`: Compute phonon frequencies
- `fit_al_force_constants()`: Fit BvK force constants to experimental data
- `compute_bandstructure()`: Compute band structure along high-symmetry path
- `compute_dos()`: Compute phonon DOS via Gaussian smearing
- `thermal_properties()`: Compute Cv, S, F from phonon DOS
- `phonon_free_energy_per_atom()`: Phonon Helmholtz free energy per atom
- `vinet_energy()`: Vinet equation of state
- `scale_force_constants()`: Volume-dependent force constants for QHA
- `qha_properties()`: Full QHA calculation

### 7.2 Figure Scripts
| Figure | Script Path |
|--------|-------------|
| Fig 1 | `reproduction/fig1_bandstructure_dos.py` |
| Fig 2 | `reproduction/fig2_thermal_properties.py` |
| Fig 3 | `reproduction/fig3_qha_properties.py` |

### 7.3 Output Directories
- Data files: `data/`
- Plot images: `plots/`

---

## 8. Reproduction Requirements

### 8.1 Computational Resources
| Requirement | Value | Notes |
|-------------|-------|-------|
| **Total time limit** | 1 hour | Fig 3 is slowest (QHA at 10 volumes) |
| **Memory limit** | 2 GB | Modest requirements |
| **CPU** | Single-threaded | No parallelization needed |

### 8.2 Numerical Precision
- Export data with 8 decimal places
- Force constant fitting: Nelder-Mead with tolerances 1e-8 (x) and 1e-10 (f)

### 8.3 Environment
```bash
pip install numpy scipy matplotlib
```

### 8.4 Execution
Run scripts from the `reproduction/` directory:
```bash
cd reproduction/
python3 fig1_bandstructure_dos.py
python3 fig2_thermal_properties.py
python3 fig3_qha_properties.py
```
