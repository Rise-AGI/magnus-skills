# Markos2006: Numerical Analysis of the Anderson Localization

## Citation

P. Markos, "Numerical Analysis of the Anderson Localization," arXiv:cond-mat/0609580, submitted to Acta Physica Slovaca (2006).

## Overview

This paper provides a comprehensive numerical tutorial on Anderson localization in disordered electron systems. It covers:
- The Anderson tight-binding model on d-dimensional lattices
- Metal-insulator transition and mobility edges
- Transfer matrix methods for conductance
- Level spacing statistics (Wigner vs Poisson)
- Weak localization corrections
- Critical regime and scaling theory
- Universal conductance fluctuations

## Key Models

### Anderson Hamiltonian (Eq. 6)
$$H = W \sum_r \varepsilon_r c_r^\dagger c_r + V \sum_{[rr']} c_r^\dagger c_{r'}$$

- Box disorder: $\varepsilon_r \in [-W/2, W/2]$ uniformly
- Gaussian disorder: $\varepsilon_r \sim N(0, W^2)$
- Hopping $V = 1$ (energy unit)

### Critical Disorder Values (3D)
- **Box disorder**: $W_c \approx 16.5$ (at band center E=0)
- **Gaussian disorder**: $W_c \approx 6.1$ (at band center E=0)

### 1D Localization
- All states localized for any disorder
- Localization length $\lambda = 1/\text{Re}(\gamma)$
- For E=1, W=1 (box): $\lambda \approx 72$ (analytical), ~70 (numerical)
- Transfer matrix: $T_n = \begin{pmatrix} E-\varepsilon_n & -1 \\ 1 & 0 \end{pmatrix}$
- Conductance: $g = 1/\cosh^2(x/2)$

### Level Spacing Statistics
- **Metallic regime**: Wigner surmise $p_1(s) = (\pi/2) s \exp(-\pi s^2/4)$ (GOE)
- **Localized regime**: Poisson $p(s) = \exp(-s)$
- **Critical regime**: universal distribution between Wigner and Poisson

### Weak Localization (2D)
- Logarithmic correction: $\delta g = -(1/\pi) \ln(L/\ell)$
- Slope $\approx -1/\pi = -0.318$ for orthogonal systems

## Parameter Ranges

| Parameter | Symbol | Typical Values |
|-----------|--------|----------------|
| Lattice size | L | 6, 8, 10, 14, 22 |
| Box disorder | W | 0 - 32 |
| Gaussian disorder | W | 0 - 8 |
| Energy | E | Band center (0) or E=1 |
| Realizations | N_stat | 200 - 10^6 |

## Computational Methods

1. **Exact diagonalization**: For L^d lattice (feasible up to L~20 in 3D)
2. **Transfer matrix**: For quasi-1D systems (efficient for 1D, strips)
3. **QR stabilization**: Essential for long transfer matrix products

## Common Pitfalls

- L=14 in 3D means 2744x2744 matrix — ~17s per diagonalization
- Statistical convergence requires many realizations (100-1000+)
- Transfer matrix becomes numerically unstable without QR decomposition
- Level spacing must be normalized: $\langle s \rangle = 1$
- Energy window selection affects level statistics quality

## Blueprints

| Blueprint | Description |
|-----------|-------------|
| [markos2006-dos](magnus:///blueprints/markos2006-dos) | Density of states for 3D Anderson model |
| [markos2006-1d-localization](magnus:///blueprints/markos2006-1d-localization) | 1D chain localization and p(x) distribution |
| [markos2006-level-statistics](magnus:///blueprints/markos2006-level-statistics) | Level spacing distribution in metallic/localized/critical regimes |

## Repository

Code branch: `markos2006` on Rise-AGI/magnus-skills
