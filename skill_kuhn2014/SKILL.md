# Schwinger Model Quantum Simulation - Kuhn, Cirac, Banuls (2014)

## Citation

S. Kuhn, J. I. Cirac, and M.-C. Banuls, "Quantum simulation of the Schwinger model: A study of feasibility," Phys. Rev. A **90**, 042305 (2014). [arXiv:1407.4995](https://arxiv.org/abs/1407.4995)

## Overview

This skill covers the computational reproduction of key results from the study of quantum simulation feasibility for the Schwinger model (QED in 1+1 dimensions). The paper investigates three questions:

1. **Finite-dimensional truncation**: How well do finite-dimensional link Hilbert spaces approximate the full Schwinger model?
2. **Adiabatic preparation**: What time scales are needed to prepare the interacting vacuum state?
3. **Noise robustness**: How does gauge invariance breaking noise affect the results?

## Models

### Truncated compact QED (cQED)
- Link operators truncated at dimension $d = 2J+1$
- Preserves the same U(1) gauge symmetry as the full Schwinger model
- $L^+|J\rangle = 0$ (hard boundary)

### $Z_d$ model
- Cyclic link operators: $L^+|J\rangle = |-J\rangle$
- Has discrete $Z_d$ gauge symmetry (different from Schwinger)
- Converges faster: already very accurate at $d = 3$

## Key Equations

| Equation | Description |
|----------|-------------|
| $W = \sum_n (L_n^z)^2 + x \sum_n (\sigma_n^+ L_n^+ \sigma_{n+1}^- + \text{h.c.})$ | Rescaled Hamiltonian (massless) |
| $x = 1/(ag)^2$ | Dimensionless lattice spacing parameter |
| $G_n = L_n^z - L_{n-1}^z - \phi_n^\dagger\phi_n + \frac{1}{2}[1-(-1)^n]$ | Gauss law generator |
| $x(t) = x_F (t/T)^3$ | Cubic adiabatic ramp |

## Key Results

- **Convergence with d**: The Zd model at $d=3$ is already extremely close to the full model. The cQED model converges more slowly but reaches good accuracy by $d=9$.
- **Adiabatic time**: Overlaps > 0.99 achieved by $T \approx 60$, remarkably independent of system size N and link dimension d.
- **Noise tolerance**: Energy remains reliable up to $\lambda \approx 10^{-3}$. Penalty energy scales as $\lambda^2 t^8 / T^6$, consistent with second-order perturbation theory.

## Computational Methods

### Exact Diagonalization (our reproduction)
- Works for $N \leq 12$ with $d \leq 9$
- Gauss law reduces Hilbert space from $2^N \cdot d^{N-1}$ to $\binom{N}{N/2}$ (approx)
- Sparse eigensolvers for ground states; dense matrix exponentials for time evolution

### MPS (original paper)
- Bond dimension $D = 50-200$
- System sizes $N = 50-200$
- Suzuki-Trotter time evolution

## Parameter Ranges

| Parameter | Typical Range | Notes |
|-----------|--------------|-------|
| $N$ | 4 - 200 | System size (even) |
| $d$ | 3 - 9 | Link dimension (odd) |
| $x$ | 0 - 100 | Coupling parameter |
| $T$ | 50 - 100 | Evolution time |
| $\lambda$ | $10^{-4}$ - $10^{-3}$ | Noise strength |

## Pitfalls

1. **Gauss law constraint**: Must enforce the physical subspace. Without it, spurious states dominate.
2. **Zero charge sector**: Restrict to total charge = 0 for comparison with analytical results.
3. **Open boundary conditions**: The paper uses OBC, not periodic. This affects the Gauss law at boundaries.
4. **cQED vs Zd normalization**: The $L^+$ operator in cQED includes a factor $i/\sqrt{l(l+1)}$ from the Schwinger boson representation. The simplified form used here drops this factor, which is consistent with the Zd-like treatment but may affect quantitative comparison.
5. **Continuum limit**: The exact energy density $-1/\pi$ requires both $N \to \infty$ and $x \to \infty$ extrapolation.

## Repository

Code: [Rise-AGI/magnus-skills](https://github.com/Rise-AGI/magnus-skills) branch `kuhn2014`

## Blueprints

- [kuhn2014-ground-state](magnus:///blueprints/kuhn2014-ground-state): Ground state energy density and gap (Figs 1-2)
- [kuhn2014-adiabatic](magnus:///blueprints/kuhn2014-adiabatic): Adiabatic preparation (Figs 3-4)
- [kuhn2014-noise](magnus:///blueprints/kuhn2014-noise): Noise effects (Figs 5-6)
