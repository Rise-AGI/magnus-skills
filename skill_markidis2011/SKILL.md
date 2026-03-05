# Energy-Conserving Particle-in-Cell Method

## Citation
Markidis, S. & Lapenta, G., "The Energy Conserving Particle-in-Cell Method",
Journal of Computational Physics 230 (2011) 7037-7052.
DOI: 10.1016/j.jcp.2011.05.033

## Repository
Branch: [markidis2011](magnus:///blueprints/markidis2011-twostream) on Rise-AGI/magnus-skills

## Overview
This paper presents a Particle-in-Cell (PIC) method that exactly conserves total energy.
The key innovation is using the **implicit midpoint rule** for both particle equations and Maxwell's
equations, solved concurrently by a **Jacobian-Free Newton-Krylov (JFNK)** solver.

## Key Equations

### Midpoint Rule for Particles
```
v_bar = v_n + (q/2m) * dt * E_bar(x_bar)
x_{n+1} = x_n + v_bar * dt
v_{n+1} = 2*v_bar - v_n
```
where bars denote time-averaged quantities at (n+1/2).

### Current Density (Critical for Energy Conservation)
```
J_bar_g = sum_p q_s * v_bar_p * W(x_g - x_bar_p) / V_g
```
The current must use the **average** position and velocity — this is what ensures energy conservation.

### Electrostatic Ampere's Law
```
E_{n+1} - E_n = -4*pi * J_bar * dt
```

### Cold Plasma Dispersion (EC-PIC)
```
tan(omega*dt/2) * sin(omega*dt/2) = (omega_pe*dt/2)^2
```
Always has real roots -> **linearly unconditionally stable** (no CFL limit on dt).

### Stability Condition (Momentum Non-Conservation)
```
v_c * dt / dx < 1.5
```
Required to avoid aliasing instability from particle self-forces.

## Key Physical Tests

### 1. Two-Stream Instability (Electrostatic)
- Two cold electron beams at +/-0.2c
- Growth rate: gamma = 0.35355 * omega_pe for k=1 mode
- EC-PIC: energy variation ~1e-4%
- Explicit PIC: energy variation ~5%

### 2. Finite Grid Instability (Electrostatic)
- Maxwellian plasma with dx >> lambda_D
- Explicit PIC: numerical heating ~2%
- EC-PIC: energy variation ~1e-7%

### 3. Weibel Instability (Electromagnetic)
- Bi-Maxwellian with anisotropy a=15 (vthy=0.4c, vthx=0.1c)
- Growth rate: gamma = 0.22 * omega_pe for k=1 Bz mode
- Energy conserved to ~1e-3%

## Parameter Ranges
| Parameter | Typical Range | Notes |
|-----------|--------------|-------|
| omega_pe*dt | 0.1-2.0 | EC-PIC stable for any value |
| dx/lambda_D | 1-12 | EC-PIC stable even for large ratios |
| N_particles | 1000-800000 | Cost linear in N |
| JFNK tolerance | 1e-5 to 1e-10 | Smaller = better energy conservation |
| Newton iterations | 2-5 avg | Convergent in all tested cases |
| Krylov iterations | 2-3 avg | GMRes solver |

## Implementation Pitfalls
1. **JFNK cost**: O(N_particles) per solver iteration - very expensive for many particles
2. **Charge conservation**: The method does NOT conserve charge by default; pseudo-current correction (Marder) can mitigate
3. **Momentum non-conservation**: Causes aliasing instability if v*dt/dx > 1.5
4. **EM-PIC CFL**: If using explicit FDTD for fields (not full EC-PIC), must satisfy c*dt/dx < 1
5. **Energy conservation precision**: Controlled by JFNK tolerance, not by dt or dx

## Normalization
All code uses CGS-like normalized units:
- omega_pe = 1, c = 1, m_e = 1
- Q = WP^2 / (QM * N/L) absorbs 4*pi factor
- Energy units: n_e * m_e * c^2 / 2
