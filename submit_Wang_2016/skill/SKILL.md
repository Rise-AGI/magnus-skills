# Runaway Electron Energy Limits in Tokamak Fields

## Domain

Relativistic plasma physics, magnetic confinement fusion, runaway electron dynamics.

## Source

Wang, Qin & Liu, "Multi-scale full-orbit analysis on phase-space behavior of runaway electrons in tokamak fields with synchrotron radiation," Phys. Plasmas **23**, 062505 (2016). DOI: 10.1063/1.4953608

## Physical Model

Runaway electrons in tokamaks are accelerated by a loop electric field E_l and decelerated by synchrotron and curvature radiation losses. The energy limit is set by the balance between these.

### Key Equations

**Relaxation (gyro-center) equations:**
```
dp_∥/dt = eE_l − D·p_∥
dp_⊥/dt =      − D·p_⊥
```
where D = P_rad · γ · m_e / p² is the radiation drag coefficient.

**Synchrotron radiation power:**
```
P_sync = e⁴B₀²p_⊥² / (6πε₀m_e⁴c³)
```

**Curvature radiation power:**
```
P_curv = e²γ⁴v_∥⁴ / (6πε₀c³R_eff²)
```

**Analytical energy limit (curvature-dominated equilibrium):**
```
γ_max⁴ = 6πε₀ · E_l · R_eff² / e
```

For default EAST-like parameters (E_l=0.2 V/m, B₀=2T, R₀=1.7m, q=2): γ_max ≈ 157, E_max ≈ 80 MeV.

### Implementation Notes

- All ODE integration uses normalized momentum p̃ = p/(m₀c) to avoid overflow
- The normalized drag coefficients are:
  - τ_s = C_SYNC · m_e · B₀² (synchrotron)
  - τ_c = e²/(6πε₀·c·m_e·R_eff²) (curvature)
- Effective curvature radius accounts for toroidal + helical field-line geometry
- The relaxation model does NOT capture the full-orbit neoclassical pitch-angle scattering (the paper's key contribution), which requires ~10¹² time steps

## Available Blueprints

| Blueprint ID | Description |
|---|---|
| `wang2016-analytical` | Analytical energy limit from curvature balance (instant) |
| `wang2016-momentum` | Full momentum evolution trajectory (Fig. 4) |
| `wang2016-energy-scan` | Parametric energy limit scan vs E_l, R₀, or q (Figs. 9/10/12) |

## Parameter Ranges

| Parameter | Symbol | Typical Range | Default |
|---|---|---|---|
| Loop electric field | E_l | 0.05–5.0 V/m | 0.2 V/m |
| Toroidal field | B₀ | 1.5–5.0 T | 2.0 T |
| Major radius | R₀ | 0.5–5.0 m | 1.7 m |
| Safety factor | q | 0.5–6.0 | 2.0 |

## Key Results

- Energy limit increases with E_l (as E_l^{1/4}) and R₀ (as R₀^{1/2})
- Energy limit is insensitive to B₀ in the curvature-dominated regime
- Safety factor affects energy through field-line helicity correction to curvature
- Synchrotron radiation controls p_⊥ decay; curvature radiation sets the p_∥ limit
