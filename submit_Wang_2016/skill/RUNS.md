# Runs Audit Trail

## Run 1: Analytical Blueprint Verification

- **Blueprint**: `wang2016-analytical`
- **Job ID**: `36f5bb8120ae12c9`
- **Parameters**: defaults (E_l=0.2, R₀=1.7, q=2.0, r=0.1)
- **Result**: γ_max = 156.64, E_max = 79.53 MeV
- **Status**: SUCCESS
- **Verification**: Matches analytical formula γ⁴ = 6πε₀·E_l·R_eff²/e

## Run 2: Momentum Evolution Blueprint

- **Blueprint**: `wang2016-momentum`
- **Job ID**: (see latest run)
- **Parameters**: defaults (E_l=0.2, B₀=2.0, R₀=1.7, q=2.0, t_max=3.5)
- **Expected**: E_max ≈ 79.5 MeV, p_∥ saturates at ~157 m₀c
- **Status**: SUBMITTED

## Local Reproduction Runs

### Figure 4: Momentum Evolution
- E_max = 79.50 MeV, max p_∥ = 156.58 m₀c
- Momentum trajectory shows p_∥ increasing while p_⊥ decays to near zero
- Energy saturates at curvature-radiation limit

### Figure 9: Energy vs Electric Field
- B₀=1.5T: E_max range [55.6, 162.2] MeV across E_l = 0.05–5.0 V/m
- B₀=2.0T: E_max range [55.9, 162.2] MeV
- Energy limit scales as E_l^{1/4}, weak B₀ dependence confirms curvature dominance

### Figure 10: Energy vs Major Radius
- E_max range [36.0, 117.1] MeV across R₀ = 0.5–5.0 m
- Energy limit scales as R₀^{1/2}

### Figure 12: Energy vs Safety Factor
- E_max range [78.5, 79.8] MeV across q = 0.5–6.0
- Very weak q dependence (relaxation model); full-orbit effects are stronger
