# Run History

## Blueprints

| Blueprint ID | Title | Verified |
|---|---|---|
| kuhn2014-ground-state | Ground State Properties (Figs 1-2) | Yes (local) |
| kuhn2014-adiabatic | Adiabatic Preparation (Figs 3-4) | No (not tested on cluster) |
| kuhn2014-noise | Noise Effects (Figs 5-6) | No (not tested on cluster) |

## Verification Runs

| Date | Blueprint | Job | Status | Notes |
|---|---|---|---|---|
| 2026-03-06 | kuhn2014-ground-state | [Job f707d996306593b6](magnus:///jobs/f707d996306593b6) | Running | Cluster verification in progress |
| 2026-03-06 | kuhn2014-ground-state | local | Success | All 6 figure scripts run successfully on local machine. Energy densities, gaps, adiabatic overlaps, and noise effects all computed correctly. |

## Local Verification Results

### Figure 1: Energy Density
- cQED and Zd models computed for N=4,6,8 and d=3,5,7
- Both models converge with increasing d
- Zd model at d=3 already matches d=5,7 closely (consistent with paper)

### Figure 2: Energy Gap
- Gap increases approximately linearly with x (consistent with paper's observation)
- Nearly independent of d for the Zd model

### Figures 3-4: Adiabatic Preparation
- Overlaps exceed 0.999 for T >= 40 (both models, N=4,6, d=3,5)
- Results independent of N and d (consistent with paper)

### Figures 5-6: Noise Effects
- Penalty energy scales as ~lambda^2 (consistent with perturbative prediction)
- Overlap degrades at lambda ~ 10^{-3} (consistent with paper)
