# Run History

## Blueprints

| Blueprint ID | Title | Verified |
|---|---|---|
| markos2006-dos | Density of States | No (not tested remotely) |
| markos2006-1d-localization | 1D Localization | Yes |
| markos2006-level-statistics | Level Statistics | No (not tested remotely) |

## Verification Runs

| Date | Blueprint | Job | Status | Notes |
|---|---|---|---|---|
| 2026-03-06 | markos2006-1d-localization | [Job 655f5e75f2cd682d](magnus:///jobs/655f5e75f2cd682d) | Success | Lz=100, n_samples=1000: <x>=4.17, var=3.74 — consistent with analytical prediction |

## Local Verification

All 6 figure scripts were tested locally with full parameters:
- fig3_dos.py: Density of states for 3D Anderson model (box and Gaussian disorder) — matches paper Fig 3
- fig13_level_spacing.py: Wigner distribution in metallic, Poisson in localized regime — matches paper Fig 13
- fig14_critical_spacing.py: Critical level spacing distribution — matches paper Fig 14
- fig19_localization_1d.py: <x> linear in Lz, slope gives lambda=69.3 (paper: 70.92, analytical: 72)
- fig20_px_distribution.py: <x>=9.99, var=13.99 (paper: <x>=10.41, var=14.45) — good agreement
- fig25_weak_localization.py: Logarithmic decrease of conductance with system size
