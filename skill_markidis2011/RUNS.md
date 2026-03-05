# Run History

## Blueprints

| Blueprint ID | Title | Verified |
|---|---|---|
| [markidis2011-twostream](magnus:///blueprints/markidis2011-twostream) | EC-PIC: Two-Stream Instability | Yes |
| [markidis2011-finitegrid](magnus:///blueprints/markidis2011-finitegrid) | EC-PIC: Finite Grid Test | No (not tested) |
| [markidis2011-weibel](magnus:///blueprints/markidis2011-weibel) | EM-PIC: Weibel Instability | Yes |

## Verification Runs

| Date | Blueprint | Job | Status | Notes |
|---|---|---|---|---|
| 2026-03-06 | markidis2011-weibel | [Job 89faad30](magnus:///jobs/89faad3085384fe5) | Success | gamma=0.23 vs theory 0.22 wpe |
| 2026-03-06 | markidis2011-twostream | [Job 01cf1ee0](magnus:///jobs/01cf1ee0efdb6824) | Success | EC-PIC dE=5e-4%, Explicit dE=199% |

## Local Verification Results

All scripts ran successfully on local environment:

| Script | Status | Key Result |
|---|---|---|
| dispersion_cold_plasma.py | Success | EC-PIC dispersion always stable |
| fig5_tolerance_energy.py | Success | Smaller tol = better energy conservation |
| fig7_finite_grid_energy.py | Success | EC-PIC: 3e-6% vs Explicit: 95% energy change |
| fig8_twostream_growth.py | Success | Growth rate ~0.42 wpe (theory: 0.35) |
| fig9_twostream_energy.py | Success | EC-PIC: 5e-4% vs Explicit: 199% |
| fig10_11_weibel.py | Success | Weibel Bz growth observed |
