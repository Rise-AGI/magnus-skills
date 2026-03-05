# Bandaru et al. (2018) — RE-MHD Fluid Model

## Citation

Bandaru, V., Hoelzl, M., Artola, F.J., Papp, G., & Huijsmans, G.T.A. (2018).
"Simulating the non-linear interaction of relativistic electrons and tokamak plasma instabilities:
implementation and validation of a fluid model." *Physical Review E*, submitted / arXiv:1809.03489.

## Domain

Magnetohydrodynamics (MHD), tokamak disruptions, runaway electrons (RE), plasma stability.

## Summary

This paper presents a self-consistent fluid model for runaway electrons (REs) coupled with
resistive MHD in the JOREK code. The model treats REs as a separate fluid species with its
own density and current, evolving under Dreicer generation (Connor-Hastie formula) and
avalanche multiplication (Rosenbluth-Putvinski). Key results include:

1. **Current quench dynamics** (Fig 3a): During a thermal quench, the rising resistivity drives
   an inductive electric field that generates REs via Dreicer and avalanche mechanisms. The RE
   current grows to partially replace the decaying thermal current.

2. **Current density profiles** (Fig 3b): Post-quench RE current is centrally peaked (due to
   on-axis E-field concentration), while thermal current is broader and reduced.

3. **Kink mode growth rate scaling** (Fig 4b): The internal kink mode growth rate gamma*tau_A
   scales as S^{-1/3} with the inverse Lundquist number. RE current fraction reduces the
   effective resistivity, shifting the scaling curve.

4. **Resistivity profile** (Fig 5): The ITER VDE simulation uses a stepped resistivity profile:
   flat in core (psi_N < 0.8), transitioning to ~3x at the LCFS, exponentially rising in vacuum.

5. **VDE current evolution** (Fig 6a): During an ITER vertical displacement event, RE avalanche
   growth can maintain or increase total plasma current even as thermal current decays resistively.

## Key Equations

- **Dreicer source** (Eq. 7): S_p = C * n_e * nu_ee * (E/E_D)^alpha * exp(-1/(4*E/E_D) - sqrt((1+Z)/(E/E_D)))
- **Avalanche source** (Eq. 8): S_s = n_r * nu_fp * (E/E_c - 1) / ln(Lambda) * sqrt(pi*phi/(3*(Z+5)))
- **Critical field**: E_c = n_e * e^3 * ln(Lambda) / (4*pi*eps0^2 * m_e * c^2)
- **Dreicer field**: E_D = n_e * e^3 * ln(Lambda) / (4*pi*eps0^2 * T_e)
- **Spitzer resistivity**: eta = eta0 * (T/T0)^{-3/2}

## Parameters

### GO-code benchmark (Section 3.2, Figs 3a-3b)
- R = 10 m, a = 1 m, B_phi = 1 T
- T0 = 1.7 keV, n0 = 1e20 m^-3, Z_eff = 1, eta0 = 1.1e-7 Ohm*m
- I_p = 0.67 MA, ln(Lambda) = 15

### Kink scaling (Section 3.3, Fig 4b)
- S^{-1} range: 1e-7 to 1e-2
- RE fractions: I_r/I_p = 0, 0.5, 1.0

### ITER VDE (Section 4, Figs 5-6)
- I_p = 14.5 MA, B_phi = 4.8 T, n_e = 5e19 m^-3
- eta_axis = 1.24e-4 Ohm*m, T = 2.35 eV, a = 2 m, R = 6.2 m

## Reproduction Code

The `reproduction/` directory contains:
- `re_fluid.py` — Core physics: PlasmaParams, Dreicer/avalanche sources, CurrentQuench0D, VDE1D, kink scaling
- `fig3a_current_quench.py` — Current quench time evolution
- `fig3b_profiles.py` — Current density profiles pre/post quench
- `fig4b_kink_scaling.py` — Kink mode growth rate vs S^{-1}
- `fig5_resistivity.py` — Resistivity profile for ITER VDE
- `fig6a_vde_currents.py` — VDE current evolution with/without REs

## Blueprints

- `bandaru2018-kink-resistivity`: Kink scaling (Fig 4b) + resistivity profile (Fig 5)
- `bandaru2018-current-quench`: Current quench evolution (Fig 3a) + profiles (Fig 3b)
- `bandaru2018-vde-currents`: ITER VDE current evolution (Fig 6a)
