# Verification Runs

## Run 1: bandaru2018-kink-resistivity

- **Blueprint**: `bandaru2018-kink-resistivity`
- **Job ID**: `933487a7b5162db2`
- **Status**: SUCCESS
- **Container**: `docker://git.pku.edu.cn/2200011363/knowledge-distiller:latest`
- **Description**: Kink mode growth rate scaling (Fig 4b) and resistivity profile (Fig 5)
- **Parameters**: f_re_1=0.0, f_re_2=0.5, f_re_3=1.0 (default)
- **Outputs**: `fig4b_kink_scaling.png`, `fig4b_kink_scaling.csv`, `fig5_resistivity.png`, `fig5_resistivity.csv`
- **Notes**: Successfully computed kink mode growth rate gamma*tau_A vs S^{-1} for three RE fractions,
  and generated the stepped resistivity profile eta/eta_axis vs psi_N.
