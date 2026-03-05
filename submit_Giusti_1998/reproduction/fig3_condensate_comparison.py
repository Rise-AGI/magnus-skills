"""
Figure 3: Chiral condensate determination from lattice data.

Reproduces the computation of the physical chiral condensate
from the lattice input data of Tables 1-4, using both methods:
- Method 1: Eq. 41 (Wilson action only, uses Z_S)
- Method 2: Eq. 42 (Both actions, uses Z_P/Z_A)

Shows comparison across beta values and actions.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, '.')
from qcd_perturbative import (
    compute_condensate_method1, compute_condensate_method2,
    delta_z_ri_msbar, evolution_coefficient, run_condensate,
    LATTICE_RUNS, RENORM_CONSTANTS, CONDENSATE_PHYSICAL,
    F_CHI, LAMBDA_QCD_MSBAR
)

plt.rcParams.update({
    'font.size': 12,
    'axes.labelsize': 14,
    'axes.titlesize': 14,
    'legend.fontsize': 9,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'text.usetex': False,
    'mathtext.fontset': 'cm',
})


def main():
    print("=" * 60)
    print("Figure 3: Chiral condensate from lattice data")
    print("=" * 60)

    mu_target = 2.0  # GeV

    # Compute condensate for all runs
    results = []
    print(f"\nChiral condensate in MSbar at mu = {mu_target} GeV:")
    print(f"{'Run':<8} {'Method':>8} {'Computed':>12} {'Paper':>12} {'Action':>8} {'beta':>6}")
    print("-" * 60)

    for run_key in ['C60a', 'C60b', 'W60', 'C62', 'W62', 'C64', 'W64']:
        run = LATTICE_RUNS[run_key]
        paper_data = CONDENSATE_PHYSICAL[run_key]

        # Method 2 (available for all)
        cond2 = compute_condensate_method2(run_key, mu_target)
        paper_val2 = paper_data['psi2'][0] if paper_data['psi2'][0] else None
        results.append({
            'run': run_key, 'method': 2, 'computed': cond2,
            'paper': paper_val2,
            'paper_stat_err': paper_data['psi2'][1] if paper_val2 else None,
            'paper_sys_err': paper_data['psi2'][2] if paper_val2 else None,
            'action': run['action'], 'beta': run['beta']
        })
        paper_str2 = f"{paper_val2:.4f}" if paper_val2 else "N/A"
        print(f"{run_key:<8} {'psi2':>8} {cond2:>12.4f} {paper_str2:>12} {run['action']:>8} {run['beta']:>6.1f}")

        # Method 1 (Wilson only)
        if run['action'] == 'Wilson':
            cond1 = compute_condensate_method1(run_key, mu_target)
            paper_val1 = paper_data['psi1'][0] if paper_data['psi1'][0] else None
            results.append({
                'run': run_key, 'method': 1, 'computed': cond1,
                'paper': paper_val1,
                'paper_stat_err': paper_data['psi1'][1] if paper_val1 else None,
                'paper_sys_err': paper_data['psi1'][2] if paper_val1 else None,
                'action': run['action'], 'beta': run['beta']
            })
            paper_str1 = f"{paper_val1:.4f}" if paper_val1 else "N/A"
            print(f"{run_key:<8} {'psi1':>8} {cond1:>12.4f} {paper_str1:>12} {run['action']:>8} {run['beta']:>6.1f}")

    # Final result from paper
    print(f"\n{'='*60}")
    print("Paper's best estimate (C62, method 2):")
    print(f"  <psibar psi>^MSbar(2 GeV) = -0.0147(8)(16)(12) GeV^3")
    print(f"  = -[245(4)(9)(7) MeV]^3")

    best = compute_condensate_method2('C62', 2.0)
    print(f"  Our computation: {best:.4f} GeV^3")

    # Run to 1 GeV
    cond_1gev = run_condensate(best, 2.0, 1.0)
    print(f"\n  RG run to 1 GeV: {cond_1gev:.4f} GeV^3")
    print(f"  Paper value at 1 GeV: 0.0124 GeV^3")

    # Plot: comparison across runs
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # Left panel: Method 2 results for all runs
    beta_labels = []
    computed_vals = []
    paper_vals = []
    paper_errs = []
    run_labels = []

    for r in results:
        if r['method'] == 2 and r['paper'] is not None:
            run_labels.append(r['run'])
            computed_vals.append(r['computed'])
            paper_vals.append(r['paper'])
            total_err = np.sqrt(r['paper_stat_err']**2 + r['paper_sys_err']**2)
            paper_errs.append(total_err)

    x = np.arange(len(run_labels))
    width = 0.35

    ax1.bar(x - width/2, computed_vals, width, label='Computed (this work)',
            color='steelblue', alpha=0.8)
    ax1.bar(x + width/2, paper_vals, width, label='Paper (Table 4)',
            color='coral', alpha=0.8, yerr=paper_errs, capsize=4)

    ax1.set_xlabel('Run')
    ax1.set_ylabel(r'$-\langle\bar{\psi}\psi\rangle^{\overline{MS}}(2\,\mathrm{GeV})$ [GeV$^3$]')
    ax1.set_title(r'Method 2: $\langle\bar{\psi}\psi\rangle_2$ comparison')
    ax1.set_xticks(x)
    ax1.set_xticklabels(run_labels, rotation=45)
    ax1.legend()
    ax1.grid(True, linestyle='--', alpha=0.3, axis='y')

    # Right panel: beta dependence for Clover and Wilson
    for action, color, marker in [('Clover', 'blue', 'o'), ('Wilson', 'red', 's')]:
        betas = []
        vals = []
        for r in results:
            if r['method'] == 2 and r['action'] == action and r['paper'] is not None:
                # Skip duplicates (C60a, C60b)
                if r['run'] == 'C60b':
                    continue
                betas.append(r['beta'])
                vals.append(r['computed'])
        if betas:
            ax2.plot(betas, vals, marker=marker, color=color, linewidth=1.5,
                    markersize=8, label=f'{action} (computed)')

    # Paper's final result band
    ax2.axhspan(0.0147 - 0.0016, 0.0147 + 0.0016, alpha=0.2, color='green',
                label='Paper final result')
    ax2.axhline(y=0.0147, color='green', linestyle='--', alpha=0.5)

    ax2.set_xlabel(r'$\beta = 6/g_0^2$')
    ax2.set_ylabel(r'$-\langle\bar{\psi}\psi\rangle^{\overline{MS}}(2\,\mathrm{GeV})$ [GeV$^3$]')
    ax2.set_title(r'$\beta$ dependence of chiral condensate')
    ax2.legend(loc='best')
    ax2.grid(True, linestyle='--', alpha=0.3)

    plt.tight_layout()
    plt.savefig('../plots/fig3_condensate_comparison.png', bbox_inches='tight')
    print("\nFigure saved to plots/fig3_condensate_comparison.png")
    plt.close()

    # Export data
    with open('../data/fig3_condensate_comparison.csv', 'w') as f:
        f.write("# Figure 3: Chiral condensate comparison\n")
        f.write("# All values are -<psibar psi>^MSbar at mu=2 GeV in GeV^3\n")
        f.write("run,method,action,beta,computed,paper_value,paper_stat_err,paper_sys_err\n")
        for r in results:
            pv = f"{r['paper']:.8f}" if r['paper'] else ""
            pse = f"{r['paper_stat_err']:.8f}" if r['paper_stat_err'] else ""
            psys = f"{r['paper_sys_err']:.8f}" if r['paper_sys_err'] else ""
            f.write(f"{r['run']},{r['method']},{r['action']},{r['beta']:.1f},"
                    f"{r['computed']:.8f},{pv},{pse},{psys}\n")
    print("Data exported to data/fig3_condensate_comparison.csv")


main()
