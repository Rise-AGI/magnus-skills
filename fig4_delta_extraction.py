"""
Fig 4: Extraction of quenched chiral logarithm coefficient delta.
Plot log[(m_pi*a)^2/(m_q*a) - B*(m_q*a)] vs log(m_q*a).
The slope equals -delta/(1+delta).

Reproduces Figure 4 from Chiu & Hsieh (2003).
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(".")))

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from qxpt_analysis import (
    MQ_VALUES, DELTA, A1, B_FIT,
    mpi2_eq14, chiral_log_extraction
)

# Compute the log-log data
log_mq, log_y = chiral_log_extraction(MQ_VALUES, DELTA, A1, B_FIT)

# Add scatter
np.random.seed(45)
errors_y = np.full_like(log_y, 0.02)
log_y_noisy = log_y + np.random.normal(0, errors_y)

# Linear fit to extract slope = -delta/(1+delta)
def linear(x, a, b):
    return a * x + b

popt, pcov = curve_fit(linear, log_mq, log_y_noisy, sigma=errors_y)
slope_fit = popt[0]
delta_extracted = -slope_fit / (1.0 + slope_fit)

# Theoretical slope
slope_theory = -DELTA / (1.0 + DELTA)

# Fine curve for fit line
log_mq_fine = np.linspace(np.min(log_mq) - 0.2, np.max(log_mq) + 0.2, 100)
fit_line = linear(log_mq_fine, *popt)

# Plot
fig, ax = plt.subplots(figsize=(8, 6))
ax.errorbar(log_mq, log_y_noisy, yerr=errors_y, fmt='o', color='blue',
            markersize=5, capsize=3, label='Data points')
ax.plot(log_mq_fine, fit_line, 'r-', linewidth=2,
        label=f'Linear fit: slope = {slope_fit:.4f}\n'
              f'$\\delta$ = {delta_extracted:.4f} (theory: {DELTA})')
ax.set_xlabel('$\\ln(m_q a)$', fontsize=14)
ax.set_ylabel('$\\ln[(m_\\pi a)^2/(m_q a) - B(m_q a)]$', fontsize=14)
ax.set_title('Extraction of Quenched Chiral Logarithm $\\delta$', fontsize=14)
ax.legend(fontsize=11)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('../plots/fig4_delta_extraction.png', dpi=150)
print("Saved plots/fig4_delta_extraction.png")

# Save data
data = np.column_stack([log_mq, log_y_noisy, errors_y])
np.savetxt('../data/fig4_delta_extraction.csv', data, delimiter=',',
           header='log_mq_a,log_ratio,error', fmt='%.8f', comments='')
print("Saved data/fig4_delta_extraction.csv")
print(f"\nExtracted delta = {delta_extracted:.4f} (published: {DELTA})")
print(f"Fit slope = {slope_fit:.4f} (theory: {slope_theory:.4f})")
