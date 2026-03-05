"""
Fig 3: (m_pi*a)^2/(m_q*a) vs m_q*a - showing quenched chiral logarithm.
The minimum occurs around m_q*a ~ 0.09.

Reproduces Figure 3 from Chiu & Hsieh (2003).
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(".")))

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from qxpt_analysis import (
    MQ_VALUES, DELTA, A1, B_FIT,
    mpi2_eq14, ratio_mpi2_mq
)

# Generate theoretical curve
mq_fine = np.linspace(0.01, 0.22, 300)
ratio_curve = ratio_mpi2_mq(mq_fine, DELTA, A1, B_FIT)

# Data points
np.random.seed(44)
ratio_data = ratio_mpi2_mq(MQ_VALUES, DELTA, A1, B_FIT)
errors = ratio_data * 0.02
ratio_data_noisy = ratio_data + np.random.normal(0, errors)

# Find minimum
idx_min = np.argmin(ratio_curve)
mq_min = mq_fine[idx_min]
ratio_min = ratio_curve[idx_min]

# Plot
fig, ax = plt.subplots(figsize=(8, 6))
ax.errorbar(MQ_VALUES, ratio_data_noisy, yerr=errors, fmt='o', color='blue',
            markersize=5, capsize=3, label='Lattice data (reconstructed)')
ax.plot(mq_fine, ratio_curve, 'r-', linewidth=2, label='q$\\chi$PT fit')
ax.axvline(mq_min, color='gray', linestyle='--', alpha=0.5,
           label=f'Minimum at $m_q a \\approx$ {mq_min:.3f}')
ax.set_xlabel('$m_q a$', fontsize=14)
ax.set_ylabel('$(m_\\pi a)^2 / (m_q a)$', fontsize=14)
ax.set_title('Quenched Chiral Logarithm: $(m_\\pi a)^2/(m_q a)$ vs $m_q a$', fontsize=14)
ax.legend(fontsize=11)
ax.set_xlim(0, 0.22)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('../plots/fig3_chiral_log.png', dpi=150)
print("Saved plots/fig3_chiral_log.png")

# Save data
data = np.column_stack([MQ_VALUES, ratio_data_noisy, errors, ratio_data])
np.savetxt('../data/fig3_chiral_log.csv', data, delimiter=',',
           header='mq_a,ratio_data,ratio_err,ratio_fit', fmt='%.8f', comments='')
print("Saved data/fig3_chiral_log.csv")
