"""
Figure 1: Bresenham algorithm path construction.
Reproduces the illustration of lattice path for C = (5, 3).
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(".")))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from lattice_vectors import bresenham_2d

# Parameters from Figure 1
cmax, cmin = 5, 3

# Generate path
path = bresenham_2d(cmax, cmin)

# Save path data
data_dir = os.path.join(os.path.dirname(os.path.abspath(".")), "data")
os.makedirs(data_dir, exist_ok=True)
with open(os.path.join(data_dir, "fig1_bresenham_path.csv"), "w") as f:
    f.write("step,x,y\n")
    for i, (x, y) in enumerate(path):
        f.write(f"{i},{x},{y}\n")

# Plot
fig, ax = plt.subplots(1, 1, figsize=(6, 5))

# Draw grid
for i in range(cmax + 2):
    ax.axvline(i, color='lightgray', linewidth=0.5)
for j in range(cmin + 2):
    ax.axhline(j, color='lightgray', linewidth=0.5)

# Draw the ideal line
ax.plot([0, cmax], [0, cmin], 'k--', linewidth=0.8, alpha=0.5, label='Ideal line')

# Draw the Bresenham path
xs = [p[0] for p in path]
ys = [p[1] for p in path]
ax.plot(xs, ys, 'b-o', markersize=8, linewidth=2, label='Bresenham path')

# Mark start and end
ax.plot(0, 0, 'gs', markersize=12, zorder=5)
ax.plot(cmax, cmin, 'rs', markersize=12, zorder=5)
ax.annotate('(0,0)', (0, 0), textcoords="offset points", xytext=(-15, -15), fontsize=10)
ax.annotate(f'({cmax},{cmin})', (cmax, cmin), textcoords="offset points", xytext=(5, 5), fontsize=10)

ax.set_xlabel('max-direction', fontsize=12)
ax.set_ylabel('min-direction', fontsize=12)
ax.set_title(f'Bresenham path for C = ({cmax}, {cmin})', fontsize=13)
ax.set_xlim(-0.5, cmax + 0.5)
ax.set_ylim(-0.5, cmin + 0.5)
ax.set_aspect('equal')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)

plots_dir = os.path.join(os.path.dirname(os.path.abspath(".")), "plots")
os.makedirs(plots_dir, exist_ok=True)
plt.tight_layout()
plt.savefig(os.path.join(plots_dir, "fig1_bresenham.png"), dpi=150)
print("Figure 1 saved.")
