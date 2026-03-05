"""
Figure 4: Integration contour in the complex v-plane.

Reproduces Fig. 1 of the paper, showing the integration contours
used for calculating the density perturbation:
- For gamma > 0: singularity at v = omega/k is in upper half-plane,
  real-axis integral + semicircle C_R captures the residue
- For gamma < 0: singularity is in lower half-plane, and the
  S*Delta correction term accounts for the pole contribution

Output: plots/fig4_contour.png
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

fig, axes = plt.subplots(1, 3, figsize=(15, 5))

for ax_idx, (gamma_sign, title) in enumerate([
    ("positive", r"$\gamma > 0$: pole above real axis"),
    ("zero", r"$\gamma = 0$: pole on real axis"),
    ("negative", r"$\gamma < 0$: pole below real axis")
]):
    ax = axes[ax_idx]

    # Draw axes
    ax.axhline(y=0, color='k', linewidth=0.8)
    ax.axvline(x=0, color='k', linewidth=0.8)

    # Integration contour along real axis
    v_real = np.linspace(-3, 3, 100)
    ax.plot(v_real, np.zeros_like(v_real), 'b-', linewidth=2.5, zorder=5)

    # Add arrows on the real axis contour
    ax.annotate("", xy=(1.5, 0), xytext=(1.0, 0),
                arrowprops=dict(arrowstyle="->", color='b', lw=2))
    ax.annotate("", xy=(-1.0, 0), xytext=(-1.5, 0),
                arrowprops=dict(arrowstyle="->", color='b', lw=2))

    # Pole position
    Omega_k = 1.5  # v = Omega/k on real axis
    if gamma_sign == "positive":
        pole_im = 0.5
    elif gamma_sign == "zero":
        pole_im = 0.0
    else:
        pole_im = -0.5

    ax.plot(Omega_k, pole_im, 'rx', markersize=12, markeredgewidth=3, zorder=10)
    ax.annotate(r"$v = \omega/k$", xy=(Omega_k, pole_im),
                xytext=(Omega_k + 0.3, pole_im + 0.3),
                fontsize=11, color='r')

    # Semicircle C_R in upper half-plane (from +inf to -inf)
    theta = np.linspace(0, np.pi, 50)
    R = 2.8
    x_semi = R * np.cos(theta)
    y_semi = R * np.sin(theta)
    ax.plot(x_semi, y_semi, 'g--', linewidth=1.5, alpha=0.6)
    ax.annotate("", xy=(-R * np.cos(np.pi/4), R * np.sin(np.pi/4)),
                xytext=(R * np.cos(np.pi/4), R * np.sin(np.pi/4)),
                arrowprops=dict(arrowstyle="->", color='g', lw=1.5, ls='--'))
    ax.text(0, R * 0.7, r"$C_R$", fontsize=12, color='g', ha='center')

    # For gamma=0, show indentation around pole
    if gamma_sign == "zero":
        theta_ind = np.linspace(0, np.pi, 30)
        r_ind = 0.3
        x_ind = Omega_k + r_ind * np.cos(theta_ind)
        y_ind = r_ind * np.sin(theta_ind)
        ax.plot(x_ind, y_ind, 'b-', linewidth=2.5, zorder=5)

    ax.set_xlabel(r"Re$(v)$", fontsize=12)
    ax.set_ylabel(r"Im$(v)$", fontsize=12)
    ax.set_title(title, fontsize=12)
    ax.set_xlim(-3.5, 3.5)
    ax.set_ylim(-1.5, 3.2)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.2)

plt.tight_layout()
plt.savefig("../plots/fig4_contour.png", dpi=150, bbox_inches="tight")
print("Saved plots/fig4_contour.png")

# Also save a minimal data file for completeness
np.savetxt("../data/fig4_contour.csv",
           np.array([[0, 0]]),
           header="This figure is a schematic diagram with no numerical data",
           fmt="%.1f")
print("Saved data/fig4_contour.csv")
