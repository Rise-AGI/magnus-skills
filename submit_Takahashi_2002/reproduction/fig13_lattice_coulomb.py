"""
Figure 13: Lattice Coulomb potential V^LC(n) vs continuum 1/r.

Computes V^LC on a lattice and compares with the continuum Coulomb potential.
"""
import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(os.path.join(os.getcwd(), __file__))) if '__file__' not in dir() else os.path.dirname(os.path.abspath(__file__)))
from three_quark_potential import lattice_coulomb_potential


def main():
    # Compute V^LC for various lattice vectors
    vectors = []
    vlc_vals = []
    r_vals = []

    # On-axis vectors
    for n in range(0, 9):
        for vec in [(n, 0, 0)]:
            r = np.sqrt(vec[0]**2 + vec[1]**2 + vec[2]**2)
            if r > 0:
                vlc = lattice_coulomb_potential(vec, N_grid=64)
                vectors.append(vec)
                vlc_vals.append(vlc)
                r_vals.append(r)

    # Off-axis vectors
    off_axis = [(1, 1, 0), (1, 1, 1), (2, 1, 0), (2, 1, 1), (2, 2, 0),
                (2, 2, 1), (2, 2, 2), (3, 1, 0), (3, 1, 1), (3, 2, 0),
                (3, 2, 1), (3, 3, 0), (3, 3, 1), (3, 3, 3),
                (4, 0, 0), (4, 1, 0), (4, 2, 0), (4, 4, 0),
                (5, 0, 0), (5, 1, 0), (6, 0, 0)]
    for vec in off_axis:
        r = np.sqrt(vec[0]**2 + vec[1]**2 + vec[2]**2)
        vlc = lattice_coulomb_potential(vec, N_grid=64)
        vectors.append(vec)
        vlc_vals.append(vlc)
        r_vals.append(r)

    # Also compute V^LC(0,0,0)
    vlc_origin = lattice_coulomb_potential((0, 0, 0), N_grid=64)

    r_vals = np.array(r_vals)
    vlc_vals = np.array(vlc_vals)

    # Continuum 1/r
    r_cont = np.linspace(0.1, 8.5, 200)
    v_cont = 1.0 / r_cont

    # Save data
    os.makedirs('../data', exist_ok=True)
    with open('../data/fig13_lattice_coulomb.csv', 'w') as f:
        f.write('# Lattice Coulomb potential V^LC(n) and continuum 1/r\n')
        f.write('# V^LC(0,0,0) = {:.8f}\n'.format(vlc_origin))
        f.write('n1,n2,n3,r,V_LC,V_continuum\n')
        for i, vec in enumerate(vectors):
            v_cont_val = 1.0 / r_vals[i]
            f.write('{},{},{},{:.8f},{:.8f},{:.8f}\n'.format(
                vec[0], vec[1], vec[2], r_vals[i], vlc_vals[i], v_cont_val))

    # Plot
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))

    # On-axis points
    on_mask = np.array([v[1] == 0 and v[2] == 0 for v in vectors])
    off_mask = ~on_mask

    ax.plot(r_cont, v_cont, 'k-', linewidth=1.5, label=r'Continuum $1/r$')
    ax.plot(r_vals[on_mask], vlc_vals[on_mask], 'ro', markersize=8,
            label=r'$V^{LC}$ (on-axis)')
    ax.plot(r_vals[off_mask], vlc_vals[off_mask], 'b^', markersize=6,
            label=r'$V^{LC}$ (off-axis)')

    # Mark the origin value
    ax.plot(0, vlc_origin, 'gs', markersize=10, label=r'$V^{LC}(\vec{0})$')

    ax.set_xlabel(r'$|\vec{n}|$ (lattice units)', fontsize=14)
    ax.set_ylabel(r'$V^{LC}(\vec{n})$', fontsize=14)
    ax.set_title('Fig. 13: Lattice Coulomb Potential', fontsize=14)
    ax.set_xlim(-0.5, 8.5)
    ax.set_ylim(0, 4.0)
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)

    os.makedirs('../plots', exist_ok=True)
    plt.tight_layout()
    plt.savefig('../plots/fig13_lattice_coulomb.png', dpi=150)
    plt.close()
    print('Fig 13 saved.')


if __name__ == '__main__':
    main()
