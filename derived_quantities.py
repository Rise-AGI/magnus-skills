"""
Derived physical quantities: quark masses, condensates, renormalization.
Computes all main results of Chiu & Hsieh (2003).
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(".")))

import numpy as np

from qxpt_analysis import (
    compute_quark_masses, compute_condensates,
    compute_topological_susceptibility, delta_from_chi_t,
    eta_prime_mass, quark_mass_ratio,
    DELTA, A1, B_FIT, FPI_F0, A_INV_GEV, A_FM, Z_S,
    M_PI_PHYS, M_K_PHYS, F_PI_PHYS, N_F
)

print("=" * 70)
print("FULL RESULTS: Chiu & Hsieh (2003)")
print("Quenched Lattice QCD with Exact Chiral Symmetry")
print("=" * 70)

# 1. Lattice parameters
print("\n--- Lattice Parameters ---")
print(f"beta = 6.0, lattice = 16^3 x 32")
print(f"a^{{-1}} = {A_INV_GEV:.3f} GeV  (published: 1.979(6) GeV)")
print(f"a = {A_FM:.4f} fm  (published: 0.0997(3) fm)")
print(f"f_pi*a = {FPI_F0:.4f}  (published: 0.0667(2))")

# 2. qXPT parameters
print("\n--- qXPT Parameters (Eq. 14 fit) ---")
print(f"delta = {DELTA:.3f}  (published: 0.164(13))")
print(f"A1 = {A1:.3f}  (published: 1.044(22))")
print(f"B = {B_FIT:.3f}  (published: 2.077(130))")

# 3. Topological analysis
print("\n--- Topological Susceptibility ---")
a4_chi_t, _, mean_Q2 = compute_topological_susceptibility()
delta_topo = delta_from_chi_t(a4_chi_t, FPI_F0)
meta = eta_prime_mass(a4_chi_t, FPI_F0, A_INV_GEV)
print(f"<Q^2> = {mean_Q2:.2f}")
print(f"a^4*chi_t = {a4_chi_t:.2e}  (published: 6.03(75)e-5)")
print(f"delta (from chi_t) = {delta_topo:.3f}  (published: 0.16(2))")
print(f"m_eta' = {meta*1000:.0f} MeV  (published: 813(51) MeV)")

# 4. Quark masses
print("\n--- Light Quark Masses ---")
qm = compute_quark_masses()
print(f"m_s/m = {qm['ms_over_m']:.2f}  (published: 22.58(23))")
print(f"m_bare = {qm['m_bare_MeV']:.1f} MeV  (published: 4.7(5) MeV)")
print(f"m_s_bare = {qm['ms_bare_MeV']:.0f} MeV  (published: 107(11) MeV)")
print(f"Z_m = 1/Z_s = {1/Z_S:.4f}")
print(f"m_{{u,d}}^MSbar(2 GeV) = {qm['m_ud_msbar_MeV']:.1f} MeV  (published: 4.1(3) MeV)")
print(f"m_s^MSbar(2 GeV) = {qm['ms_msbar_MeV']:.0f} MeV  (published: 92(9) MeV)")

# 5. Condensates
print("\n--- Chiral and Quark-Gluon Condensates ---")
cond = compute_condensates()
print(f"<qbar q> = {cond['qqbar_GeV3']:.4f} GeV^3  (published: -0.0134(2) GeV^3)")
print(f"<qbar q>_MSbar(2 GeV) = -({cond['qqbar_msbar_MeV_cuberoot']:.0f} MeV)^3  (published: -(250(3) MeV)^3)")
print(f"g<qbar sigma F q> = {cond['qgq_GeV5']:.5f} GeV^5  (published: -0.0170(4) GeV^5)")
print(f"g<qbar sigma F q>_MSbar = -({cond['qgq_msbar_MeV_fifthroot']:.0f} MeV)^5  (published: -(434(4) MeV)^5)")
print(f"M_0^2 = {cond['M02_GeV2']:.2f} GeV^2  (published: 0.98(2) GeV^2)")

# Save summary to file
with open('../data/derived_quantities.csv', 'w') as f:
    f.write("quantity,value,unit,published_value\n")
    f.write(f"a_inverse,{A_INV_GEV:.4f},GeV,1.979(6)\n")
    f.write(f"a,{A_FM:.5f},fm,0.0997(3)\n")
    f.write(f"delta_pion_fit,{DELTA:.4f},,0.164(13)\n")
    f.write(f"delta_topology,{delta_topo:.4f},,0.16(2)\n")
    f.write(f"m_eta_prime,{meta*1000:.1f},MeV,813(51)\n")
    f.write(f"ms_over_m,{qm['ms_over_m']:.3f},,22.58(23)\n")
    f.write(f"m_ud_msbar,{qm['m_ud_msbar_MeV']:.2f},MeV,4.1(3)\n")
    f.write(f"m_s_msbar,{qm['ms_msbar_MeV']:.1f},MeV,92(9)\n")
    f.write(f"qqbar_msbar,{cond['qqbar_msbar_MeV_cuberoot']:.1f},MeV_cuberoot,250(3)\n")
    f.write(f"qgq_msbar,{cond['qgq_msbar_MeV_fifthroot']:.1f},MeV_fifthroot,434(4)\n")
    f.write(f"M02,{cond['M02_GeV2']:.3f},GeV2,0.98(2)\n")
print("\nSaved data/derived_quantities.csv")
