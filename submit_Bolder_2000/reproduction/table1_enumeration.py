"""
Table I: R-vector enumeration using the All-R Approach (ARA).
Enumerates all lattice vectors in the range 10 <= R <= 12*sqrt(3) with |Ci| <= 12,
and compares standard approach vs ARA statistics.
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(".")))
import numpy as np

from lattice_vectors import enumerate_r_vectors, classify_vectors

# Parameters from the paper (Section II)
R_min = 10.0
R_max = 12.0 * np.sqrt(3)
C_max_limit = 12

print(f"Enumerating vectors: {R_min:.1f} <= R <= {R_max:.4f}, |C_i| <= {C_max_limit}")

# Enumerate all vectors
vectors_by_R2 = enumerate_r_vectors(R_min, R_max, C_max_limit)

# Classify standard vs ARA
stats = classify_vectors(vectors_by_R2)

print(f"\nTable I Results:")
print(f"{'Method':<12} {'#R values':<12} {'#R vectors':<14} {'#R vectors/R-entry':<20}")
print(f"{'standard':<12} {stats['standard']['n_R']:<12} {stats['standard']['n_vectors']:<14} {stats['standard']['avg_per_R']:<20.1f}")
print(f"{'ARA':<12} {stats['ARA']['n_R']:<12} {stats['ARA']['n_vectors']:<14} {stats['ARA']['avg_per_R']:<20.1f}")

# Paper values: standard: 21, 302, 14.4; ARA: 175, 11486, 65.6
print(f"\nPaper values: standard: 21, 302, 14.4; ARA: 175, 11486, 65.6")

# Save detailed enumeration
data_dir = os.path.join(os.path.dirname(os.path.abspath(".")), "data")
os.makedirs(data_dir, exist_ok=True)

# Save R-value statistics
sorted_R2 = sorted(vectors_by_R2.keys())
with open(os.path.join(data_dir, "table1_r_vectors.csv"), "w") as f:
    f.write("R_squared,R,n_vectors\n")
    for R2 in sorted_R2:
        R = np.sqrt(R2)
        n = len(vectors_by_R2[R2])
        f.write(f"{R2},{R:.8f},{n}\n")

# Save summary
with open(os.path.join(data_dir, "table1_summary.csv"), "w") as f:
    f.write("method,n_R_values,n_R_vectors,avg_vectors_per_R\n")
    f.write(f"standard,{stats['standard']['n_R']},{stats['standard']['n_vectors']},{stats['standard']['avg_per_R']:.1f}\n")
    f.write(f"ARA,{stats['ARA']['n_R']},{stats['ARA']['n_vectors']},{stats['ARA']['avg_per_R']:.1f}\n")

print("Table I data saved.")
