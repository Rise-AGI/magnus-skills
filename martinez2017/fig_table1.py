"""
Reproduce Table 1 from Martinez-Fuentes & Herrera-Velazquez (2017).

Computes the distance rho_d at which the linearized (Yukawa) approximation
deviates from the exact numerical solution by 1%, 5%, and 10% for various
values of the plasma parameter N_D.

Outputs: ../data/table1.csv
"""

import sys
import os
import csv

sys.path.insert(0, os.path.dirname(os.path.abspath(".")))
from debye_solver import compute_table1


def main():
    print("Computing Table 1 data...")
    results = compute_table1()

    outdir = os.path.join(os.path.dirname(os.path.abspath(".")), "data")
    os.makedirs(outdir, exist_ok=True)
    outpath = os.path.join(outdir, "table1.csv")

    with open(outpath, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["N_D", "rho_d_1pct", "rho_d_5pct", "rho_d_10pct"])
        for row in results:
            writer.writerow([
                f"{row['N_D']:.1e}",
                f"{row['rho_1pct']:.6e}",
                f"{row['rho_5pct']:.6e}",
                f"{row['rho_10pct']:.6e}",
            ])

    print(f"Table 1 saved to {outpath}")

    # Print table
    print(f"\n{'N_D':>12s} {'rho_d(1%)':>14s} {'rho_d(5%)':>14s} {'rho_d(10%)':>14s}")
    print("-" * 58)
    for row in results:
        print(f"{row['N_D']:>12.1e} {row['rho_1pct']:>14.6e} {row['rho_5pct']:>14.6e} {row['rho_10pct']:>14.6e}")


if __name__ == "__main__":
    main()
