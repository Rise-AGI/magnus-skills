"""
Master script: runs all figure reproductions for Bret (2009).
Generates data/ and plots/ directories.
"""
import subprocess
import sys
import os
import time

scripts = [
    "fig1_parallel_growth.py",
    "fig5_hierarchy.py",
]

os.chdir(os.path.dirname(os.path.abspath(__file__)) if "__file__" in dir() else ".")
os.environ["MPLBACKEND"] = "Agg"

total_start = time.time()
results = {}

for script in scripts:
    print(f"\n{'='*60}")
    print(f"Running {script}...")
    print(f"{'='*60}")
    start = time.time()
    ret = subprocess.run([sys.executable, script], capture_output=False)
    elapsed = time.time() - start
    results[script] = {"returncode": ret.returncode, "time": elapsed}
    print(f"  -> {'SUCCESS' if ret.returncode == 0 else 'FAILED'} ({elapsed:.1f}s)")

total = time.time() - total_start
print(f"\n{'='*60}")
print(f"All scripts completed in {total:.1f}s")
for s, r in results.items():
    status = "OK" if r["returncode"] == 0 else "FAIL"
    print(f"  {s}: {status} ({r['time']:.1f}s)")
