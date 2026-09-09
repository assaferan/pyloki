"""Paired-replicate driver for the poly_basis comparison (see HANDOFF.md section 11/12).

Each replicate regenerates the data (unseeded RNG => independent realization) and runs
both bases on that same realization, so realization variance cancels.

Usage: python replicate_basis.py N   (appends to replicates_basis.jsonl)
"""

import json
import pathlib
import subprocess
import sys

import numpy as np

BASE = pathlib.Path("/Users/assaferan/Documents/GitHub/pyloki/scratch/econ_experiment")
PY = "/Users/assaferan/Documents/GitHub/pyloki/.venv/bin/python"
OUT = BASE / "replicates_basis.jsonl"
KEYS = ("best_score", "best_score_ep", "closest_score", "closest_score_ep",
        "n_cands", "elapsed", "pred_leaves")

n_reps = int(sys.argv[1]) if len(sys.argv) > 1 else 8


def run(cmd):
    env = {"PATH": "/usr/bin:/bin", "PYTHONPATH": str(BASE)}
    r = subprocess.run(cmd, capture_output=True, text=True, env=env, check=False)
    if r.returncode != 0:
        tail = (r.stdout + r.stderr).strip().splitlines()[-5:]
        print(f"  FAILED rc={r.returncode}: {' | '.join(tail)}", flush=True)
    return r.returncode


done = 0
for rep in range(n_reps):
    print(f"--- replicate {rep + 1}/{n_reps} ---", flush=True)
    if run([PY, str(BASE / "make_data_snap.py")]) != 0:
        continue
    if any(run([PY, str(BASE / "run_basis.py"), b]) != 0 for b in ("taylor", "chebyshev")):
        continue
    row = {"rep": rep}
    for b in ("taylor", "chebyshev"):
        z = np.load(BASE / f"basis_{b}.npz")
        for k in KEYS:
            row[f"{b}_{k}"] = float(z[k])
    with OUT.open("a") as f:
        f.write(json.dumps(row) + "\n")
    done += 1
    print(f"  closest: taylor={row['taylor_closest_score']:.4f} "
          f"cheby={row['chebyshev_closest_score']:.4f} "
          f"delta={row['chebyshev_closest_score'] - row['taylor_closest_score']:+.4f}",
          flush=True)

print(f"\n{done}/{n_reps} replicates completed -> {OUT}")
