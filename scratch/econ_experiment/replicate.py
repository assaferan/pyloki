"""Paired-replicate driver for the end-to-end SNR comparison.

Each replicate regenerates the data (PulseSignalConfig uses an unseeded
default_rng, so every call is an independent noise realization) and runs BOTH
variants on that same realization. Pairing removes realization-to-realization
variance, which is far larger than the effect under test.

Usage: python replicate.py N   (appends to replicates.jsonl)
"""

import json
import pathlib
import subprocess
import sys

import numpy as np

BASE = pathlib.Path("/Users/assaferan/Documents/GitHub/pyloki/scratch/econ_experiment")
PY = "/Users/assaferan/Documents/GitHub/pyloki/.venv/bin/python"
OUT = BASE / "replicates.jsonl"
KEYS = ("best_score", "best_score_ep", "closest_score", "closest_score_ep", "n_cands")

n_reps = int(sys.argv[1]) if len(sys.argv) > 1 else 5


def run(cmd, env_path):
    env = {"PATH": "/usr/bin:/bin", "PYTHONPATH": env_path}
    r = subprocess.run(cmd, capture_output=True, text=True, env=env, check=False)
    if r.returncode != 0:
        tail = (r.stdout + r.stderr).strip().splitlines()[-6:]
        print(f"  FAILED rc={r.returncode}: {' | '.join(tail)}")
    return r.returncode


done = 0
for rep in range(n_reps):
    print(f"--- replicate {rep + 1}/{n_reps} ---", flush=True)
    if run([PY, str(BASE / "make_data_snap.py")], str(BASE)) != 0:
        continue
    if run([PY, str(BASE / "run_snr.py"), "econ"], str(BASE)) != 0:
        continue
    if run(
        [PY, str(BASE / "run_snr.py"), "naive"],
        f"{BASE / 'variant_naive'}:{BASE}",
    ) != 0:
        continue

    row = {"rep": rep}
    for variant in ("econ", "naive"):
        z = np.load(BASE / f"snr_{variant}.npz")
        for k in KEYS:
            row[f"{variant}_{k}"] = float(z[k])
    with OUT.open("a") as f:
        f.write(json.dumps(row) + "\n")
    done += 1
    print(
        f"  closest: naive={row['naive_closest_score']:.4f} "
        f"econ={row['econ_closest_score']:.4f} "
        f"delta={row['econ_closest_score'] - row['naive_closest_score']:+.4f}",
        flush=True,
    )

print(f"\n{done}/{n_reps} replicates completed -> {OUT}")
