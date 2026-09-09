"""Paired analysis of replicates.jsonl."""

import json
import pathlib

import numpy as np
from scipy import stats

BASE = pathlib.Path("/Users/assaferan/Documents/GitHub/pyloki/scratch/econ_experiment")
rows = [json.loads(x) for x in (BASE / "replicates.jsonl").read_text().splitlines() if x.strip()]
n = len(rows)
print(f"paired replicates: n = {n}\n")
if n < 2:
    raise SystemExit("need >= 2 replicates")

print(f"{'metric':<20} {'naive mean':>11} {'econ mean':>11} {'mean delta':>11} "
      f"{'sd(delta)':>10} {'p (2-sided)':>12} {'econ wins':>10}")
print("-" * 90)
for key in ("closest_score", "closest_score_ep", "best_score", "best_score_ep", "n_cands"):
    a = np.array([r[f"naive_{key}"] for r in rows])
    b = np.array([r[f"econ_{key}"] for r in rows])
    d = b - a
    if np.allclose(d, 0):
        print(f"{key:<20} {a.mean():>11.4f} {b.mean():>11.4f} {0.0:>11.4f} "
              f"{0.0:>10.4f} {'--':>12} {'identical':>10}")
        continue
    t, p = stats.ttest_rel(b, a)
    wins = int((d > 0).sum())
    print(f"{key:<20} {a.mean():>11.4f} {b.mean():>11.4f} {d.mean():>+11.4f} "
          f"{d.std(ddof=1):>10.4f} {p:>12.4g} {f'{wins}/{n}':>10}")

print("\nPer-replicate delta in closest-to-true score (econ - naive):")
d = np.array([r["econ_closest_score"] - r["naive_closest_score"] for r in rows])
print("  " + "  ".join(f"{x:+.3f}" for x in d))
print(f"\n  mean {d.mean():+.4f}   sd {d.std(ddof=1):.4f}   "
      f"sd of mean {d.std(ddof=1) / np.sqrt(n):.4f}")
if n >= 3:
    t, p = stats.ttest_rel(
        [r["econ_closest_score"] for r in rows],
        [r["naive_closest_score"] for r in rows],
    )
    verdict = "significant" if p < 0.05 else "NOT significant"
    print(f"  paired t-test: t={t:.3f}, p={p:.4g} -> {verdict} at alpha=0.05")
