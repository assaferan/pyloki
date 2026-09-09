"""Paired analysis of replicates_basis.jsonl (poly_basis taylor vs chebyshev)."""

import json
import pathlib

import numpy as np
from scipy import stats

BASE = pathlib.Path(__file__).resolve().parent
import sys
SRC = sys.argv[1] if len(sys.argv) > 1 else "replicates_basis.jsonl"
rows = [json.loads(x) for x in
        (BASE / SRC).read_text().splitlines() if x.strip()]
n = len(rows)
print(f"paired replicates: n = {n}\n")
if n < 2:
    raise SystemExit("need >= 2 replicates")

print(f"{'metric':<20} {'taylor':>11} {'chebyshev':>11} {'mean delta':>11} "
      f"{'sd(delta)':>10} {'p (2-sided)':>12} {'cheby wins':>11}")
print("-" * 92)
for key in ("closest_score", "closest_score_ep", "best_score", "best_score_ep",
            "n_cands", "elapsed"):
    a = np.array([r[f"taylor_{key}"] for r in rows])
    b = np.array([r[f"chebyshev_{key}"] for r in rows])
    d = b - a
    t, p = stats.ttest_rel(b, a)
    print(f"{key:<20} {a.mean():>11.4f} {b.mean():>11.4f} {d.mean():>+11.4f} "
          f"{d.std(ddof=1):>10.4f} {p:>12.4g} {f'{int((d > 0).sum())}/{n}':>11}")

d = np.array([r["chebyshev_closest_score"] - r["taylor_closest_score"] for r in rows])
print("\nper-replicate delta in closest-to-true score (chebyshev - taylor):")
print("  " + "  ".join(f"{x:+.3f}" for x in d))
sem = d.std(ddof=1) / np.sqrt(n)
tc = stats.t.ppf(0.975, n - 1)
print(f"\n  mean {d.mean():+.4f}   sd {d.std(ddof=1):.4f}   sem {sem:.4f}")
print(f"  95% CI [{d.mean() - tc * sem:+.4f}, {d.mean() + tc * sem:+.4f}]")
t, p = stats.ttest_rel([r["chebyshev_closest_score"] for r in rows],
                       [r["taylor_closest_score"] for r in rows])
print(f"  paired t-test: t={t:.3f}, p={p:.4g} -> "
      f"{'SIGNIFICANT' if p < 0.05 else 'not significant'} at alpha=0.05")

ta = np.array([r["taylor_elapsed"] for r in rows])
cb = np.array([r["chebyshev_elapsed"] for r in rows])
print(f"\n  cost: chebyshev is {cb.mean() / ta.mean():.2f}x the taylor runtime "
      f"({ta.mean():.1f}s -> {cb.mean():.1f}s)")

# Which metric is actually usable? "closest surviving candidate" moves around because
# which candidates clear thresholding depends on the noise, so report the spread of each
# metric alongside the effect, and the replicates needed to resolve it.
print("\nmetric stability (sd of the paired delta, and n needed for 80% power")
print("to resolve an effect of 0.5 in score at alpha=0.05):")
for key in ("closest_score", "closest_score_ep", "best_score", "best_score_ep"):
    a = np.array([r[f"taylor_{key}"] for r in rows])
    b = np.array([r[f"chebyshev_{key}"] for r in rows])
    sd = (b - a).std(ddof=1)
    n_need = ((1.96 + 0.84) * sd / 0.5) ** 2
    print(f"  {key:<20} sd={sd:>7.4f}   n~{n_need:>6.0f}")
