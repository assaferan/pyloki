"""Paired analysis of replicates_tiling.jsonl (basis x tiling_strategy)."""

import json
import pathlib

import numpy as np
from scipy import stats

BASE = pathlib.Path(__file__).resolve().parent
rows = [json.loads(x) for x in
        (BASE / "replicates_tiling.jsonl").read_text().splitlines() if x.strip()]
n = len(rows)
keys = ("taylor_aggressive", "taylor_conservative",
        "chebyshev_aggressive", "chebyshev_conservative")
v = {k: np.array([r[k] for r in rows]) for k in keys}

print(f"paired replicates: n = {n}   (final-stage max EP score, ref_seg=3)\n")
print(f"{'config':>24} {'mean':>8} {'sd':>7}")
print("-" * 41)
for k in keys:
    print(f"{k:>24} {v[k].mean():>8.3f} {v[k].std(ddof=1):>7.3f}")


def paired(a, b, label):
    d = v[b] - v[a]
    t, p = stats.ttest_rel(v[b], v[a])
    tc = stats.t.ppf(0.975, n - 1)
    sem = d.std(ddof=1) / np.sqrt(n)
    print(f"{label:<46} {d.mean():+7.3f}  95%CI [{d.mean()-tc*sem:+.3f},"
          f" {d.mean()+tc*sem:+.3f}]  p={p:.4g}  {int((d>0).sum())}/{n}")


print("\npaired contrasts:")
paired("taylor_aggressive", "taylor_conservative",
       "cost of the default 'aggressive' in TAYLOR")
paired("chebyshev_aggressive", "chebyshev_conservative",
       "cost of the default 'aggressive' in CHEBYSHEV")
paired("taylor_aggressive", "chebyshev_aggressive",
       "cheby - taylor, both 'aggressive' (the sec 12/13 setup)")
paired("taylor_conservative", "chebyshev_conservative",
       "cheby - taylor, both 'conservative' (corrected)")
