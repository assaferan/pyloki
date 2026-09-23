"""Measure the recovery margin of each example test across noise realisations.

Reports error/tolerance per predicate per seed. A ratio >= 1 is a failure; ratios
just under 1 are the realisations worth pinning in a sweep, because they are where
the tolerance is actually being tested.
"""
import logging, sys, tempfile, warnings
from pathlib import Path
logging.disable(logging.CRITICAL); warnings.filterwarnings("ignore")
import numpy as np

SEEDS = [int(x) for x in sys.argv[2:]] or list(range(1, 21))
which = sys.argv[1]


def ratios_accel(seed):
    from tests import test_example_ep_accel as M
    with tempfile.TemporaryDirectory() as td:
        r = M._run_search(seed, Path(td))
    b = r["data"].loc[r["data"]["score"].idxmax()]
    tf = r["true_freq"]
    return {
        "accel": abs(float(b["accel"]) - M.ACCEL) / M.ACCEL_TOL,
        "freq": abs(float(b["freq"]) - tf) / M.FREQ_TOL,
        "accel/daccel": abs(float(b["accel"]) - M.ACCEL) / float(b["daccel"]),
    }


def ratios_circular(seed):
    from tests import test_example_ep_circular as M
    with tempfile.TemporaryDirectory() as td:
        r = M._run_search(seed, Path(td))
    b = r["data"].loc[r["data"]["score"].idxmax()]
    out = {}
    for p, (tol, _) in M.TOLERANCES.items():
        out[p] = abs(float(b[p]) - r["truth"][p]) / tol
        out[p + "/unc"] = abs(float(b[p]) - r["truth"][p]) / max(float(b[f"d{p}"]), 1e-30)
    return out


def ratios_ffa(seed):
    from tests import test_example_ffa as M
    out = {}
    for case in M.CASES:
        r = M._run_case(case, seed)
        worst = 0
        for uf in M.BACKENDS:
            for name, off in M._index_offsets(r["pgrams"][uf], r["truth"]).items():
                worst = max(worst, abs(off))
        out[case.name + " maxoff"] = worst / M.INDEX_TOL
    return out


fn = {"accel": ratios_accel, "circular": ratios_circular, "ffa": ratios_ffa}[which]
rows = {}
for seed in SEEDS:
    try:
        rows[seed] = fn(seed)
    except AssertionError as e:
        rows[seed] = {"EXCEPTION": str(e)[:60]}
    print(f"{which} seed={seed:3d} " +
          "  ".join(f"{k}={v:.3f}" if isinstance(v, float) else f"{k}={v}"
                    for k, v in rows[seed].items()), flush=True)

keys = [k for k in next(iter(rows.values())) if isinstance(rows[SEEDS[0]].get(k), float)]
print(f"\n=== {which}: error / tolerance, {len(SEEDS)} seeds ===")
for k in keys:
    vals = {s: rows[s][k] for s in SEEDS if k in rows[s]}
    order = sorted(vals, key=lambda s: -vals[s])
    fails = [s for s in order if vals[s] >= 1.0]
    print(f"{k:22s} max={vals[order[0]]:.3f} (seed {order[0]})  "
          f"top5={[(s, round(vals[s],3)) for s in order[:5]]}  FAIL={fails}")
