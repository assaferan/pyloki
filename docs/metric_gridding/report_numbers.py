"""Single source of truth for every figure quoted in `04_upstream_report.md`.

Why this exists. Twice in one afternoon a sentence in the report contradicted a table
in the same document that was itself correct: the 1.89 nearest-template figure sat three
sections from a 7.5x headline it disagreed with, and a claim that `quadrature` and
`conservative` "describe trees the code would refuse to build" was written over a table
showing Chebyshev+`quadrature` at exactly the limit, where it builds. The measurements
survived review; the prose over them did not.

So no figure in the report is typed by hand. Everything is computed here, written to
`report_numbers.json`, and checked by `tests/test_report_numbers.py`, which asserts that
the committed JSON still reproduces *and* that every number appearing in the report's
tables is one this script produced. Claims that are really booleans -- "does this
configuration fit `branch_max`?" -- are computed as booleans rather than left to prose.

    python docs/metric_gridding/report_numbers.py            # cheap groups only
    python docs/metric_gridding/report_numbers.py --full     # everything, ~12 min

The expensive groups (the 90-cell comparisons) are only recomputed under `--full`; the
test checks the cheap ones live and the expensive ones against the committed JSON.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import phase3_config as P  # noqa: E402
from amplitude_loss import loss_for_leaves  # noqa: E402
from nearest_template import build_sets, min_excursion  # noqa: E402
from nearest_template_cheby import (  # noqa: E402
    build_sets_cheby,
    cheby_basis,
    min_excursion as min_excursion_cheby,
)
from pyloki.utils import psr_utils  # noqa: E402
from pyloki.utils.misc import C_VAL  # noqa: E402

JSON_PATH = HERE / "report_numbers.json"

F0 = 1.0 / P.PERIOD
PO = P.POLY_ORDER
NBINS = P.NBINS
STAGES_CMP = list(range(4, 63, 4))
STAGES_EARLY = list(range(1, 11))     # where 99% of detection losses occur
STAGES_LATE = list(range(32, 63, 4))  # where ~0% do
# The window where covering leaves are ACTUALLY lost, from survival_profile.py on 40
# real pruning runs: first-loss levels 13,13,17,19,24,27,36,41,41,44,47. None before 13.
STAGES_DECISION = [13, 17, 19, 24, 27, 36, 41, 44, 47]
STAGES_AMP = list(range(6, 63, 8))
N_SIGNALS = 6
SEED = 20260917
DUCY = (0.05, 0.10, 0.20)

BASES = {
    "taylor": (build_sets, min_excursion, None),
    "chebyshev": (build_sets_cheby, min_excursion_cheby, cheby_basis),
}


def _truths() -> np.ndarray:
    cfg = P.make_config("aggressive")
    dp = np.array(cfg.get_dparams_actual(cfg.niters_ffa), dtype=float)
    dp[-1] *= C_VAL / F0
    rng = np.random.default_rng(SEED)
    return np.array([0.001, 0.05, 1.0, 0.0]) + rng.uniform(
        -0.5, 0.5, size=(N_SIGNALS, PO)) * dp


# --- cheap groups ---------------------------------------------------------------

def corner_forms() -> dict:
    """Nominal cell-corner phase error, in units of eta/N_b, for both grids."""
    out = {"taylor": {}, "chebyshev": {}, "closed_form": {}}
    for k in range(2, 9):
        step = psr_utils.poly_taylor_step_d_vec(
            k, 16.78, NBINS, 1.0, np.array([F0]), t_ref=0, use_cheby=True)[0]
        asc = step[::-1]
        ph = sum((d / 2.0) * 16.78 ** (i + 1) / math.factorial(i + 1)
                 for i, d in enumerate(asc))
        out["taylor"][k] = round(float(F0 * ph / C_VAL * NBINS), 6)
        cstep = psr_utils.poly_cheb_step_vec(k, NBINS, 1.0, np.array([F0]))[0]
        out["chebyshev"][k] = round(
            float((F0 / C_VAL) * np.sum(cstep / 2.0) * NBINS), 6)
        out["closed_form"][k] = {"taylor": 2.0 ** (k - 1) - 0.5, "chebyshev": k / 2.0}
    return out


def covering_radius(k_max: int = PO, span: int = 2) -> dict:
    """Distance from the cell corner to the NEAREST grid point, not its own centre."""
    x = np.linspace(0.0, 1.0, 4001)
    coef = [2.0 ** j for j in range(k_max)]

    def sup(a):
        return float(np.abs(sum(a[j] * coef[j] * x ** (j + 1)
                                for j in range(k_max))).max())

    own = sup([0.5] * k_max)
    half = min(sup(list(a)) for a in itertools.product([0.5, -0.5], repeat=k_max))
    wide = min(sup([0.5 - n for n in ns])
               for ns in itertools.product(range(-span, span + 1), repeat=k_max))
    return {"own_cell": round(own, 6),
            "best_half_step": round(half, 6),
            "best_wide_lattice": round(wide, 6),
            "lattice_span": span}


def branch_max_counts() -> dict:
    """Max per-axis child count per basis and strategy, and whether it BUILDS.

    The guard is strict (`num_points > branch_max`), so a count equal to `branch_max`
    is fine. Recorded as a boolean so no sentence has to restate it.
    """
    out = {}
    truth = np.array([0.001, 0.05, 1.0, 0.0])
    for basis, (bs, _, _) in BASES.items():
        out[basis] = {}
        for st in ("aggressive", "quadrature", "conservative"):
            cfg = P.make_config(st)
            worst = 0
            for S in (10, 30, 62):
                sets, _, _ = bs(cfg, st, PO, F0, truth, S, n_seed=1)
                for q in sets[1:]:
                    for j in range(PO):
                        worst = max(worst, len(np.unique(np.round(q[:, j], 12))))
            out[basis][st] = {"max_per_axis": int(worst),
                              "branch_max": int(cfg.branch_max),
                              "builds": bool(worst <= cfg.branch_max)}
    return out


# --- expensive groups -----------------------------------------------------------

def nearest_template(basis: str) -> dict:
    bs, me, _ = BASES[basis]
    truths = _truths()
    cfg = P.make_config("aggressive")
    tol = cfg.eta / cfg.nbins
    agg, n_exact = {}, 0
    for ti, tr in enumerate(truths):
        for S in STAGES_CMP:
            sets, tau, sc = bs(cfg, "aggressive", PO, F0, tr, S, n_seed=1)
            v, ex, _ = me(sets, tau, sc, PO, tol, max_nodes=2_000_000)
            agg[(ti, S)] = (v, ex)
            n_exact += ex
    vals = np.array([v for v, _ in agg.values()])
    out = {"cells": len(agg),
           "aggressive": {"median": round(float(np.median(vals)), 6),
                          "min": round(float(vals.min()), 6),
                          "max": round(float(vals.max()), 6),
                          "exact_cells": int(n_exact)}}
    for st in ("quadrature", "conservative"):
        cfg2 = P.make_config(st)
        tol2 = cfg2.eta / cfg2.nbins
        better = worse = unres = 0
        ratios = []
        for ti, tr in enumerate(truths):
            for S in STAGES_CMP:
                a, a_ex = agg[(ti, S)]
                if not a_ex:
                    unres += 1
                    continue
                sets, tau, sc = bs(cfg2, st, PO, F0, tr, S, n_seed=1)
                v, ex, _ = me(sets, tau, sc, PO, tol2, max_nodes=1_000_000,
                              init_best=a)
                if v < a - 1e-12:
                    better += 1
                    ratios.append(a / max(v, 1e-12))
                elif ex:
                    worse += 1
                else:
                    unres += 1
        out[st] = {"closer": better, "not_closer": worse, "unresolved": unres,
                   "median_proven_gain": (round(float(np.median(ratios)), 6)
                                          if ratios else None)}
    return out


def stage_split(basis: str) -> dict:
    """The same comparison, split by whether the stage decides detection.

    A survival model on the calibrated ladder puts 99% of first failures in stages 1-10,
    so a gain averaged over stages 4-60 is mostly measured where nothing is at stake.
    This records both windows so the report cannot quote the aggregate as though it were
    the detection-relevant number.
    """
    bs, me, _ = BASES[basis]
    truths = _truths()
    cfg_a, cfg_q = P.make_config("aggressive"), P.make_config("quadrature")
    tol = cfg_a.eta / cfg_a.nbins
    out = {}
    for label, stages in (("early", STAGES_EARLY), ("late", STAGES_LATE),
                          ("decision", STAGES_DECISION)):
        better = worse = unres = 0
        ratios, aggv = [], []
        for tr in truths:
            for S in stages:
                sets, tau, sc = bs(cfg_a, "aggressive", PO, F0, tr, S, n_seed=1)
                v, ex, _ = me(sets, tau, sc, PO, tol, max_nodes=2_000_000)
                if not ex:
                    unres += 1
                    continue
                aggv.append(v)
                s2, t2, c2 = bs(cfg_q, "quadrature", PO, F0, tr, S, n_seed=1)
                v2, ex2, _ = me(s2, t2, c2, PO, tol, max_nodes=2_000_000, init_best=v)
                if v2 < v - 1e-12:
                    better += 1
                    ratios.append(v / max(v2, 1e-12))
                elif ex2:
                    worse += 1
                else:
                    unres += 1
        out[label] = {"stages": [stages[0], stages[-1]],
                      "aggressive_median": round(float(np.median(aggv)), 6),
                      "closer": better, "not_closer": worse, "unresolved": unres,
                      "median_proven_gain": (round(float(np.median(ratios)), 6)
                                             if ratios else None)}
    return out


def amplitude(basis: str) -> dict:
    """Median loss per strategy, and the paired advantage.

    Both statistics are recorded because they are not interchangeable and the report
    previously quoted one alongside the other two as though they combined: the median of
    the per-cell ratio is not the ratio of the two medians. The per-cell median is the
    one to lead with (it is paired, so cell-to-cell scatter cancels), but a reader who
    divides the quoted medians must be told why they get a different number.
    """
    bs, me, bf = BASES[basis]
    truths = _truths()
    per_cell = {st: {(d, f): [] for d in DUCY for f in ("boxcar", "matched")}
                for st in ("aggressive", "quadrature")}
    for st in ("aggressive", "quadrature"):
        cfg = P.make_config(st)
        tol = cfg.eta / cfg.nbins
        for tr in truths:
            for S in STAGES_AMP:
                sets, tau, sc = bs(cfg, st, PO, F0, tr, S, n_seed=1)
                _, _, _, top = me(sets, tau, sc, PO, tol,
                                  max_nodes=300_000, keep=30)
                for d in DUCY:
                    for f in ("boxcar", "matched"):
                        per_cell[st][(d, f)].append(
                            loss_for_leaves(top, tau, F0, PO, NBINS, d,
                                            filt=f, basis_fn=bf))
    out = {"cells": len(truths) * len(STAGES_AMP)}
    for st in ("aggressive", "quadrature"):
        out[st] = {f"ducy_{d:.2f}": {
            f: round(float(np.median(per_cell[st][(d, f)])), 6)
            for f in ("boxcar", "matched")} for d in DUCY}
    adv = {}
    for d in DUCY:
        adv[f"ducy_{d:.2f}"] = {}
        for f in ("boxcar", "matched"):
            la = np.array(per_cell["aggressive"][(d, f)])
            lq = np.array(per_cell["quadrature"][(d, f)])
            paired = np.median((1 - lq) / (1 - la)) - 1.0
            of_medians = (1 - np.median(lq)) / (1 - np.median(la)) - 1.0
            adv[f"ducy_{d:.2f}"][f] = {
                "per_cell_median": round(float(paired), 6),
                "ratio_of_medians": round(float(of_medians), 6)}
    out["quadrature_advantage"] = adv
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--full", action="store_true",
                    help="recompute the 90-cell comparisons too (~12 min)")
    ap.add_argument("--only", choices=["amplitude", "nearest", "split"], default=None,
                    help="with --full, recompute just one expensive group")
    args = ap.parse_args()

    data = {
        "config": {"tobs_s": round(P.NSAMPS * P.TSAMP, 4), "nsegments": 64,
                   "poly_order": PO, "nbins": NBINS, "eta": P.ETA,
                   "period_s": P.PERIOD, "signals": N_SIGNALS, "seed": SEED,
                   "stages_comparison": STAGES_CMP, "stages_amplitude": STAGES_AMP},
        "corner_forms": corner_forms(),
        "covering_radius": covering_radius(),
        "branch_max": branch_max_counts(),
    }
    if args.full and args.only == "split":
        prev = json.loads(JSON_PATH.read_text()) if JSON_PATH.exists() else {}
        data.update({k: v for k, v in prev.items() if k not in data})
        for basis in BASES:
            data[f"stage_split_{basis}"] = stage_split(basis)
        JSON_PATH.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
        print(json.dumps(data["stage_split_taylor"], indent=2))
        print(json.dumps(data["stage_split_chebyshev"], indent=2))
        return
    if args.full:
        prev = json.loads(JSON_PATH.read_text()) if JSON_PATH.exists() else {}
        for basis in BASES:
            if args.only in (None, "nearest"):
                data[f"nearest_template_{basis}"] = nearest_template(basis)
            else:
                data[f"nearest_template_{basis}"] = prev[f"nearest_template_{basis}"]
            if args.only in (None, "nearest"):
                data[f"stage_split_{basis}"] = stage_split(basis)
            elif f"stage_split_{basis}" in prev:
                data[f"stage_split_{basis}"] = prev[f"stage_split_{basis}"]
            if args.only in (None, "amplitude"):
                data[f"amplitude_{basis}"] = amplitude(basis)
            else:
                data[f"amplitude_{basis}"] = prev[f"amplitude_{basis}"]
    else:
        if JSON_PATH.exists():
            prev = json.loads(JSON_PATH.read_text())
            for k, v in prev.items():
                data.setdefault(k, v)
        print("(cheap groups recomputed; expensive groups carried over from JSON)")

    JSON_PATH.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
    print(json.dumps(data, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
