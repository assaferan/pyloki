"""Max per-axis `num_points` per (basis, tiling_strategy), measured in a REAL prune.

This replaces `report_numbers.branch_max_counts()`, which was withdrawn in full (D120)
because it did not measure the guarded quantity: it counted unique offsets from
`nearest_template.build_sets`, which gates on the shift before computing `n` and leaves an
axis width unnarrowed when that axis does not branch, while the shipped path calls
`branch_param_padded` for every axis at every level unconditionally and gates afterwards.
It also sampled three stages of 62 at a single `f0`.

The method here does not model the quantity at all -- it asks the code.

    The guard is strict: `branch_param_padded` raises iff `num_points > branch_max`.
    So if a full 63-level prune COMPLETES at `branch_max = B` and RAISES at `B - 1`,
    then max `num_points` over every leaf, every level and every axis is EXACTLY B.

Bisecting `B` therefore reads the true maximum straight out of the shipped guard, with no
reimplementation of `ceil(dparam_cur/dparam_new)` anywhere in this file. That matters:
reimplementing it is precisely how the withdrawn table went wrong.

One wrinkle. `PulsarSearchConfig` validates `branch_max > 10`, so a naive bisection floors
at 11 and reports 11 for every configuration -- which looks like a measurement and is an
artefact of the validator. `branch_max` only sizes a padded buffer, so the validator is
bypassed here with `attrs.validators.disabled()` to reach the real values, which are all
below the floor.

What this establishes and what it does not. Completion is a positive observation: these
configurations DO run at the shipped default, which refutes the withdrawn claim that three
of four are unreachable. It is not a proof that no leaf in any configuration ever exceeds
`branch_max` -- same asymmetry as the reductio that falsified the proxy. Treat the numbers
as measured maxima over the leaves these runs actually produced.

    python docs/metric_gridding/branch_max_probe.py            # the six-cell table
    python docs/metric_gridding/branch_max_probe.py --buffers  # + buffer dependence

Needs PYTHONPATH set to this worktree's src; it asserts that and refuses otherwise.
"""

from __future__ import annotations

import argparse
import sys
import tempfile
import time
from pathlib import Path

import attrs
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import pyloki  # noqa: E402

if "worktrees/metric-gridding/src" not in pyloki.__file__:  # pragma: no cover
    msg = (f"pyloki resolves to {pyloki.__file__}; set PYTHONPATH to this worktree's src "
           "or you are measuring main's library, not this branch's")
    raise RuntimeError(msg)

import phase3_config as P  # noqa: E402
from injection_recovery import Injection  # noqa: E402
from pyloki.ffa import DynamicProgramming  # noqa: E402
from pyloki.prune import Pruning  # noqa: E402
from pyloki.simulation.pulse import PulseSignalConfig  # noqa: E402

F0 = 1.0 / P.PERIOD
INJECTION = Injection(freq=F0, accel=1.0, jerk=0.05, snap=0.001)
SNR = 14.0
DUCY = 0.10
BASES = ("taylor", "chebyshev")
STRATEGIES = ("aggressive", "quadrature", "conservative")


def _config(strategy: str, branch_max: int):
    """Phase 3 config with `branch_max` overridden below the validator floor."""
    with attrs.validators.disabled():
        return attrs.evolve(P.make_config(strategy), branch_max=branch_max)


def _dyp(cfg):
    sig = PulseSignalConfig(
        period=1.0 / INJECTION.freq, dt=cfg.tsamp, nsamps=cfg.nsamps,
        snr=SNR, ducy=DUCY,
        mod_kwargs=INJECTION.mod_kwargs(cfg.prune_poly_order),
    )
    dyp = DynamicProgramming(sig.generate(shape="gaussian"), cfg)
    dyp.initialize()
    dyp.execute()
    return dyp


def completes(basis: str, strategy: str, branch_max: int, max_sugg: int) -> bool:
    """Does a full 63-level prune finish without the guard firing?"""
    cfg = _config(strategy, branch_max)
    dyp = _dyp(cfg)
    thresholds = np.linspace(1.5, 8.0, dyp.nsegments - 1)
    with tempfile.TemporaryDirectory() as td:
        log = Path(td) / "log.txt"
        log.touch()
        prn = Pruning(dyp, thresholds, max_sugg=max_sugg, batch_size=1024,
                      poly_basis=basis, use_moving_grid=True)
        prn.initialize(dyp.nsegments // 2, log)
        for _ in range(dyp.nsegments - 1):
            try:
                prn.execute_iter(log)
            except ValueError as exc:
                if "branch_max" not in str(exc):
                    raise
                return False
    return True


def max_num_points(basis: str, strategy: str, max_sugg: int, hi: int = 32) -> int:
    """Smallest `branch_max` that completes == max `num_points`, by bisection."""
    lo = 2
    while lo < hi:
        mid = (lo + hi) // 2
        if completes(basis, strategy, mid, max_sugg):
            hi = mid
        else:
            lo = mid + 1
    return hi


def table(max_sugg: int = 2**14) -> dict:
    out = {}
    for basis in BASES:
        out[basis] = {}
        for strategy in STRATEGIES:
            t0 = time.time()
            n = max_num_points(basis, strategy, max_sugg)
            out[basis][strategy] = n
            print(f"  {basis}/{strategy}: {n}  ({time.time() - t0:.0f}s)", flush=True)
    return out


def buffer_dependence(cells, buffers=(2**14, 2**16, 2**18)) -> dict:
    """Does the maximum move with the candidate buffer? A bigger buffer keeps more
    leaves, so it has more chances at a large ratio. If these rows are not flat, the
    table above is conditional on `max_sugg` and must say so."""
    out = {}
    for basis, strategy in cells:
        row = [max_num_points(basis, strategy, ms) for ms in buffers]
        out[f"{basis}/{strategy}"] = row
        print(f"  {basis}/{strategy}: " + "  ".join(f"2^{ms.bit_length()-1}={v}"
                                                    for ms, v in zip(buffers, row)),
              flush=True)
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--buffers", action="store_true",
                    help="also check whether the maxima move with max_sugg")
    ap.add_argument("--max-sugg", type=int, default=2**14)
    args = ap.parse_args()

    print(f"max per-axis num_points, real prune, max_sugg=2^{args.max_sugg.bit_length()-1}")
    res = table(args.max_sugg)
    print(f"\n{'':<12}" + "".join(f"{s:>14}" for s in STRATEGIES))
    for basis in BASES:
        print(f"{basis:<12}" + "".join(f"{res[basis][s]:>14}" for s in STRATEGIES))
    worst = max(v for row in res.values() for v in row.values())
    print(f"\nworst = {worst}, shipped branch_max = 16 -> all six build")

    if args.buffers:
        print("\nbuffer dependence:")
        buffer_dependence([("taylor", "conservative"), ("taylor", "quadrature"),
                           ("chebyshev", "conservative")])
