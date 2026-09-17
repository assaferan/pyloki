"""Pilot for the paired `tiling_strategy` campaign: pairing, cost, and p_disc.

NOT a production campaign. This exists to establish three things the power calculation
in `injection_power.py` cannot establish on its own:

1. **that the two arms can be paired at all.** The library never seeds its noise --
   `simulation/pulse.py:338` is a bare `np.random.default_rng()` and there is no
   `seed`/`rng` parameter on `PulseSignalConfig`, `PulseSignalConfig.generate`,
   `DynamicProgramming`, `Pruning` or `prune_dyp_tree`. So the only route to an
   identical noise realisation in both arms is to generate the `TimeSeries` once and
   persist `(ts_e, ts_v, dt)`. `--make` does that; `--arm` consumes it. Everything
   downstream of the time series is deterministic (verified: repeated runs on the same
   saved series give bit-identical score vectors), so this is sufficient.

2. **the measured per-pair wall clock**, which is what turns a required n into a
   feasibility verdict. Numba JIT dominates the first run in a process, so the cost that
   matters is the *marginal* one, reported here per repetition.

3. **a real p_disc**, the discordance rate the power calculation can only model.

Two things this pilot deliberately does NOT do: it does not try to measure the effect
(n is far too small, by design), and it does not report a p-value. The analysis is
pre-committed in `05_injection_design.md` and a pilot is not a test.

Usage (one arm per process -- building several configs in one interpreter is unstable,
per the note in `run_injections.py`):

    python injection_pilot.py --make  --n 24 --snr 12 --dir /tmp/pilot
    python injection_pilot.py --arm aggressive --dir /tmp/pilot
    python injection_pilot.py --arm quadrature --dir /tmp/pilot
    python injection_pilot.py --report --dir /tmp/pilot
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import phase3_config as P  # noqa: E402
from injection_recovery import Injection  # noqa: E402

from pyloki.core import metric  # noqa: E402
from pyloki.ffa import DynamicProgramming  # noqa: E402
from pyloki.io.timeseries import TimeSeries  # noqa: E402
from pyloki.periodogram import ScatteredPeriodogram  # noqa: E402
from pyloki.prune import prune_dyp_tree  # noqa: E402
from pyloki.simulation.pulse import PulseSignalConfig  # noqa: E402
from pyloki.utils.misc import C_VAL  # noqa: E402

SCHEME = str(HERE / "schemes" / "cheby_{}.npz")
ARMS = ("aggressive", "quadrature")
DUCY = 0.10
M_RECOVER = 1.0
INJECTION = Injection(freq=P.F0, accel=1.0, jerk=0.05, snap=0.001)


def make(outdir: Path, n: int, snr: float) -> None:
    """Generate and persist n noise+signal realisations, shared by both arms."""
    outdir.mkdir(parents=True, exist_ok=True)
    cfg = P.make_config("aggressive")
    for i in range(n):
        sig = PulseSignalConfig(
            period=1.0 / INJECTION.freq, dt=cfg.tsamp, nsamps=cfg.nsamps,
            snr=snr, ducy=DUCY,
            mod_kwargs=INJECTION.mod_kwargs(cfg.prune_poly_order))
        tim = sig.generate(shape="gaussian")
        np.savez(outdir / f"tim_{i:04d}.npz",
                 ts_e=tim.ts_e, ts_v=tim.ts_v, dt=tim.dt, snr=snr)
        print(f"made {i + 1}/{n}", flush=True)
    (outdir / "meta.json").write_text(json.dumps({"n": n, "snr": snr, "ducy": DUCY}))


def _recovered(df, cfg) -> tuple[bool, float, float]:
    """Recovery judged in the metric of the full accumulated baseline.

    The same criterion `injection_recovery.py` uses, and for the same reason: a fixed
    tolerance in Hz would favour whichever grid happens to be finer in frequency. The
    periodogram reports physical Taylor parameters whatever the internal basis, so this
    works unchanged for `poly_basis="chebyshev"`.
    """
    if len(df) == 0:
        return False, float("inf"), float("nan")
    po = cfg.prune_poly_order
    tobs = cfg.nsamps * cfg.tsamp
    g = metric.poly_phase_metric(0.0, -tobs / 2, tobs / 2, po, INJECTION.freq,
                                 cfg.nbins, cfg.ducy_max)
    truth = INJECTION.truth_vector(po)
    names = cfg.param_names
    cand = np.column_stack([df[n].to_numpy() for n in names[:-1]])
    vel = C_VAL * (1.0 - df[names[-1]].to_numpy() / INJECTION.freq)
    delta = np.column_stack([cand, vel]) - truth
    mism = np.einsum("ni,ij,nj->n", delta, g, delta)
    b = int(np.argmin(mism))
    return bool(mism[b] <= M_RECOVER), float(mism[b]), float(df["score"].to_numpy()[b])


def run_arm(outdir: Path, arm: str, max_sugg: int) -> None:
    logging.disable(logging.INFO)
    cfg = P.make_config(arm)
    thr = np.load(SCHEME.format(arm))["thresholds"]
    rows = []
    for f in sorted(outdir.glob("tim_*.npz")):
        z = np.load(f)
        tim = TimeSeries(z["ts_e"], z["ts_v"], float(z["dt"]))
        t0 = time.perf_counter()
        dyp = DynamicProgramming(tim, cfg)
        dyp.initialize()
        dyp.execute()
        with tempfile.TemporaryDirectory() as td:
            rf = prune_dyp_tree(dyp, thr, n_runs=1, max_sugg=max_sugg, outdir=td,
                                file_prefix="p", poly_basis="chebyshev", n_workers=1,
                                use_moving_grid=True)
            df = ScatteredPeriodogram.load(str(rf)).data
        secs = time.perf_counter() - t0
        rec, mism, score = _recovered(df, cfg)
        rows.append({"file": f.name, "recovered": rec, "mismatch": mism,
                     "score": score, "ncand": len(df), "seconds": secs,
                     # If the candidate buffer saturates, the realised cut is
                     # max(stage threshold, top-K, median) and is no longer the
                     # scheme value (world_tree.py:527-551). Silent in the log, so
                     # it is recorded here and checked in --report.
                     "saturation": len(df) / max_sugg})
        print(f"{arm} {f.name} rec={rec} mism={mism:.4g} score={score:.3f} "
              f"ncand={len(df)} sat={len(df) / max_sugg:.3f} {secs:.2f}s", flush=True)
    (outdir / f"arm_{arm}.json").write_text(json.dumps(rows, indent=1))


def report(outdir: Path) -> None:
    res = {a: json.loads((outdir / f"arm_{a}.json").read_text()) for a in ARMS}
    a, b = res["aggressive"], res["quadrature"]
    assert [r["file"] for r in a] == [r["file"] for r in b], "pairing broken"
    n = len(a)
    ra = np.array([r["recovered"] for r in a])
    rb = np.array([r["recovered"] for r in b])
    n01 = int((~ra & rb).sum())     # quadrature only
    n10 = int((ra & ~rb).sum())     # aggressive only
    disc = n01 + n10
    # marginal (>= 1 run) is fine for a pilot; the first run in a process pays JIT
    sec_a = np.median([r["seconds"] for r in a[1:]] or [a[0]["seconds"]])
    sec_b = np.median([r["seconds"] for r in b[1:]] or [b[0]["seconds"]])
    sat = max(r["saturation"] for r in (*a, *b))
    print(f"\npairs                 {n}")
    print(f"recovered aggressive  {int(ra.sum())}/{n}")
    print(f"recovered quadrature  {int(rb.sum())}/{n}")
    print(f"discordant            {disc}/{n}   "
          f"(n01={n01} quad-only, n10={n10} agg-only)")
    print(f"p_disc (pilot)        {disc / n:.3f}   "
          f"[Wilson 95%: {_wilson(disc, n)[0]:.3f}, {_wilson(disc, n)[1]:.3f}]")
    print(f"marginal seconds/pair {sec_a + sec_b:.2f}  "
          f"(aggressive {sec_a:.2f} + quadrature {sec_b:.2f})")
    warn = "*** max_sugg BINDS -- realised cut is not the scheme value ***"
    print(f"max buffer saturation {sat:.3f}  {warn if sat > 0.95 else 'ok'}")


def _wilson(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    if n == 0:
        return 0.0, 1.0
    p = k / n
    d = 1 + z**2 / n
    c = (p + z**2 / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / d
    return max(c - h, 0.0), min(c + h, 1.0)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", type=Path, required=True)
    ap.add_argument("--make", action="store_true")
    ap.add_argument("--arm", choices=ARMS)
    ap.add_argument("--report", action="store_true")
    ap.add_argument("--n", type=int, default=24)
    ap.add_argument("--snr", type=float, default=12.0)
    ap.add_argument("--max-sugg", type=int, default=2**14)
    args = ap.parse_args()
    if args.make:
        make(args.dir, args.n, args.snr)
    elif args.arm:
        run_arm(args.dir, args.arm, args.max_sugg)
    elif args.report:
        report(args.dir)


if __name__ == "__main__":
    main()
