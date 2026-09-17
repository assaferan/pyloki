"""Does the overflow ratchet actually fire, and can a stricter ladder stop it?

Two questions that saturation cannot answer, now that `threshold_eff` exists.

**Q1 — does it fire?** Median saturation 0.768 (§6.6) says the buffer is not full at the
*end* of a level. The ratchet fires on overflow *during* one and never relaxes, so a
run can end below the buffer having spent most of the level ratcheted. "Saturation < 1"
is therefore consistent with heavy ratcheting and with none. `threshold_eff` separates
them: it records the cut actually applied per level, so `eff > thresh` counts the levels
where the scheme's threshold was not what ran.

**Q2 — can a stricter ladder fix it?** §6.6's option 4. The candidate count is set by
how much the scheme admits, and both shipped ladders target `P_d` = 0.1031. A ladder
backtracked at a smaller `P_d` has higher thresholds at every stage and admits fewer
leaves, so `quadrature` might run its whole search at the nominal scheme. Both arms move
to the stricter operating point together, so the comparison stays a fair equal-`P_d`
one; what changes is which operating point the campaign speaks about.

The success criterion is **not** saturation. It is `eff == thresh` at every level: the
search ran the ladder it was given. That is the property the campaign needs and the one
§6.6 showed `quadrature` lacks at every affordable buffer.

Ladder generation backtracks several `P_d` targets out of ONE `DynamicThresholdScheme`
run, so the strict and baseline ladders are matched -- same scheme object, same
stochastic draw (`thresholding.py` is unseeded, so ladders from separate runs are not
comparable).

    python ratchet_probe.py --ladders --arm quadrature --pds 0.1031,0.03,0.01
    python ratchet_probe.py --dir D --arm quadrature --max-sugg 1048576 --scheme S.npz
    python ratchet_probe.py --report D
"""

from __future__ import annotations

import argparse
import json
import logging
import re
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
from pyloki.detection import schemes  # noqa: E402
from pyloki.detection.thresholding import DynamicThresholdScheme  # noqa: E402
from pyloki.ffa import DynamicProgramming  # noqa: E402
from pyloki.io.timeseries import TimeSeries  # noqa: E402
from pyloki.periodogram import ScatteredPeriodogram  # noqa: E402
from pyloki.prune import prune_dyp_tree  # noqa: E402
from pyloki.utils.misc import C_VAL  # noqa: E402

DUCY = 0.10
M_RECOVER = 1.0
REF_DUCY = 0.1
TARGET_SNR = 10.0
INJECTION = Injection(freq=P.F0, accel=1.0, jerk=0.05, snap=0.001)
LINE = re.compile(r"score thresh:\s*([-\d.]+), eff:\s*([-\d.naN]+)")


# ------------------------------------------------------------------ ladders

def make_ladders(arm: str, pds: list[float], outdir: Path) -> None:
    """One DynamicThresholdScheme run, backtracked at several P_d targets.

    The expensive part is `run()`; `backtrack_best` is a path choice over the states it
    produced. Taking every ladder from one run makes them matched -- and since
    `DynamicThresholdScheme.__init__` is unseeded, ladders from separate runs would not
    be.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    cfg = P.make_config(arm)
    bp = cfg.generate_branching_pattern(kind="poly_chebyshev_moving",
                                        ref_seg=P.nsegments() // 2)
    print(f"{arm}: prod(B) = {np.prod(np.asarray(bp)):.3e}, {len(bp)} stages", flush=True)
    t0 = time.perf_counter()
    dyn = DynamicThresholdScheme(
        np.asarray(bp, dtype=np.float64), ref_ducy=REF_DUCY, nbins=cfg.nbins,
        ntrials=1024, nprobs=30, prob_min=0.01, snr_final=TARGET_SNR,
        nthresholds=100, ducy_max=cfg.ducy_max, wtsp=cfg.wtsp,
        beam_width=2.5, mode="improved")
    dyn.run(thres_neigh=11)
    print(f"  scheme run: {time.perf_counter() - t0:.1f} s", flush=True)
    with tempfile.TemporaryDirectory() as td:
        analyser = schemes.DynamicThresholdSchemeAnalyser.from_file(dyn.save(outdir=td))
    for pd in pds:
        best = analyser.backtrack_best(min_probs=[pd])
        si = best[0] if isinstance(best, (tuple, list)) else best
        thr = np.asarray(si.thresholds, dtype=float)
        succ = np.asarray(si.get_info("success_h1_cumul"), dtype=float)
        path = outdir / f"cheby_{arm}_pd{pd:g}.npz"
        np.savez(path, thresholds=thr, succ_h1=succ, p_d=float(succ[-1]),
                 bp=np.asarray(bp), requested_pd=pd)
        print(f"  P_d req {pd:<7g} -> realised {succ[-1]:.4f}   "
              f"thresholds {thr.min():.2f}..{thr.max():.2f}   {path.name}", flush=True)


# -------------------------------------------------------------------- probe

def _recovered(df, cfg) -> tuple[bool, float]:
    if len(df) == 0:
        return False, float("inf")
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
    return bool(mism.min() <= M_RECOVER), float(mism.min())


def probe(outdir: Path, arm: str, max_sugg: int, scheme: Path, tag: str) -> None:
    logging.disable(logging.INFO)
    cfg = P.make_config(arm)
    thr = np.load(scheme)["thresholds"]
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
            logs = list(Path(td).rglob("*_log.txt"))
            text = max(logs, key=lambda q: q.stat().st_size).read_text() if logs else ""
        pairs = [(float(a), float(b)) for a, b in LINE.findall(text)
                 if b.lower() != "nan"]
        ratcheted = [(t, e) for t, e in pairs if e > t + 1e-9]
        rec, mism = _recovered(df, cfg)
        rows.append({
            "file": f.name, "levels": len(pairs), "levels_ratcheted": len(ratcheted),
            "max_excess": max((e - t for t, e in ratcheted), default=0.0),
            "first_ratchet_level": (pairs.index(ratcheted[0]) + 1) if ratcheted else None,
            "ncand": len(df), "saturation": len(df) / max_sugg,
            "recovered": rec, "mismatch": mism,
            "seconds": time.perf_counter() - t0,
        })
        r = rows[-1]
        print(f"{arm} {f.name} ratcheted {r['levels_ratcheted']:2d}/{r['levels']:2d} "
              f"max_excess {r['max_excess']:5.2f} ncand {r['ncand']:7d} "
              f"sat {r['saturation']:.3f} rec={rec} {r['seconds']:.1f}s", flush=True)
    (outdir / f"probe_{tag}.json").write_text(json.dumps(rows, indent=1))


def report(outdir: Path) -> None:
    print(f"\n{'cell':<34}{'runs':>5}{'ratcheted levels':>18}{'clean runs':>12}"
          f"{'med ncand':>11}{'rec':>6}")
    for p in sorted(outdir.glob("probe_*.json")):
        rows = json.loads(p.read_text())
        lv = np.array([r["levels_ratcheted"] for r in rows])
        tot = np.array([r["levels"] for r in rows])
        clean = int((lv == 0).sum())
        print(f"{p.stem[len('probe_'):]:<34}{len(rows):>5}"
              f"{f'{lv.sum()}/{tot.sum()} ({lv.mean():.1f}/run)':>18}"
              f"{f'{clean}/{len(rows)}':>12}"
              f"{np.median([r['ncand'] for r in rows]):>11.0f}"
              f"{sum(r['recovered'] for r in rows):>4}/{len(rows)}")
    print("\n'clean runs' = runs where eff == thresh at EVERY level, i.e. the search ran"
          "\nthe ladder it was given. That is the property the campaign needs.")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", type=Path)
    ap.add_argument("--arm", choices=("aggressive", "quadrature"))
    ap.add_argument("--max-sugg", type=int, default=2**18)
    ap.add_argument("--scheme", type=Path)
    ap.add_argument("--tag", type=str)
    ap.add_argument("--ladders", action="store_true")
    ap.add_argument("--pds", type=str, default="0.1031,0.03,0.01")
    ap.add_argument("--out", type=Path, default=HERE / "schemes")
    ap.add_argument("--report", type=Path)
    a = ap.parse_args()
    if a.ladders:
        make_ladders(a.arm, [float(x) for x in a.pds.split(",")], a.out)
    elif a.report:
        report(a.report)
    else:
        tag = a.tag or f"{a.arm}_{a.max_sugg}_{a.scheme.stem}"
        probe(a.dir, a.arm, a.max_sugg, a.scheme, tag)


if __name__ == "__main__":
    main()
