"""Cost-at-matched-tolerance variant of run_basis.py, with eta exposed.

The plain basis comparison is confounded: chebyshev explores ~2x the leaves, so a
higher score could simply be "more trials" rather than a better-conditioned grid.
The efficiency claim in HANDOFF.md section 11 is specifically that chebyshev reaches a
given tolerance MORE CHEAPLY. To test that, run taylor with eta reduced until its
worst-case cell error matches chebyshev at eta=1 (2.0/7.5 = 0.267 at poly_order=4)
and compare cost and score against chebyshev at eta=1.

Original docstring follows.

Recovered-SNR comparison between poly_basis="taylor" and "chebyshev".

Unlike run_snr.py this needs no variant_naive: both arms run the same source and
differ only in the enumeration basis (and the matching branching-pattern kind).

Motivation (HANDOFF.md section 11): the Taylor step is coarsened by 2**k, which pushes
the worst-case sup-norm phase error per grid cell to 7.5x the nominal eta/nbins at
poly_order=4, against 2.0x for Chebyshev. If that matters, the Chebyshev arm should
recover a higher score for the injected signal, at the cost of exploring more leaves.

Usage: python run_basis.py <taylor|chebyshev>
"""

import argparse
import shutil
import sys
import time

import numpy as np

from pyloki.config import ParamLimits, PulsarSearchConfig
from pyloki.detection import thresholding
from pyloki.ffa import DynamicProgramming
from pyloki.io.timeseries import TimeSeries
from pyloki.periodogram import ScatteredPeriodogram
from pyloki.prune import prune_dyp_tree

BASE = "/Users/assaferan/Documents/GitHub/pyloki/scratch/econ_experiment"

ap = argparse.ArgumentParser()
ap.add_argument("basis", choices=["taylor", "chebyshev"])
ap.add_argument("--eta", type=float, default=1.0)
ap.add_argument("--n-runs", type=int, default=1)
ap.add_argument("--max-sugg-pow", type=int, default=16)
ap.add_argument("--snap-frac", type=float, default=0.05)
ap.add_argument("--branch-max", type=int, default=16)
args = ap.parse_args()

d = np.load(f"{BASE}/data_snap.npz")
ts_e, ts_v, dt = d["ts_e"], d["ts_v"], float(d["dt"])
accel, jerk, snap = float(d["accel"]), float(d["jerk"]), float(d["snap"])
freq, nsamps = float(d["freq"]), int(d["nsamps"])
nseg, nbins = int(d["nseg"]), int(d["nbins"])
brute_div, accel_widen = int(d["brute_div"]), float(d["accel_widen"])
tobs = nsamps * dt
tim_data = TimeSeries(ts_e, ts_v, dt)

bseg_ffa = nsamps // nseg
bseg_brute = nsamps // brute_div
half = abs(snap) * args.snap_frac
p = ParamLimits.from_upper(
    (freq - 1.0, freq + 1.0), [snap, jerk, accel], (snap - half, snap + half), tobs
)
limits = np.array(p.limits, dtype=np.float64, copy=True)
if accel_widen != 1.0:
    lo, hi = limits[2]
    mid, halfw = 0.5 * (lo + hi), 0.5 * (hi - lo) * accel_widen
    limits[2] = (mid - halfw, mid + halfw)

search_cfg = PulsarSearchConfig(
    nsamps=nsamps, tsamp=dt, nbins=nbins, eta=args.eta, param_limits=limits,
    bseg_brute=bseg_brute, bseg_ffa=bseg_ffa, prune_poly_order=4,
    ducy_max=0.5, wtsp=1.2, use_fourier=True,
    tiling_strategy="aggressive", branch_max=args.branch_max,
)
dyp = DynamicProgramming(tim_data, search_cfg)
dyp.initialize()
dyp.execute()

kind = f"poly_{args.basis}_moving"
branching_pattern = search_cfg.generate_branching_pattern(
    kind=kind, ref_seg=dyp.nsegments // 2
)
cum = np.cumprod(np.asarray(branching_pattern, dtype=np.float64))
print(f"[{args.basis}] eta={args.eta}")
print(f"[{args.basis}] kind={kind} pattern={np.array2string(np.asarray(branching_pattern), precision=2)}")
print(f"[{args.basis}] predicted total leaves={cum.sum():.4g} peak={cum.max():.4g}",
      flush=True)

thresholds = np.linspace(1.5, 6.0, len(branching_pattern))
thresholding.evaluate_scheme(
    thresholds, branching_pattern, ref_ducy=0.1, nbins=nbins,
    ntrials=1024, snr_final=9.0, ducy_max=0.5, wtsp=1.2,
)

outdir = f"{BASE}/results_basiseta_{args.basis}/"
shutil.rmtree(outdir, ignore_errors=True)
t0 = time.time()
result_file = prune_dyp_tree(
    dyp, thresholds, n_runs=args.n_runs, max_sugg=2**args.max_sugg_pow,
    outdir=outdir, file_prefix=f"basiseta_{args.basis}",
    poly_basis=args.basis, n_workers=1, use_moving_grid=True,
)
elapsed = time.time() - t0
print(f"[{args.basis}] Pruning {elapsed:.1f}s", flush=True)

pgram = ScatteredPeriodogram.load(result_file)
data = pgram.data
print(f"[{args.basis}] TRUE: snap={snap:.4g} jerk={jerk:g} accel={accel:g} "
      f"freq={freq:.10f}")
print(f"[{args.basis}] n_candidates={len(data)}")
if len(data) == 0:
    print(f"[{args.basis}] NO CANDIDATES SURVIVED")
    sys.exit(0)

best = data.loc[data["score"].idxmax()]
scale_a = max(abs(accel), 1.0)
dist = ((data["accel"] - accel) / scale_a) ** 2 + ((data["freq"] - freq) / 1e-4) ** 2
closest = data.loc[dist.idxmin()]
print(f"[{args.basis}] BEST    score={best['score']:.4f} score_ep={best['score_ep']:.4f} "
      f"accel={best.get('accel', float('nan')):.4f} freq={best.get('freq', float('nan')):.8f}")
print(f"[{args.basis}] CLOSEST score={closest['score']:.4f} "
      f"score_ep={closest['score_ep']:.4f} "
      f"accel={closest.get('accel', float('nan')):.4f} "
      f"freq={closest.get('freq', float('nan')):.8f}")
np.savez(
    f"{BASE}/basis_eta_{args.basis}.npz",
    best_score=float(best["score"]), best_score_ep=float(best["score_ep"]),
    closest_score=float(closest["score"]), closest_score_ep=float(closest["score_ep"]),
    n_cands=len(data), elapsed=elapsed, pred_leaves=float(cum.sum()),
)
print(f"[{args.basis}] DONE", flush=True)
