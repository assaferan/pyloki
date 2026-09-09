import pathlib
"""End-to-end SNR comparison at poly_order=4: naive truncation vs economization.

Both variants read the SAME data_snap.npz. The variant is selected by which pyloki
tree is on sys.path -- variant_naive/ is a copy of src/pyloki with the three
economize_taylor_params calls stripped, which is exactly the pre-fix behaviour
(resolve only ever reads the low 3 coefficients, so not economizing IS naive
slicing). The runner asserts which one it actually loaded, so a PYTHONPATH mistake
cannot silently produce a bogus "no difference" result.

Usage: python run_snr.py <econ|naive> [--n-runs N] [--max-sugg-pow P]
"""

import argparse
import inspect
import shutil
import sys
import time

import numpy as np

BASE = str(pathlib.Path(__file__).resolve().parent)

ap = argparse.ArgumentParser()
ap.add_argument("variant", choices=["econ", "naive"])
ap.add_argument("--n-runs", type=int, default=1)
ap.add_argument("--max-sugg-pow", type=int, default=16)
ap.add_argument("--snap-frac", type=float, default=0.05)
ap.add_argument("--branch-max", type=int, default=16)
ap.add_argument("--batch-size", type=int, default=1024)
args = ap.parse_args()

from pyloki.config import ParamLimits, PulsarSearchConfig  # noqa: E402
from pyloki.core import taylor  # noqa: E402
from pyloki.detection import thresholding  # noqa: E402
from pyloki.ffa import DynamicProgramming  # noqa: E402
from pyloki.io.timeseries import TimeSeries  # noqa: E402
from pyloki.periodogram import ScatteredPeriodogram  # noqa: E402
from pyloki.prune import prune_dyp_tree  # noqa: E402

# --- Guard: verify the loaded tree matches the requested variant -------------
src = inspect.getsource(taylor.poly_taylor_resolve_batch.py_func)
has_econ = "economize_taylor_params" in src
want_econ = args.variant == "econ"
import pyloki  # noqa: E402

print(f"[{args.variant}] pyloki from : {pyloki.__file__}")
print(f"[{args.variant}] economize in resolve: {has_econ}")
if has_econ != want_econ:
    print(
        f"FATAL: variant '{args.variant}' wants economize={want_econ} but the loaded "
        f"tree has {has_econ}. Set PYTHONPATH to variant_naive/ for 'naive', or leave "
        f"it unset for 'econ'."
    )
    sys.exit(2)

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
    nsamps=nsamps, tsamp=dt, nbins=nbins, eta=1, param_limits=limits,
    bseg_brute=bseg_brute, bseg_ffa=bseg_ffa, prune_poly_order=4,
    ducy_max=0.5, wtsp=1.2, use_fourier=True,
    tiling_strategy="aggressive", branch_max=args.branch_max,
)
dyp = DynamicProgramming(tim_data, search_cfg)
dyp.initialize()
t0 = time.time()
dyp.execute()
print(f"[{args.variant}] FFA execute {time.time() - t0:.1f}s  "
      f"nsegments={dyp.nsegments}  tseg={dyp.tseg:.4f}s  "
      f"grid={np.asarray(dyp.param_grid_count).tolist()}", flush=True)

branching_pattern = search_cfg.generate_branching_pattern(
    kind="poly_taylor_moving", ref_seg=dyp.nsegments // 2
)
thresholds = np.linspace(1.5, 6.0, len(branching_pattern))
thresholding.evaluate_scheme(
    thresholds, branching_pattern, ref_ducy=0.1, nbins=nbins,
    ntrials=1024, snr_final=9.0, ducy_max=0.5, wtsp=1.2,
)

outdir = f"{BASE}/results_{args.variant}/"
shutil.rmtree(outdir, ignore_errors=True)

t0 = time.time()
result_file = prune_dyp_tree(
    dyp, thresholds, n_runs=args.n_runs, max_sugg=2**args.max_sugg_pow, batch_size=args.batch_size,
    outdir=outdir, file_prefix=f"snr_{args.variant}",
    poly_basis="taylor", n_workers=1, use_moving_grid=True,
)
print(f"[{args.variant}] Pruning {time.time() - t0:.1f}s", flush=True)

pgram = ScatteredPeriodogram.load(result_file)
data = pgram.data
print(f"[{args.variant}] TRUE: snap={snap:.4g} jerk={jerk:g} accel={accel:g} "
      f"freq={freq:.10f}")
print(f"[{args.variant}] n_candidates={len(data)}")
if len(data) == 0:
    print(f"[{args.variant}] NO CANDIDATES SURVIVED")
    sys.exit(0)

best = data.loc[data["score"].idxmax()]
print(f"[{args.variant}] BEST score={best['score']:.4f} score_ep={best['score_ep']:.4f} "
      f"accel={best.get('accel', float('nan')):.4f} "
      f"freq={best.get('freq', float('nan')):.8f}")

# Candidate closest to truth, in units of the search's own grid steps.
dacc = np.asarray(dyp.param_grid_count).ravel()
scale_a = max(abs(accel), 1.0)
dist = ((data["accel"] - accel) / scale_a) ** 2 + ((data["freq"] - freq) / 1e-4) ** 2
closest = data.loc[dist.idxmin()]
print(f"[{args.variant}] CLOSEST score={closest['score']:.4f} "
      f"score_ep={closest['score_ep']:.4f} "
      f"accel={closest.get('accel', float('nan')):.4f} "
      f"freq={closest.get('freq', float('nan')):.8f}")
np.savez(
    f"{BASE}/snr_{args.variant}.npz",
    best_score=float(best["score"]), best_score_ep=float(best["score_ep"]),
    closest_score=float(closest["score"]), closest_score_ep=float(closest["score_ep"]),
    n_cands=len(data), grid=dacc,
)
print(f"[{args.variant}] DONE", flush=True)
