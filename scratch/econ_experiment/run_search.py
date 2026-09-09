import sys
import time

import numpy as np

from pyloki.config import ParamLimits, PulsarSearchConfig
from pyloki.detection import thresholding
from pyloki.ffa import DynamicProgramming
from pyloki.io.timeseries import TimeSeries
from pyloki.periodogram import ScatteredPeriodogram
from pyloki.prune import prune_dyp_tree

variant = sys.argv[1]  # "economized" or "naive"
outdir = f"/Users/assaferan/Documents/GitHub/pyloki/scratch/econ_experiment/results_{variant}/"

d = np.load("/Users/assaferan/Documents/GitHub/pyloki/scratch/econ_experiment/data.npz")
ts_e, ts_v, dt = d["ts_e"], d["ts_v"], float(d["dt"])
period, accel, jerk, freq, nsamps = (
    float(d["period"]), float(d["accel"]), float(d["jerk"]), float(d["freq"]), int(d["nsamps"]),
)
tim_data = TimeSeries(ts_e, ts_v, dt)

eta = 1
nbins = 64
tobs = nsamps * dt
p = ParamLimits.from_upper((freq - 1.0, freq + 1.0), [jerk, accel], (-8.0, 8.0), tobs)
print(f"[{variant}] Param limits: {p.limits}")
bseg_brute = nsamps // 8
bseg_ffa = nsamps // 2

search_cfg = PulsarSearchConfig(
    nsamps=nsamps,
    tsamp=dt,
    nbins=nbins,
    eta=eta,
    param_limits=p.limits,
    bseg_brute=bseg_brute,
    bseg_ffa=bseg_ffa,
    prune_poly_order=3,
    ducy_max=0.5,
    wtsp=1.2,
    use_fourier=True,
    tiling_strategy="aggressive",
    branch_max=16,
)
dyp = DynamicProgramming(tim_data, search_cfg)
dyp.initialize()
print(f"[{variant}] nsegments = {dyp.nsegments}")

t0 = time.time()
dyp.execute()
print(f"[{variant}] FFA execute took {time.time() - t0:.1f}s", flush=True)

branching_pattern = search_cfg.generate_branching_pattern(kind="poly_taylor_moving", ref_seg=dyp.nsegments // 2)

thresholds = np.linspace(1.5, 6.0, len(branching_pattern))
thresh_state = thresholding.evaluate_scheme(
    thresholds,
    branching_pattern,
    ref_ducy=0.1,
    nbins=nbins,
    ntrials=1024,
    snr_final=9.0,
    ducy_max=0.5,
    wtsp=1.2,
)

t0 = time.time()
result_file = prune_dyp_tree(
    dyp,
    thresholds,
    n_runs=2,
    max_sugg=2**16,
    outdir=outdir,
    file_prefix=f"test_{variant}",
    poly_basis="taylor",
    n_workers=1,
    use_moving_grid=True,
)
print(f"[{variant}] Pruning took {time.time() - t0:.1f}s", flush=True)

pgram = ScatteredPeriodogram.load(result_file)
print(f"[{variant}] True params: jerk={jerk:.3f}, accel={accel:.3f}, freq={freq:.10f}")
print(f"[{variant}] N candidates: {len(pgram.data)}")
print(pgram.get_summary_cands(10, score_type="score", run_id=None))

data = pgram.data
if len(data) > 0:
    best = data.loc[data["score"].idxmax()]
    print(f"[{variant}] BEST: score={best['score']:.3f}, score_ep={best['score_ep']:.3f}, "
          f"accel={best.get('accel', float('nan')):.4f}, freq={best.get('freq', float('nan')):.8f}, "
          f"jerk={best.get('jerk', float('nan')):.4f}")
    dist = (data["accel"] - accel) ** 2 + (data["freq"] - freq) ** 2 * 1e6
    closest = data.loc[dist.idxmin()]
    print(f"[{variant}] CLOSEST-TO-TRUE: score={closest['score']:.3f}, score_ep={closest['score_ep']:.3f}, "
          f"accel={closest.get('accel', float('nan')):.4f}, freq={closest.get('freq', float('nan')):.8f}, "
          f"jerk={closest.get('jerk', float('nan')):.4f}")
else:
    print(f"[{variant}] NO CANDIDATES SURVIVED")
print(f"[{variant}] DONE", flush=True)
