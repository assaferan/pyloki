"""Derive a calibrated threshold ladder per basis, once, and cache it.

run_basis.py originally passed a hardcoded np.linspace(1.5, 6.0, nstages), inherited
from run_search.py. That is badly off: at this config the calibrated final-stage
threshold is ~3.6 against the 6.0 being used, which discards exactly the marginal
candidates where a finer grid would show an advantage.

The ladder is derived the way examples/optimal_thresholds.ipynb builds its "constant"
scheme: survival probability per stage = 1 / branching_factor, so the survivor
population stays constant. That is basis-adaptive (each basis gets the ladder suited to
its own branching pattern) and equalizes computational load between the arms, which also
addresses the leaf-count confound.

Derived once and cached because determine_scheme uses an unseeded RNG; re-deriving per
replicate would inject threshold noise into a paired comparison.
"""

import numpy as np

from pyloki.config import ParamLimits, PulsarSearchConfig
from pyloki.detection import thresholding
from pyloki.ffa import DynamicProgramming
from pyloki.io.timeseries import TimeSeries

BASE = "/Users/assaferan/Documents/GitHub/pyloki/scratch/econ_experiment"
NTRIALS = 8192  # high, for a stable ladder

d = np.load(f"{BASE}/data_snap.npz")
dt = float(d["dt"])
accel, jerk, snap = float(d["accel"]), float(d["jerk"]), float(d["snap"])
freq, nsamps = float(d["freq"]), int(d["nsamps"])
nseg, nbins = int(d["nseg"]), int(d["nbins"])
brute_div, accel_widen = int(d["brute_div"]), float(d["accel_widen"])
tobs = nsamps * dt

half = abs(snap) * 0.05
p = ParamLimits.from_upper(
    (freq - 1.0, freq + 1.0), [snap, jerk, accel], (snap - half, snap + half), tobs
)
limits = np.array(p.limits, dtype=np.float64, copy=True)
lo, hi = limits[2]
mid, halfw = 0.5 * (lo + hi), 0.5 * (hi - lo) * accel_widen
limits[2] = (mid - halfw, mid + halfw)

search_cfg = PulsarSearchConfig(
    nsamps=nsamps, tsamp=dt, nbins=nbins, eta=1, param_limits=limits,
    bseg_brute=nsamps // brute_div, bseg_ffa=nsamps // nseg, prune_poly_order=4,
    ducy_max=0.5, wtsp=1.2, use_fourier=True,
    tiling_strategy="aggressive", branch_max=16,
)
dyp = DynamicProgramming(TimeSeries(d["ts_e"], d["ts_v"], dt), search_cfg)
dyp.initialize()
dyp.execute()
ref_seg = dyp.nsegments // 2

print(f"old (hardcoded, both arms): {np.linspace(1.5, 6.0, 3)}")
for basis in ("taylor", "chebyshev"):
    bp = np.asarray(
        search_cfg.generate_branching_pattern(
            kind=f"poly_{basis}_moving", ref_seg=ref_seg
        ),
        dtype=np.float64,
    )
    probs = 1.0 / bp
    scheme = thresholding.determine_scheme(
        probs, bp, ref_ducy=0.1, nbins=nbins, ntrials=NTRIALS,
        snr_final=9.0, ducy_max=0.5, wtsp=1.2,
    )
    th = np.asarray(scheme.thresholds, dtype=np.float64)
    np.save(f"{BASE}/thresholds_{basis}.npy", th)
    print(f"{basis:>10}: pattern={bp}  probs={np.round(probs, 4)}")
    print(f"{'':>10}  thresholds={np.round(th, 4)}  -> thresholds_{basis}.npy")
