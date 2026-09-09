import pathlib
"""Localize the EP-score deficit: is the Chebyshev drop tied to the LAST STAGE, or to
the specific segment added there?

Section 13 found Chebyshev leads on per-stage max score at stages 1-2 and then DROPS at
stage 3 (11.08 -> 10.26) where Taylor rises (10.92 -> 11.39). With nsegments=4 and
ref_seg=2 the addition order is 2,1,3,0, so stage 3 adds segment 0. Varying ref_seg
changes which segment lands last: if the drop follows the stage index it is a
stage-count effect; if it follows a particular segment it is positional.

Usage: python diag_stages.py            (all four ref_segs, both bases)
"""

import re
import shutil
import subprocess
import sys
import tempfile

import numpy as np

from pyloki.config import ParamLimits, PulsarSearchConfig
from pyloki.detection import thresholding
from pyloki.ffa import DynamicProgramming
from pyloki.io.timeseries import TimeSeries
from pyloki.prune import prune_dyp_tree

BASE = str(pathlib.Path(__file__).resolve().parent)

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

LINE = re.compile(r"Prune level:\s*(\d+), seg_idx:\s*(\d+).*?max:\s*([0-9.]+)")

print(f"nsegments={dyp.nsegments}  tseg={dyp.tseg:.4f}s  poly_order=4\n")
rows = {}
for ref_seg in range(dyp.nsegments):
    for basis in ("taylor", "chebyshev"):
        bp = np.asarray(
            search_cfg.generate_branching_pattern(
                kind=f"poly_{basis}_moving", ref_seg=ref_seg
            ),
            dtype=np.float64,
        )
        scheme = thresholding.determine_scheme(
            1.0 / bp, bp, ref_ducy=0.1, nbins=nbins, ntrials=4096,
            snr_final=9.0, ducy_max=0.5, wtsp=1.2,
        )
        thresholds = np.asarray(scheme.thresholds, dtype=np.float64)
        tmp = tempfile.mkdtemp()
        try:
            prune_dyp_tree(
                dyp, thresholds, ref_segs=[ref_seg], max_sugg=2**16,
                batch_size=1024, outdir=tmp, file_prefix="diag",
                poly_basis=basis, n_workers=1, use_moving_grid=True,
            )
            log = next(p for p in __import__("pathlib").Path(tmp).glob("*log.txt"))
            found = LINE.findall(log.read_text())
            rows[(ref_seg, basis)] = [
                (int(lv), int(si), float(mx)) for lv, si, mx in found
            ]
        finally:
            shutil.rmtree(tmp, ignore_errors=True)

print("per-stage max score (stage: seg_idx -> max)")
print(f"{'ref_seg':>8} {'basis':>10} | " + "  ".join(f"{'stage ' + str(i):>16}" for i in range(dyp.nsegments)))
print("-" * 100)
for ref_seg in range(dyp.nsegments):
    for basis in ("taylor", "chebyshev"):
        r = rows.get((ref_seg, basis), [])
        cells = "  ".join(f"{f'seg{si}: {mx:6.2f}':>16}" for _, si, mx in r)
        print(f"{ref_seg:>8} {basis:>10} | {cells}")
    t = rows.get((ref_seg, "taylor"), [])
    c = rows.get((ref_seg, "chebyshev"), [])
    if t and c:
        print(f"{'':>8} {'delta':>10} | " + "  ".join(
            f"{f'{cc[2] - tt[2]:+6.2f}':>16}" for tt, cc in zip(t, c, strict=False)))
    print()
