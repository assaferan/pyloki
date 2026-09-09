"""Test whether the Chebyshev EP-score deficit is caused by tiling_strategy.

diag_stages.py showed Taylor's final-stage score is stable across ref_seg
(11.75/11.73/11.39/11.30) while Chebyshev swings (10.11/11.72/10.26/9.32). That points
at error propagation rather than the basis algebra.

tiling_strategy="aggressive" keeps only the diagonal of the transform matrix
(shift_taylor_errors / shift_cheby_full: `errors * abs(diag(t_mat))`). For the Taylor
re-centering matrix the diagonal is all ones, so aggressive tracking is harmless. The
Chebyshev interval-change matrix has non-unit diagonal and heavy off-diagonal mixing, so
diagonal-only tracking should badly underestimate the leaf uncertainty -- the grid then
under-refines and the true signal falls outside the leaf's claimed cell, costing
coherence on later additions.

Run at ref_seg=3, where the deficit is largest (-1.98).
"""

import re
import shutil
import tempfile
from pathlib import Path

import numpy as np

from pyloki.config import ParamLimits, PulsarSearchConfig
from pyloki.detection import thresholding
from pyloki.ffa import DynamicProgramming
from pyloki.io.timeseries import TimeSeries
from pyloki.prune import prune_dyp_tree

BASE = "/Users/assaferan/Documents/GitHub/pyloki/scratch/econ_experiment"
REF_SEG = 3
BRANCH_MAX = 128
LINE = re.compile(r"Prune level:\s*(\d+), seg_idx:\s*(\d+).*?max:\s*([0-9.]+)")

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

print(f"ref_seg={REF_SEG}  poly_order=4\n")
print(f"{'basis':>10} {'tiling':>13} | " + "  ".join(f"{'stage ' + str(i):>14}" for i in range(4)))
print("-" * 82)

for basis in ("taylor", "chebyshev"):
    for tiling in ("aggressive", "quadrature", "conservative"):
        cfg = PulsarSearchConfig(
            nsamps=nsamps, tsamp=dt, nbins=nbins, eta=1, param_limits=limits,
            bseg_brute=nsamps // brute_div, bseg_ffa=nsamps // nseg,
            prune_poly_order=4, ducy_max=0.5, wtsp=1.2, use_fourier=True,
            tiling_strategy=tiling, branch_max=BRANCH_MAX,
        )
        dyp = DynamicProgramming(TimeSeries(d["ts_e"], d["ts_v"], dt), cfg)
        dyp.initialize()
        dyp.execute()
        bp = np.asarray(
            cfg.generate_branching_pattern(kind=f"poly_{basis}_moving", ref_seg=REF_SEG),
            dtype=np.float64,
        )
        scheme = thresholding.determine_scheme(
            1.0 / bp, bp, ref_ducy=0.1, nbins=nbins, ntrials=4096,
            snr_final=9.0, ducy_max=0.5, wtsp=1.2,
        )
        tmp = tempfile.mkdtemp()
        try:
            prune_dyp_tree(
                dyp, np.asarray(scheme.thresholds, dtype=np.float64),
                ref_segs=[REF_SEG], max_sugg=2**16, batch_size=1024,
                outdir=tmp, file_prefix="diag", poly_basis=basis,
                n_workers=1, use_moving_grid=True,
            )
            log = next(Path(tmp).glob("*log.txt"))
            found = LINE.findall(log.read_text())
            cells = "  ".join(
                f"{f'seg{si}: {float(mx):6.2f}':>14}" for _, si, mx in found
            )
            print(f"{basis:>10} {tiling:>13} | {cells}")
        finally:
            shutil.rmtree(tmp, ignore_errors=True)
