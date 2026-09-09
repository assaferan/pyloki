"""Replicate the tiling finding across fresh noise realizations.

diag_tiling.py (one realization, ref_seg=3) indicated the Chebyshev EP-score deficit is
an artifact of tiling_strategy="aggressive" -- the codebase default (config.py:415) --
and that with proper error propagation Chebyshev overtakes Taylor. Single realizations
have misled twice in this thread already, so this repeats it on independent data.

Reports the final-stage max EP score per (basis, tiling), paired within each realization.
"""

import json
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

BASE = Path("/Users/assaferan/Documents/GitHub/pyloki/scratch/econ_experiment")
PY = "/Users/assaferan/Documents/GitHub/pyloki/.venv/bin/python"
OUT = BASE / "replicates_tiling.jsonl"
REF_SEG, BRANCH_MAX = 3, 128
LINE = re.compile(r"Prune level:\s*(\d+), seg_idx:\s*(\d+).*?max:\s*([0-9.]+)")

n_reps = int(sys.argv[1]) if len(sys.argv) > 1 else 4


def one(basis, tiling):
    """Run one config on the CURRENT data_snap.npz, return final-stage max score."""
    import numpy as np

    from pyloki.config import ParamLimits, PulsarSearchConfig
    from pyloki.detection import thresholding
    from pyloki.ffa import DynamicProgramming
    from pyloki.io.timeseries import TimeSeries
    from pyloki.prune import prune_dyp_tree

    d = np.load(BASE / "data_snap.npz")
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
    cfg = PulsarSearchConfig(
        nsamps=nsamps, tsamp=dt, nbins=nbins, eta=1, param_limits=limits,
        bseg_brute=nsamps // brute_div, bseg_ffa=nsamps // nseg, prune_poly_order=4,
        ducy_max=0.5, wtsp=1.2, use_fourier=True,
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
            dyp, np.asarray(scheme.thresholds, dtype=np.float64), ref_segs=[REF_SEG],
            max_sugg=2**16, batch_size=1024, outdir=tmp, file_prefix="rep",
            poly_basis=basis, n_workers=1, use_moving_grid=True,
        )
        log = next(Path(tmp).glob("*log.txt"))
        found = LINE.findall(log.read_text())
        return float(found[-1][2])
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    if len(sys.argv) > 2:  # worker mode: one config, print the number
        print(one(sys.argv[2], sys.argv[3]))
        raise SystemExit(0)

    combos = [(b, t) for b in ("taylor", "chebyshev")
              for t in ("aggressive", "conservative")]
    env = {"PATH": "/usr/bin:/bin", "PYTHONPATH": str(BASE)}
    for rep in range(n_reps):
        print(f"--- replicate {rep + 1}/{n_reps} ---", flush=True)
        r = subprocess.run([PY, str(BASE / "make_data_snap.py")],
                           capture_output=True, text=True, env=env, check=False)
        if r.returncode != 0:
            print("  data gen FAILED", flush=True)
            continue
        row = {"rep": rep}
        ok = True
        for basis, tiling in combos:
            pr = subprocess.run(
                [PY, str(BASE / "replicate_tiling.py"), "1", basis, tiling],
                capture_output=True, text=True, env=env, check=False,
            )
            if pr.returncode != 0:
                tail = (pr.stdout + pr.stderr).strip().splitlines()[-3:]
                print(f"  {basis}/{tiling} FAILED: {' | '.join(tail)}", flush=True)
                ok = False
                break
            row[f"{basis}_{tiling}"] = float(pr.stdout.strip().splitlines()[-1])
        if not ok:
            continue
        with OUT.open("a") as f:
            f.write(json.dumps(row) + "\n")
        print("  " + "  ".join(f"{k}={v:.2f}" for k, v in row.items()
                               if k != "rep"), flush=True)
    print(f"\ndone -> {OUT}")
