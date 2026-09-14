"""Phase 3 step 3: injection-recovery per strategy, with recalibrated thresholds.

One strategy per process (building several configs in one interpreter kills it), so
this is a CLI: `python run_injections.py <strategy> <n_rep> <snr...>`.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from injection_recovery import Injection, run_one  # noqa: E402
from phase3_config import F0, make_config  # noqa: E402

SCHEMES = Path(__file__).parent / "schemes"
KWARGS = {"metric": {"m_max": 0.00208139, "defer": 16.0}}
# A signal the base grid cannot resolve on its own: pruning has to find all four axes.
INJECTION = Injection(freq=F0, accel=1.0, jerk=0.05, snap=0.001)


def main() -> None:
    strategy = sys.argv[1]
    n_rep = int(sys.argv[2])
    snrs = [float(x) for x in sys.argv[3:]]

    cfg = make_config(strategy, **KWARGS.get(strategy, {}))
    thresholds = np.load(SCHEMES / f"{strategy}.npz")["thresholds"]

    for snr in snrs:
        t0 = time.perf_counter()
        recovered, mism, scores = 0, [], []
        for _ in range(n_rep):
            r = run_one(cfg, INJECTION, snr, max_sugg=2**14, thresholds=thresholds)
            recovered += int(r.recovered)
            mism.append(r.best_mismatch)
            scores.append(r.best_score)
        print(
            f"RESULT {strategy:>13} snr={snr:5.1f} "
            f"recovered={recovered}/{n_rep} "
            f"median_m={np.median(mism):11.4g} "
            f"median_score={np.median(scores):7.3f} "
            f"secs={time.perf_counter() - t0:6.1f}",
            flush=True,
        )


if __name__ == "__main__":
    main()
