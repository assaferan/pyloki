"""Phase 3 step 2: recalibrate the threshold scheme per strategy at P_d = 0.1.

B(s) differs between strategies, so a scheme tuned for one is meaningless for another --
demonstrated in D31, where a shared placeholder ramp made the metric look like an
outright recovery failure when it was only thresholded away. Nothing about sensitivity
may be quoted until this has been run for every strategy being compared.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from phase3_config import make_config, nsegments  # noqa: E402

from pyloki.detection import thresholding  # noqa: E402

REF_DUCY = 0.1
TARGET_SNR = 10.0
P_D = 0.1

# Matches conservative's worst-case leaf mismatch at this config, so the strategies
# are compared at equal sensitivity rather than equal cost.
METRIC_M_MAX = 0.00208139
METRIC_DEFER = 16.0

STRATEGIES = {
    "aggressive": {},
    "conservative": {},
    "metric": {"m_max": METRIC_M_MAX, "defer": METRIC_DEFER},
}


def run() -> dict[str, dict[str, float]]:
    """Per-strategy threshold scheme, one strategy per process.

    Uses `determine_scheme` with `probs = 1 / B(s)` -- "expect one survivor per branch"
    -- and NOT the Viterbi optimiser, which is broken on this base: it empties every
    state at stage 2, including on the example notebook's own branching pattern, in both
    modes. See D32. These schemes are per-strategy and comparable, but they are not
    P_d-optimal, so anything derived from them is provisional.

    Run one strategy per process: constructing several configs in one interpreter kills
    it silently.
    """
    import numpy as np

    ref = nsegments() // 2
    out = {}
    for name, kw in STRATEGIES.items():
        cfg = make_config(name, **kw)
        bp = cfg.generate_branching_pattern(kind="poly_taylor_moving", ref_seg=ref)
        st = thresholding.determine_scheme(
            np.clip(1.0 / bp, 1e-6, 1.0), bp, REF_DUCY, cfg.nbins,
            ntrials=1024, snr_final=TARGET_SNR, ducy_max=0.3, wtsp=cfg.wtsp,
        )
        out[name] = {
            "prod_b": float(np.prod(bp)),
            "p_d": float(st.get_info("success_h1_cumul")[-1]),
            "log2_complexity": float(np.log2(st.get_info("complexity_cumul")[-1])),
            "log2_cost": float(np.log2(st.get_info("cost")[-1])),
            "thresholds": st.thresholds,
        }
    return out


if __name__ == "__main__":
    import numpy as np

    for name, r in run().items():
        print(f"{name:>13}: prodB={r['prod_b']:.3g} P_d={r['p_d']:.4g} "
              f"log2(complexity)={r['log2_complexity']:.2f} "
              f"log2(cost)={r['log2_cost']:.2f} "
              f"thr=[{np.nanmin(r['thresholds']):.2f},"
              f"{np.nanmax(r['thresholds']):.2f}]")
