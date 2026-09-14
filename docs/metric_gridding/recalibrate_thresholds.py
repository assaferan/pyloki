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

from pyloki.detection import schemes, thresholding  # noqa: E402
from pyloki.detection.thresholding import DynamicThresholdScheme  # noqa: E402

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


def run_viterbi(outdir: str = "docs/metric_gridding/schemes") -> dict[str, dict]:
    """The real thing: Viterbi-optimised scheme per strategy at P_d = 0.1.

    Unblocked by the fix for upstream issue #8 -- before it, `run()` wrote zero-filled
    state records and `backtrack_best` could not find a path (D32).

    Run one strategy per process: building several configs in one interpreter kills it.
    """
    import tempfile

    import numpy as np

    ref = nsegments() // 2
    out = {}
    for name, kw in STRATEGIES.items():
        cfg = make_config(name, **kw)
        bp = cfg.generate_branching_pattern(kind="poly_taylor_moving", ref_seg=ref)
        dyn = DynamicThresholdScheme(
            bp, ref_ducy=REF_DUCY, nbins=cfg.nbins, ntrials=1024, nprobs=30,
            prob_min=0.05, snr_final=TARGET_SNR, nthresholds=100,
            ducy_max=0.3, wtsp=cfg.wtsp, beam_width=2.5, mode="improved",
        )
        dyn.run(thres_neigh=11)
        with tempfile.TemporaryDirectory() as td:
            analyser = schemes.DynamicThresholdSchemeAnalyser.from_file(
                dyn.save(outdir=td),
            )
        best = analyser.backtrack_best(min_probs=[P_D])
        si = best[0] if isinstance(best, (tuple, list)) else best
        out[name] = {
            "prod_b": float(np.prod(bp)),
            "thresholds": si.thresholds,
            "p_d": float(si.get_info("success_h1_cumul")[-1]),
            "log2_complexity": float(np.log2(si.get_info("complexity_cumul")[-1])),
            "log2_cost": float(np.log2(si.get_info("cost")[-1])),
        }
    return out


def run() -> dict[str, dict[str, float]]:
    """Non-optimised fallback, kept for comparison with the pre-fix numbers (D33)."""
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
