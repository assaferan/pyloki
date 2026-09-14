"""The Phase 3 config, in one place so every experiment uses the same one.

Deviation from metric_PLAN.md, agreed 2026-09-14 (D29): the plan specifies a
circular-orbit search, but `"metric"` is refused on the circular basis (D16) because
`core/metric.py` builds `g` on the Taylor kinematic basis. Phase 3 therefore runs a
**Taylor analogue** at `poly_order = 4`, in the regime where the branching comparison
showed the mechanism actually biting. The paper's Figures 8/11/12 are circular-orbit
panels, so these reproduce them in spirit, not literally.
"""

from __future__ import annotations

import numpy as np

from pyloki.config import ParamLimits, PulsarSearchConfig

PERIOD = 0.007          # 7 ms spin period, as the plan specifies
F0 = 1.0 / PERIOD
TSAMP = 64e-6
NBINS = 64
ETA = 1.0
DUCY_MAX = 0.2
POLY_ORDER = 4

# 268.4 s over 64 segments. Chosen by measurement, not preference (D30):
#   - the plan's 18 min / 128 segments needs a 3.8 TB FFA fold array -- not runnable
#     on this machine; 537 s / 128 segments needs 5.0 GB and is a spot-check at best;
#   - below ~268 s the FFA base grid is degenerate. At 67 s it is [1, 1, 1, 12]: only
#     frequency is gridded, and every strategy recovers the same candidate, so nothing
#     is discriminated. At 268 s it is [1, 1, 5, 467] and pruning refines all four axes.
# FFA here costs 2.9 s and 0.1 GB, so an injection grid is affordable.
NSAMPS = 2**22
BSEG_FFA = 2**16


def make_config(strategy: str, m_max: float = 0.2, defer: float = 1.0,
                nsamps: int = NSAMPS, bseg_ffa: int = BSEG_FFA,
                poly_order: int = POLY_ORDER) -> PulsarSearchConfig:
    tobs = nsamps * TSAMP
    limits = ParamLimits.from_upper(
        (F0 * (1 - 1e-3), F0 * (1 + 1e-3)),
        [0.3, 6.0, 500.0][-(poly_order - 1):],
        (-8.0, 8.0),
        tobs,
    )
    return PulsarSearchConfig(
        nsamps=nsamps, tsamp=TSAMP, nbins=NBINS, eta=ETA,
        param_limits=limits.limits, bseg_brute=bseg_ffa // 8, bseg_ffa=bseg_ffa,
        prune_poly_order=poly_order, ducy_max=DUCY_MAX, wtsp=1.5, use_fourier=True,
        tiling_strategy=strategy, branch_max=16,
        m_max=m_max, metric_defer_factor=defer, metric_branch_max=5_000_000,
    )


def nsegments(nsamps: int = NSAMPS, bseg_ffa: int = BSEG_FFA) -> int:
    return int(np.ceil(nsamps / bseg_ffa))
