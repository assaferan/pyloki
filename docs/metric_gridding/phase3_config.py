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

# 67.1 s over 64 segments. Short of the plan's 18 min, but the same regime: this is
# where prod B(s) showed the metric ahead by ~6 orders (see DECISIONS "Regime
# dependence"), and it keeps an injection-recovery grid affordable.
NSAMPS = 2**20
BSEG_FFA = 2**14


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
