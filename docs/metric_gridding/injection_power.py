"""Power calculation for the injection campaign, committed before any production run.

The question the campaign exists to answer is whether tiling redundancy keeps signals
alive through pruning: a signal near a cell edge may sit in a leaf that is thresholded
away, and a strategy claiming more territory keeps more such leaves alive. That is NOT
the 0.5% amplitude difference measured in D55 -- a thresholded-away leaf loses the signal
entirely -- so the effect is bimodal, and powering against a mean amplitude shift would
understate it. The right quantity is the fraction of signals whose score passes close
enough to the ladder for a 0.5% amplitude difference to change the outcome.

Model, stated so it can be argued with:

- the signal-carrying leaf accumulates segments; at stage `s` its score is
  `a * (1 - loss_s) * sqrt(s+1) + (1/sqrt(s+1)) * sum_{i<=s} n_i` with `n_i ~ N(0,1)`
  i.i.d. per segment, so scores are correlated across stages the way accumulating the
  same data makes them;
- `a = snr / sqrt(nseg)` is the per-segment amplitude;
- `loss_s` is the measured fractional S/N loss of the nearest template at stage `s`,
  from `nearest_template.py` and `amplitude_loss.py`, per strategy;
- the signal is recovered iff its score clears the calibrated threshold at every stage;
- **paired**: the two arms share the noise draw `n_i` and the injected signal, so only
  `loss_s` differs. That is what makes McNemar on discordant pairs the right test.

What it ignores, and why each would have to be argued separately: competing candidates
and the beam (so this is survival of the true leaf, not the full search), a single
`ducy`, Gaussian noise, and the assumption that both arms are calibrated to equal on-grid
`P_d` so that a difference can only come from off-grid loss.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import phase3_config as P  # noqa: E402
from amplitude_loss import loss_for_leaves  # noqa: E402
from nearest_template import build_sets, min_excursion  # noqa: E402
from nearest_template_cheby import (  # noqa: E402
    build_sets_cheby,
    cheby_basis,
    min_excursion as min_excursion_cheby,
)
from pyloki.utils.misc import C_VAL  # noqa: E402

BASES = {
    "taylor": (build_sets, min_excursion, None),
    "chebyshev": (build_sets_cheby, min_excursion_cheby, cheby_basis),
}


def stage_losses(basis, strategy, truths, stages, ducy=0.10, keep=20,
                 max_nodes=200_000):
    """Median-over-signals fractional S/N loss at each stage."""
    bs, me, bf = BASES[basis]
    cfg = P.make_config(strategy)
    tol = cfg.eta / cfg.nbins
    out = np.zeros(len(stages))
    for k, S in enumerate(stages):
        vals = []
        for tr in truths:
            sets, tau, sc = bs(cfg, strategy, P.POLY_ORDER, 1.0 / P.PERIOD, tr,
                               S, n_seed=1)
            _, _, _, top = me(sets, tau, sc, P.POLY_ORDER, tol,
                              max_nodes=max_nodes, keep=keep)
            vals.append(loss_for_leaves(top, tau, 1.0 / P.PERIOD, P.POLY_ORDER,
                                        P.NBINS, ducy, basis_fn=bf))
        out[k] = float(np.median(vals))
    return out


def survival(snr, losses, thresholds, noise):
    """Paired survival: does the leaf clear every threshold? `noise` is (n_mc, nseg).

    `losses` has one entry per stage 0..nseg-1; `thresholds` has one per *prune level*,
    i.e. stages 1..nseg-1, so the comparison skips stage 0 (nothing has branched yet).
    """
    nseg = noise.shape[1]
    a = snr / np.sqrt(nseg)
    s = np.arange(1, nseg + 1)
    walk = np.cumsum(noise, axis=1) / np.sqrt(s)[None, :]
    score = a * (1.0 - losses)[None, :] * np.sqrt(s)[None, :] + walk
    return (score[:, 1:] >= thresholds[None, :]).all(axis=1)


def mcnemar_n(p01, p10, alpha=0.05, power=0.80):
    """Pairs needed to detect the asymmetry between discordant rates."""
    from scipy.stats import norm
    pdisc = p01 + p10
    if pdisc <= 0 or np.isclose(p01, p10):
        return np.inf
    za, zb = norm.ppf(1 - alpha / 2), norm.ppf(power)
    psi = abs(p10 - p01) / pdisc
    n_disc = (za + zb * np.sqrt(1 - psi**2)) ** 2 / psi**2
    return float(np.ceil(n_disc / pdisc))
