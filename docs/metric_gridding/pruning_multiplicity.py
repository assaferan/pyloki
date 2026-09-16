"""Does tiling redundancy put a *closer* template near the signal?

`sensitivity_loss.py` (D39) sampled offsets uniformly inside a leaf's own cell and
measured the phase excursion from *that* cell's centre. That is the right quantity only
if the signal is carried by exactly one leaf. The box strategies over-claim their
parent width under transport, so neighbouring parents' claimed intervals can overlap and
their children can too -- and `conservative` generates 1e20 more leaves than
`aggressive` for the same final cell size, which has to show up as overlap somewhere.

Where overlap exists, the signal is covered by several leaves and the one that scores
highest is the *nearest*, not the one whose cell happens to contain it. Under pruning
that is what decides survival. So the quantity that matters is

    min over tracked leaves of  max_t |dPhi(t)|   in units of eta / N_b,

which is what this module measures, against the own-cell number D39 reported.

Method: follow the true signal down the tree using the shipped `poly_taylor_branch_batch`
and `poly_taylor_transform_batch`, keeping every leaf within `window` claimed widths of the truth on every axis (so the
eventual nearest leaf is never pruned away by the bookkeeping -- `window=1` already keeps
the full ring of neighbouring cells, and the reported numbers are checked for stability
against a wider window) and transporting the truth with the same `T(delta_t)` as the
leaves.
"""

from __future__ import annotations

import numpy as np

from pyloki.core import metric, taylor
from pyloki.utils import psr_utils, transforms
from pyloki.utils.misc import C_VAL
from pyloki.utils.snail import MiddleOutScheme

from sensitivity_loss import phase_excursion

MAX_TRACKED = 20_000


def follow_truth(
    cfg,
    strategy: str,
    poly_order: int,
    f0: float,
    truth_kin: np.ndarray,
    window: float = 1.0,
    n_seed: int = 1,
    n_time: int = 256,
) -> dict[str, np.ndarray]:
    """Track leaves near `truth_kin` and report nearest/own-cell excursion per stage."""
    nseg = int(np.ceil(cfg.nsamps / cfg.bseg_ffa))
    scheme = MiddleOutScheme(nseg, nseg // 2, cfg.tseg_ffa, stride=1)
    tol = cfg.eta / cfg.nbins

    dparam = np.array(cfg.get_dparams_actual(cfg.niters_ffa), dtype=float)
    dparam[-1] *= C_VAL / f0

    # Seed a BLOCK of base cells around the truth, not just the one containing it. The
    # real search seeds the whole base grid, so a signal that drifts out of one family's
    # descendants is picked up by a neighbour's. Seeding a single cell would charge that
    # drift to the strategy and overstate the loss (it produced a spurious 86x outlier
    # for `aggressive`).
    centre = np.round(truth_kin / dparam) * dparam
    offs = np.array(np.meshgrid(*[np.arange(-n_seed, n_seed + 1)] * poly_order,
                                indexing="ij")).reshape(poly_order, -1).T
    leaves = np.zeros((len(offs), poly_order + 2, 2))
    leaves[:, :-2, 0] = centre + offs * dparam
    leaves[:, :-2, 1] = dparam
    leaves[:, -1, 0] = f0
    truth = np.zeros(poly_order + 1)
    truth[:poly_order] = truth_kin

    nearest, own, mult, ntrack, mism = [], [], [], [], []
    for lvl in range(1, nseg):
        ref_cur, t_half = scheme.get_current_coord(lvl, moving_grid=True)
        ref_next, _ = scheme.get_coord(lvl)
        delta_t = ref_next - ref_cur

        leaves, _ = taylor.poly_taylor_branch_batch(
            leaves, (ref_cur, t_half), cfg.nbins, cfg.eta, poly_order, cfg.branch_max,
        )

        # Excursion of every tracked leaf centre from the truth, over the leaf's window.
        tau = np.linspace(delta_t - t_half, delta_t + t_half, n_time)
        delta = leaves[:, :-2, 0] - truth[:poly_order]
        exc = phase_excursion(delta, tau, f0, poly_order) / tol

        half = leaves[:, :-2, 1] / 2.0
        inside = np.all(np.abs(delta) <= half + 1e-12 * np.abs(half), axis=1)
        nearest.append(float(exc.min()))
        # Metric mismatch (fractional S/N loss, 2nd order) at the nearest leaf, so the
        # phase-error gain can be converted into an amplitude number.
        g = metric.poly_phase_metric(
            0.0, delta_t - t_half, delta_t + t_half, poly_order, f0,
            cfg.nbins, cfg.ducy_max,
        )
        d_near = delta[int(np.argmin(exc))]
        mism.append(float(d_near @ g @ d_near))
        own.append(float(exc[inside].min()) if inside.any() else np.nan)
        mult.append(int(inside.sum()))
        ntrack.append(int(len(leaves)))

        # Keep a neighbourhood of the truth so the nearest leaf is never lost.
        keep = np.all(np.abs(delta) <= window * leaves[:, :-2, 1], axis=1)
        if not keep.any():
            keep = exc <= np.sort(exc)[0]
        leaves = np.ascontiguousarray(leaves[keep])
        if len(leaves) > MAX_TRACKED:
            order = np.argsort(exc[keep])[:MAX_TRACKED]
            leaves = np.ascontiguousarray(leaves[order])

        leaves = taylor.poly_taylor_transform_batch(
            leaves, (ref_next, t_half), (ref_cur, t_half), strategy, np.empty(0),
        )
        truth = transforms.shift_taylor_params(truth, delta_t)

    return {
        "nearest": np.asarray(nearest),
        "own": np.asarray(own),
        "multiplicity": np.asarray(mult),
        "n_tracked": np.asarray(ntrack),
        "mismatch": np.asarray(mism),
    }
