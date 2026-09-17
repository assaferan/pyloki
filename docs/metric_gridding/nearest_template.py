"""Exact nearest-template search over the pruning tree, by branch and bound.

Replaces `pruning_multiplicity.py`, which was unsound (D47): it materialised leaves and
kept the best `MAX_TRACKED` by *current* excursion, discarding leaves that would have
become nearest after later branching. A minimum obtained that way is not a minimum, and
it was not monotone in the search width.

The fix is structural rather than a bigger cap. `branch_param_padded` computes a child's
offset from its parent as a function of `(dparam_cur, dparam_new)` only, and every leaf
at a stage shares those, so **the child offsets are identical for every parent**
(verified in `test_offsets_are_parent_independent`). Writing `O_s` for that per-stage
offset set (a Cartesian product of per-axis arithmetic progressions) and `T_s` for the
re-centring matrix, the full leaf set at stage `S` is therefore *exactly*

    L_S  =  M_0 . Seed  (+)  M_1 . O_1  (+)  ...  (+)  O_S

a Minkowski sum of small finite sets, with no approximation anywhere. The leaf count is
the product of the set sizes -- the 1e32 of the `conservative` strategy -- but the
*representation* is tiny, one set of at most a few hundred points per stage.

Minimising the sup-norm phase error over a Minkowski sum is still combinatorial (the sup
over `t` couples the axes), so it is done by branch and bound with an admissible bound:
having fixed the choices `i < m`, whatever the remaining sets contribute at time `t` is
bounded by `sum_{i >= m} max_{a in A_i} |<a, b(t)>|`, so

    |<delta, b(t)>|  >=  |<p, b(t)>| - R_m(t)

for the partial sum `p`, and the max of that over `t` (floored at zero) lower-bounds the
final objective. Nothing is ever discarded except against this bound, so the result is
the true minimum. If the node budget is exhausted the search reports failure rather than
a number -- the one thing the old method would not do.
"""

from __future__ import annotations

import math

import numpy as np

from pyloki.utils import psr_utils, transforms
from pyloki.utils.misc import C_VAL
from pyloki.utils.snail import MiddleOutScheme


def _basis(tau: np.ndarray, poly_order: int) -> np.ndarray:
    """(n_time, poly_order) matrix of t^k/k!, columns in leaf order [d_kmax .. d_1]."""
    orders = np.arange(poly_order, 0, -1)
    return (tau[:, None] ** orders[None, :]) / _fact(orders)[None, :]


def _fact(k: np.ndarray) -> np.ndarray:
    return np.array([float(math.factorial(int(i))) for i in k])


def build_sets(
    cfg,
    strategy: str,
    poly_order: int,
    f0: float,
    truth_kin: np.ndarray,
    target_stage: int,
    n_seed: int = 1,
) -> tuple[list[np.ndarray], np.ndarray, float]:
    """Exact offset sets whose Minkowski sum is the leaf set at `target_stage`.

    Returns (sets, tau, scale) with every set expressed in the frame at `target_stage`,
    immediately after branching -- the same instant `pruning_multiplicity` measured.
    """
    nseg = int(np.ceil(cfg.nsamps / cfg.bseg_ffa))
    scheme = MiddleOutScheme(nseg, nseg // 2, cfg.tseg_ffa, stride=1)

    dparam = np.array(cfg.get_dparams_actual(cfg.niters_ffa), dtype=float)
    dparam[-1] *= C_VAL / f0

    # Seed: a block of base cells about the truth, as offsets from the truth.
    centre = np.round(truth_kin / dparam) * dparam - truth_kin
    axes = [centre[j] + np.arange(-n_seed, n_seed + 1) * dparam[j] for j in range(poly_order)]
    seed = np.array(np.meshgrid(*axes, indexing="ij")).reshape(poly_order, -1).T
    sets = [np.ascontiguousarray(seed)]

    w = dparam.copy()
    tau = None
    for lvl in range(1, target_stage + 1):
        ref_cur, t_half = scheme.get_current_coord(lvl, moving_grid=True)
        ref_next, _ = scheme.get_coord(lvl)
        delta_t = ref_next - ref_cur

        step = psr_utils.poly_taylor_step_d_vec(
            poly_order, t_half, cfg.nbins, cfg.eta, np.array([f0]), t_ref=0)
        shift = psr_utils.poly_taylor_shift_d_vec(
            w[None, :], step, t_half, cfg.nbins, np.array([f0]), t_ref=0)[0]

        axes = []
        for j in range(poly_order):
            if shift[j] >= cfg.eta - 1e-12:
                n = max(1, int(np.ceil(w[j] / step[0, j] - 1e-12)))
                act = w[j] / n
                axes.append(-w[j] / 2.0 + act / 2.0 + np.arange(n) * act)
                w[j] = act
            else:
                axes.append(np.zeros(1))
        if max(len(a) for a in axes) > 1:
            o = np.array(np.meshgrid(*axes, indexing="ij")).reshape(poly_order, -1).T
            sets.append(np.ascontiguousarray(o))

        tau = np.linspace(delta_t - t_half, delta_t + t_half, 256)
        if lvl == target_stage:
            break
        # Transport every set already collected into the next frame.
        n_p = poly_order + 1
        powers = np.tril(np.arange(n_p)[:, None] - np.arange(n_p))
        t_mat = delta_t**powers / _fact(powers.ravel()).reshape(powers.shape) * np.tril(
            np.ones_like(powers))
        t_kin = t_mat[:poly_order, :poly_order]
        sets = [s @ t_kin.T for s in sets]
        w = transforms.shift_taylor_errors(
            np.concatenate([w, [0.0]])[None, :], delta_t, strategy)[0][:poly_order]

    return sets, tau, f0 / C_VAL


def min_excursion(
    sets: list[np.ndarray],
    tau: np.ndarray,
    scale: float,
    poly_order: int,
    tol: float,
    max_nodes: int = 4_000_000,
    init_best: float = np.inf,
    keep: int = 0,
) -> tuple[float, bool, int] | tuple[float, bool, int, list]:
    """Min over the Minkowski sum of `sets` of max_t |dPhi| / tol.

    Returns `(value, exact, nodes)`. When `exact` is True the search completed and
    `value` is the true minimum. When it is False the node budget ran out and `value` is
    the best leaf actually found, which is still a **valid upper bound** -- every
    incumbent is a real leaf. That one-sided guarantee is enough to order two strategies
    when one side is exact, and it is the honest thing to report when the search cannot
    finish; what it never does is pass off a truncated search as a minimum.
    """
    basis = _basis(tau, poly_order) * (scale / tol)
    proj = [s @ basis.T for s in sets]                    # (m_i, n_time) each
    # Branch on the most decisive sets first. `sets` is permuted identically so that a
    # delta can be reconstructed from the choices.
    order = np.argsort([-p.__abs__().max() for p in proj])
    proj = [proj[i] for i in order]
    offs = [sets[i] for i in order]
    radii = [np.abs(p).max(axis=0) for p in proj]         # r_i(t)
    suffix = [np.zeros_like(radii[0]) for _ in range(len(proj) + 1)]
    for i in range(len(proj) - 1, -1, -1):
        suffix[i] = suffix[i + 1] + radii[i]

    # Seeding the incumbent turns optimisation into a decision problem: with
    # `init_best = v`, the search either returns something strictly below `v`
    # (proving the set beats it) or completes without improving (proving it does
    # not). Both verdicts are sound, and both are far cheaper than the true min.
    best = float(init_best)
    nodes = 0
    found: list = []

    def descend(idx: int, partial: np.ndarray, delta: np.ndarray) -> None:
        nonlocal best, nodes
        if nodes > max_nodes:
            return
        if idx == len(proj):
            val = float(np.abs(partial).max())
            if keep:
                found.append((val, delta.copy()))
                if len(found) > 4 * keep:
                    found.sort(key=lambda z: z[0])
                    del found[keep:]
            best = min(best, val)
            return
        cand = partial[None, :] + proj[idx]
        acand = np.abs(cand)
        lo = np.maximum(acand - suffix[idx + 1][None, :], 0.0).max(axis=1)
        # Sort by the bound (required for the `break` below to be valid), but break the
        # ties -- and at shallow depth the bound is 0 for nearly every child -- on the
        # partial objective. Without the tie-break the first descent is arbitrary, which
        # leaves a useless incumbent and so no pruning at all.
        for k in np.lexsort((acand.max(axis=1), lo)):
            nodes += 1
            if nodes > max_nodes:
                return
            if lo[k] >= best:
                break          # sorted, so every later child is pruned too
            descend(idx + 1, cand[k], delta + offs[idx][k])

    descend(0, np.zeros(len(tau)), np.zeros(sets[0].shape[1]))
    if keep:
        found.sort(key=lambda z: z[0])
        return best, nodes <= max_nodes, nodes, found[:keep]
    return best, nodes <= max_nodes, nodes
