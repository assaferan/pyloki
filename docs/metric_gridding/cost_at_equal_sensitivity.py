"""Cost at equal worst-case leaf mismatch, against BOTH box baselines.

`aggressive` is cheap but its recorded region under-reports the true one after a shear
-- the coverage gap the project targets. `conservative` tracks the honest AABB and is
gap-free but expensive. The tiling dilemma is exactly the distance between them, so the
metric strategy has to be placed against both, not against one labelled "the box".
"""

from __future__ import annotations

import numpy as np

from pyloki.core import metric, taylor
from pyloki.utils import psr_utils, transforms
from pyloki.utils.snail import MiddleOutScheme

F0 = 1.0 / 0.007


def box_reference(cfg, nseg, ref, po, strategy):
    """prod B(s) and worst leaf mismatch for a box strategy.

    Replicates `generate_bp_poly_taylor`'s per-leaf loop so the leaf's own dparam can be
    priced in the metric; validated to reproduce `generate_branching_pattern` exactly
    for `aggressive`.
    """
    sch = MiddleOutScheme(nseg, ref, cfg.tseg_ffa, stride=1)
    dp = np.array(cfg.get_dparams_actual(cfg.niters_ffa), dtype=float)
    dp[-1] *= 3e8 / F0
    vec = np.zeros(po + 1)
    vec[:po] = dp
    nb, worst = 1.0, 0.0
    for lvl in range(1, nseg):
        rc, th = sch.get_current_coord(lvl, moving_grid=True)
        rn, _ = sch.get_coord(lvl)
        dt = rn - rc
        g = metric.poly_phase_metric(0.0, -th, th, po, F0, cfg.nbins, cfg.ducy_max)
        h = vec[:po] / 2
        worst = max(worst, float(h @ g @ h))
        dn = psr_utils.poly_taylor_step_d_vec(
            po, th, cfg.nbins, cfg.eta, np.array([F0]), t_ref=0)
        sh = psr_utils.poly_taylor_shift_d_vec(
            vec[None, :po], dn, th, cfg.nbins, np.array([F0]), t_ref=0)[0]
        npts = np.array([
            max(1, int(np.ceil(vec[j] / dn[0, j] - 1e-12)))
            if sh[j] >= cfg.eta - 1e-12 else 1
            for j in range(po)], dtype=float)
        nb *= np.prod(npts)
        vec[:po] /= npts
        vec = transforms.shift_taylor_errors(vec[None, :], dt, strategy)[0]
    return nb, worst


def metric_cost(cfg_maker, nseg, ref, po, target_m, defers=(1.5, 2, 3, 4, 6, 8, 12, 16)):
    """Cheapest prod B(s) reaching `target_m` worst-case, over the deferral factor."""
    best = None
    for R in defers:
        m_max = target_m / R**2
        c = cfg_maker(m_max, float(R))
        seed = c.metric_seed_region(
            c.get_dparams_actual(c.niters_ffa, use_cheby_coarsening=False),
            c.get_param_arr(c.get_dparams(c.niters_ffa, use_cheby_coarsening=False)))
        sch = MiddleOutScheme(nseg, ref, c.tseg_ffa, stride=1)
        region, tot, ev = seed, 1.0, 0
        try:
            for lvl in range(1, nseg):
                rc, th = sch.get_current_coord(lvl, moving_grid=True)
                rn, _ = sch.get_coord(lvl)
                dt = rn - rc
                off, _, region = taylor.metric_branch_tables(
                    region, th, dt, c.nbins, c.metric_ducy, po, m_max,
                    c.metric_branch_max, float(R))
                tot *= len(off)
                ev += len(off) > 1
                region, _ = taylor.metric_transform_region(region, dt, po)
        except ValueError:
            continue
        if best is None or tot < best[1]:
            best = (R, tot, ev)
    return best
