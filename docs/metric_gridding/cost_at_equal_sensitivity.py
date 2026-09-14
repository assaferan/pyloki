"""Cost at equal sensitivity: prod B(s) and the worst mismatch a leaf actually carries."""
import numpy as np
from pyloki.core import metric, taylor
from pyloki.utils.snail import MiddleOutScheme
from pyloki.config import ParamLimits, PulsarSearchConfig

F0 = 1 / 0.007
nsamps, tsamp, bseg = 2**18, 64e-6, 2**13
tobs = nsamps * tsamp
nseg = nsamps // bseg
lim = ParamLimits.from_upper((F0*(1-1e-3), F0*(1+1e-3)), [6.0, 500.0], (-8.0, 8.0), tobs)

def cfg_for(m_max):
    return PulsarSearchConfig(
        nsamps=nsamps, tsamp=tsamp, nbins=64, eta=1.0, param_limits=lim.limits,
        bseg_brute=bseg//8, bseg_ffa=bseg, prune_poly_order=3, ducy_max=0.2,
        wtsp=1.5, use_fourier=True, tiling_strategy="metric",
        branch_max=16, m_max=m_max)

def run(m_max, defer):
    cfg = cfg_for(m_max)
    seed = cfg.metric_seed_region(
        cfg.get_dparams_actual(cfg.niters_ffa, use_cheby_coarsening=False),
        cfg.get_param_arr(cfg.get_dparams(cfg.niters_ffa, use_cheby_coarsening=False)))
    sch = MiddleOutScheme(nseg, nseg//2, cfg.tseg_ffa, stride=1)
    region, tot, worst, nb = seed, 1.0, 0.0, 0
    for lvl in range(1, nseg):
        rc, th = sch.get_current_coord(lvl, moving_grid=True)
        rn, _ = sch.get_coord(lvl); dt = rn - rc
        g = metric.poly_phase_metric(0.0, dt-th, dt+th, 3, 1.0, cfg.nbins, cfg.metric_ducy)
        oh = metric.region_overhang(region, g, m_max)
        # the mismatch this leaf actually carries at this level
        worst = max(worst, min(oh, defer)**2 * m_max)
        off, _, region = taylor.metric_branch_tables(
            region, th, dt, cfg.nbins, cfg.metric_ducy, 3, m_max, 5_000_000, defer)
        if len(off) > 1:
            nb += 1
        tot *= len(off)
        region, _ = taylor.metric_transform_region(region, dt, 3)
    return tot, worst, nb

print("box strategy reference: prod B(s) = 27, worst leaf mismatch = 0.024\n")
print(f"{'m_max':>9} {'R':>5} {'prod B(s)':>12} {'worst m':>10} {'events':>7}")
for m_max in (0.2, 0.05, 0.02, 0.005, 0.002):
    for defer in (1.0, 1.5, 2.0, 3.0):
        try:
            tot, worst, nb = run(m_max, defer)
        except ValueError as e:
            print(f"{m_max:>9g} {defer:>5g} {'(cap hit)':>12} {'':>10} {'':>7}")
            continue
        flag = "  <-- <= box on both" if tot <= 27 and worst <= 0.024 else ""
        print(f"{m_max:>9g} {defer:>5g} {tot:>12.4g} {worst:>10.4g} {nb:>7d}{flag}")
