"""Figure 7 analogue: cumulative prod B(s) per stage, metric vs the box strategies."""
import sys
import numpy as np
from pyloki.config import ParamLimits, PulsarSearchConfig

PERIOD = 0.007
F0 = 1.0 / PERIOD

def make_cfg(strategy, nsamps, tsamp, bseg_ffa, poly_order, nbins=64,
             m_max=0.2, ducy_max=0.2, freq_pad=1e-3):
    tobs = nsamps * tsamp
    limits = ParamLimits.from_upper(
        (F0 * (1 - freq_pad), F0 * (1 + freq_pad)), [6.0, 500.0], (-8.0, 8.0), tobs,
    )
    return PulsarSearchConfig(
        nsamps=nsamps, tsamp=tsamp, nbins=nbins, eta=1.0,
        param_limits=limits.limits,
        bseg_brute=bseg_ffa // 8, bseg_ffa=bseg_ffa,
        prune_poly_order=poly_order, ducy_max=ducy_max, wtsp=1.5,
        use_fourier=True, tiling_strategy=strategy, branch_max=16, m_max=m_max,
    )

def run(nsamps, tsamp, bseg_ffa, poly_order, m_max=0.2):
    tobs = nsamps * tsamp
    nseg = int(np.ceil(nsamps / bseg_ffa))
    ref = nseg // 2
    print(f"\n=== T_obs={tobs:.1f}s  {nseg} segments of {tobs/nseg:.2f}s  "
          f"poly_order={poly_order}  m_max={m_max} ===")
    out = {}
    for strat in ("aggressive", "quadrature", "conservative", "metric"):
        cfg = make_cfg(strat, nsamps, tsamp, bseg_ffa, poly_order, m_max=m_max)
        bp = cfg.generate_branching_pattern(kind="poly_taylor_moving", ref_seg=ref)
        out[strat] = np.asarray(bp)
    # like-for-like: the box strategies without the 2**k Chebyshev coarsening (D5)
    cfg = make_cfg("aggressive", nsamps, tsamp, bseg_ffa, poly_order, m_max=m_max)
    out["aggressive (no cheby)"] = np.asarray(cfg.generate_branching_pattern(
        kind="poly_taylor_moving", ref_seg=ref, use_cheby_coarsening=False))

    print(f"{'strategy':>22} | {'total prod B(s)':>16} | {'first refine':>12} | B(s) head")
    for k, bp in out.items():
        tot = float(np.prod(bp))
        first = next((i + 1 for i, c in enumerate(bp) if c > 1.0 + 1e-9), None)
        head = np.array2string(bp[:10], precision=1, max_line_width=200)
        print(f"{k:>22} | {tot:16.4g} | {str(first):>12} | {head}")
    return out

if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "small"
    if which == "small":
        run(2**18, 64e-6, 2**13, 3)     # 16.8 s, 32 segments
    elif which == "mid":
        run(2**21, 64e-6, 2**14, 3)     # 134 s, 128 segments
    elif which == "full":
        run(2**24, 64e-6, 2**17, 3)     # 1074 s (17.9 min), 128 segments
