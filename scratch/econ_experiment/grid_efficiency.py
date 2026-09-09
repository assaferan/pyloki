"""Measure what the Appendix D Chebyshev coarsening already buys.

The QUESTIONS.md note that started this thread claimed "additional unexploited
headroom (a coarser grid, or lower k_max, for the same tolerance)". The coarsening
factor is implemented in psr_utils.poly_taylor_step_f as `2**k * dparams_f`, and is
ON BY DEFAULT (use_cheby=True). This quantifies it so the "unexploited" claim can
be checked rather than assumed.
"""

import numpy as np

from pyloki.config import ParamLimits, PulsarSearchConfig

# The config validated in HANDOFF.md section 6.
NSAMPS, DT, PERIOD, NBINS, NSEG, BRUTE_DIV = 2**20, 64e-6, 0.007, 64, 4, 128
ACCEL, JERK, SNAP = 500.0, 6.0, 17.336

freq = 1.0 / PERIOD
tobs = NSAMPS * DT
half = abs(SNAP) * 0.05
p = ParamLimits.from_upper(
    (freq - 1.0, freq + 1.0), [SNAP, JERK, ACCEL], (SNAP - half, SNAP + half), tobs
)
cfg = PulsarSearchConfig(
    nsamps=NSAMPS, tsamp=DT, nbins=NBINS, eta=1, param_limits=p.limits,
    bseg_brute=NSAMPS // BRUTE_DIV, bseg_ffa=NSAMPS // NSEG, prune_poly_order=4,
    ducy_max=0.5, wtsp=1.2, use_fourier=True,
    tiling_strategy="aggressive", branch_max=16,
)

print(f"config: tobs={tobs:.3f}s nbins={NBINS} poly_order=4 nseg={NSEG}")
print(f"param_names: {cfg.param_names}\n")

print("=" * 72)
print("Grid step sizes at the final FFA level (dparams, reverse order)")
print("=" * 72)
lvl = cfg.niters_ffa
for flag in (False, True):
    d = np.asarray(cfg.get_dparams(lvl, flag))
    n = np.asarray(cfg.get_param_grid_count(lvl, flag))
    print(f"  use_cheby={str(flag):5s} dparams={np.array2string(d, precision=4)}")
    print(f"  {'':17s}grid_count={n.tolist()}  product={np.prod(n.astype(float)):.4g}")

d_off = np.asarray(cfg.get_dparams(lvl, False))
d_on = np.asarray(cfg.get_dparams(lvl, True))
print(f"\n  step ratio (on/off) = {np.array2string(d_on / d_off, precision=3)}")
print("  expected 2**k for the k-th derivative, k ascending from the last entry")

print("\n" + "=" * 72)
print("Branching pattern: total leaves explored over the whole prune")
print("=" * 72)
for kind, flag in (
    ("poly_taylor_moving", False),
    ("poly_taylor_moving", True),
    ("poly_chebyshev_moving", True),
):
    bp = np.asarray(cfg.generate_branching_pattern(kind=kind, ref_seg=NSEG // 2,
                                                   use_cheby_coarsening=flag))
    # leaves at stage i = product of branching factors up to i; total = sum over stages
    cum = np.cumprod(bp.astype(float))
    print(f"  {kind:22s} cheby={str(flag):5s} pattern={np.array2string(bp, precision=2)}")
    print(f"  {'':22s} {'':11s} total leaves={cum.sum():.4g}  peak={cum.max():.4g}")

bp_off = np.cumprod(np.asarray(cfg.generate_branching_pattern(
    kind="poly_taylor_moving", ref_seg=NSEG // 2, use_cheby_coarsening=False)).astype(float)).sum()
bp_on = np.cumprod(np.asarray(cfg.generate_branching_pattern(
    kind="poly_taylor_moving", ref_seg=NSEG // 2, use_cheby_coarsening=True)).astype(float)).sum()
print(f"\n  coarsening already saves a factor of {bp_off / bp_on:.4g} in total leaves")
