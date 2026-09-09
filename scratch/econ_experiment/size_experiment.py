"""Screen a candidate poly_order=4 search config for the preconditions an
end-to-end SNR comparison needs, WITHOUT paying for the pruning stage.

An SNR difference between naive truncation and Chebyshev economization is only
possible if all of these hold:
  1. poly_order >= 4, so two orders are dropped and both parities are spanned
     (at poly_order=3 the accel cell is provably identical -- see
     economize_taylor_params).
  2. The accel grid has more than one point, or the resolved cell cannot differ.
  3. The naive resolve error is a meaningful fraction of a phase bin at the real
     segment scale, half_width_add = 0.5 * tseg_ffa (snail.py:249).
  4. nsegments is moderate: it sets the number of pruning stages.

Usage: python size_experiment.py [--snap S] [--nseg N] [--dt DT] [--nsamps-pow P]
"""

import argparse

import numpy as np

import econ_metrics as em
from pyloki.config import ParamLimits, PulsarSearchConfig
from pyloki.ffa import DynamicProgramming
from pyloki.simulation.pulse import PulseSignalConfig

ap = argparse.ArgumentParser()
ap.add_argument("--period", type=float, default=0.007)
ap.add_argument("--dt", type=float, default=64e-6)
ap.add_argument("--nsamps-pow", type=int, default=15)
ap.add_argument("--accel", type=float, default=500.0)
ap.add_argument("--jerk", type=float, default=6.0)
ap.add_argument("--snap", type=float, default=0.0)
ap.add_argument("--target-bins", type=float, default=0.0,
                help="if set, solve snap so naive error = this many phase bins")
ap.add_argument("--nseg", type=int, default=16, help="target pruning stages")
ap.add_argument("--brute-div", type=int, default=0, help="bseg_brute = nsamps//this")
ap.add_argument("--nbins", type=int, default=64)
ap.add_argument("--accel-widen", type=float, default=1.0,
                help="scale the accel search half-range about its centre.\n                      Decouples the accel grid width (verdict 2) from snap\n                      (verdict 5): from_upper derives the accel range from\n                      the snap bracket, which ties the two together.")
ap.add_argument("--snap-frac", type=float, default=0.05,
                help="snap bracket half-width as a fraction of the true snap; the\n                      FFA init requires every non-freq param to fit ONE cell\n                      (ffa.py:429), so this bracket is centred on the true snap")
ap.add_argument("--snr", type=float, default=12.0)
args = ap.parse_args()

nsamps = 2**args.nsamps_pow
tobs = nsamps * args.dt
freq = 1.0 / args.period
nbins = args.nbins
bin_cycles = 1.0 / nbins

# tseg_ffa = tobs / nseg  ->  bseg_ffa = nsamps // nseg
bseg_ffa = nsamps // args.nseg
bseg_brute = nsamps // (args.brute_div or args.nseg * 4)  # must be < bseg_ffa

tseg_ffa_pre = bseg_ffa * args.dt
if args.target_bins > 0:
    args.snap = em.snap_for_phase_budget(
        args.jerk, args.accel, freq, tseg_ffa_pre, args.target_bins / nbins
    )

print("=" * 74)
print("Candidate config")
print("=" * 74)
print(f"period={args.period}s  freq={freq:.4f}Hz  dt={args.dt:g}s  nsamps={nsamps}")
print(f"tobs={tobs:.4f}s  nbins={nbins}  (1 bin = {bin_cycles:.6g} cycles)")
print(f"injected: accel={args.accel:g} jerk={args.jerk:g} snap={args.snap:g}")
print(f"bseg_ffa={bseg_ffa} -> tseg_ffa={bseg_ffa * args.dt:.4f}s   bseg_brute={bseg_brute}")

# --- Precondition 3, checked analytically first (free) ----------------------
tseg_ffa = bseg_ffa * args.dt
d_true = np.array([args.snap, args.jerk, args.accel, 0.0, 0.0])
e_n, e_e = em.phase_error_pair(d_true, freq, tseg_ffa)
snap_1bin = em.snap_for_phase_budget(args.jerk, args.accel, freq, tseg_ffa, bin_cycles)
print("\n" + "=" * 74)
print("Precondition 3: resolve error at the real segment scale")
print("=" * 74)
print(f"naive sup-norm phase error = {e_n:.6e} cycles = {e_n / bin_cycles:.4g} bins")
print(f"econ  sup-norm phase error = {e_e:.6e} cycles = {e_e / bin_cycles:.4g} bins")
if e_n > 0:
    print(f"econ/naive ratio           = {e_e / e_n:.4f}")
print(f"snap needed for 1 full bin at this tseg_ffa = {snap_1bin:.6g} m/s^4")
verdict3 = e_n >= 0.25 * bin_cycles
print(f"VERDICT 3: {'PASS' if verdict3 else 'FAIL'} "
      f"(want naive error >= 0.25 bin; got {e_n / bin_cycles:.4g})")

# --- Build the search and check preconditions 2 and 4 -----------------------
mod_kwargs = {"acc": args.accel, "jerk": args.jerk, "snap": args.snap}
cfg = PulseSignalConfig(
    period=args.period, dt=args.dt, nsamps=nsamps, snr=args.snr,
    ducy=0.1, mod_kwargs=mod_kwargs,
)
tim_data = cfg.generate(shape="gaussian")

half = max(abs(args.snap) * args.snap_frac, 1e-9)
snap_lo, snap_hi = args.snap - half, args.snap + half
p = ParamLimits.from_upper(
    (freq - 1.0, freq + 1.0),
    [args.snap, args.jerk, args.accel],
    (snap_lo, snap_hi),
    tobs,
)
print(f"snap bracket             : ({snap_lo:.6g}, {snap_hi:.6g})")

limits = np.array(p.limits, dtype=np.float64, copy=True)
if args.accel_widen != 1.0:
    lo, hi = limits[2]
    mid, halfw = 0.5 * (lo + hi), 0.5 * (hi - lo) * args.accel_widen
    limits[2] = (mid - halfw, mid + halfw)
    print(f"accel range widened x{args.accel_widen:g} -> "
          f"({limits[2, 0]:.6g}, {limits[2, 1]:.6g})")
search_cfg = PulsarSearchConfig(
    nsamps=nsamps, tsamp=args.dt, nbins=nbins, eta=1, param_limits=limits,
    bseg_brute=bseg_brute, bseg_ffa=bseg_ffa, prune_poly_order=4,
    ducy_max=0.5, wtsp=1.2, use_fourier=True,
    tiling_strategy="aggressive", branch_max=16,
)
dyp = DynamicProgramming(tim_data, search_cfg)
try:
    dyp.initialize()
except ValueError as exc:
    print(f"\nFFA INIT REJECTED: {exc}")
    print("  -> narrow --snap-frac, or lower --nseg / --nsamps-pow")
    print("\nOVERALL: NOT READY")
    raise SystemExit(1) from None
nseg_pre = dyp.nsegments
tseg_pre = dyp.tseg
dyp.execute()
nseg_actual = dyp.nsegments

print("\n" + "=" * 74)
print("Preconditions 2 and 4: grid and stage count")
print("=" * 74)
print(f"param_names              : {search_cfg.param_names}")
print(f"param_grid_count         : {dyp.param_grid_count}")
print(f"dyp.tseg (post-execute)  : {dyp.tseg:.6f}s  (expected tseg_ffa={tseg_ffa:.6f})")
print(f"nsegments  pre-execute   : {nseg_pre}   (tseg={tseg_pre:.6f}s = tseg_brute)")
print(f"nsegments post-execute   : {nseg_actual}  (pruning stages)")
print(f"total FFA grid sizes     : {np.asarray(dyp.param_grid_count).prod():.4g} cells per segment")
print(f"param_limits             : {limits.tolist()}")

vc = em.max_v_over_c(d_true, tobs)
print("\n" + "=" * 74)
print("Precondition 5: first-order Doppler validity")
print("=" * 74)
print(f"peak |v|/c over the observation = {vc:.4g}")
print(f"  implied 2nd-order phase error  ~ {freq * tobs * vc**2:.4g} cycles "
      f"({freq * tobs * vc**2 * nbins:.4g} bins)")
# The bar that matters is relative: the neglected 2nd-order Doppler term must be
# well below the 1st-order truncation error we are trying to measure, or it swamps
# the very effect under test. This ratio grows as nseg**7 / (freq * tseg), so many
# short segments are far worse than few long ones.
second = freq * tobs * vc**2
ratio5 = second / e_n if e_n > 0 else float("inf")
print(f"  2nd-order / naive resolve error = {ratio5:.4g}")
verdict5 = ratio5 <= 0.2
print(f"VERDICT 5: {'PASS' if verdict5 else 'FAIL'} "
      f"(want 2nd-order <= 0.2x the effect under test)")

n_accel = int(np.asarray(dyp.param_grid_count).ravel()[2])
verdict2 = n_accel > 1
verdict4 = 4 <= nseg_actual <= 128
print(f"\nVERDICT 2: {'PASS' if verdict2 else 'FAIL'} (accel grid points = {n_accel}, want > 1)")
print(f"VERDICT 4: {'PASS' if verdict4 else 'FAIL'} (nsegments = {nseg_actual}, want 4-128)")

print("\n" + "=" * 74)
allv = verdict2 and verdict3 and verdict4 and verdict5
print(f"OVERALL: {'ALL PRECONDITIONS MET -- worth running the pruning' if allv else 'NOT READY'}")
print("=" * 74)
