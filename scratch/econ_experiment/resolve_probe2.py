import pathlib
"""Sup-norm PHASE error of the resolve step: naive truncation vs Chebyshev economization.

Two fixes relative to the earlier version of this script (kept as .bak):

1. The sup-norm comparison was broken by a numpy SHAPE bug (not a units bug, as the
   handoff note claimed). `phase_of_exact_model` was called with a (1,1) array as
   `delta_t`; `shift_taylor_params` broadcasts that silently, so the per-sample result
   carried a stray axis and `phi_exact` came out (401,1) instead of (401,). The final
   `phi_naive - phi_exact` then broadcast to a (401,401) OUTER DIFFERENCE, so
   `np.max(np.abs(...))` returned the full phase *spread* across the segment
   (= freq*tseg/2 after midpoint alignment, hence the nonsensical ~149.8 cycles for
   BOTH methods). Fixed here by evaluating the exact polynomial directly, the same way
   tests/test_transforms.py does, which removes the shift-in-a-loop pattern entirely.

2. Added a sweep over segment half-width, because at this dataset's native scale the
   error is ~1e-7 cycles and the accel grid is degenerate (see the warning printed
   below), so the native-scale result cannot demonstrate anything about grid selection.
"""

import math

import numpy as np

from pyloki.config import ParamLimits, PulsarSearchConfig
from pyloki.ffa import DynamicProgramming
from pyloki.io.timeseries import TimeSeries
from pyloki.utils import transforms
from pyloki.utils.misc import C_VAL

BASE = str(pathlib.Path(__file__).resolve().parent)
d = np.load(f"{BASE}/data.npz")
ts_e, ts_v, dt = d["ts_e"], d["ts_v"], float(d["dt"])
freq_true, nsamps = float(d["freq"]), int(d["nsamps"])
tim_data = TimeSeries(ts_e, ts_v, dt)

# Snap-order test: two orders dropped (snap AND jerk) when resolving onto the 2-param
# (accel, freq) FFA grid -- the regime where Chebyshev economization is mathematically
# expected to differ from naive truncation even for the accel (even-order) coefficient.
accel_true, jerk_true, snap_true = 500.0, 6.0, 20.0
tobs = nsamps * dt

p = ParamLimits.from_upper(
    (freq_true - 1.0, freq_true + 1.0),
    [snap_true, jerk_true, accel_true],
    (-50.0, 50.0),
    tobs,
)
search_cfg = PulsarSearchConfig(
    nsamps=nsamps, tsamp=dt, nbins=64, eta=1, param_limits=p.limits,
    bseg_brute=nsamps // 8, bseg_ffa=nsamps // 2, prune_poly_order=4,
    ducy_max=0.5, wtsp=1.2, use_fourier=True, tiling_strategy="aggressive", branch_max=16,
)
dyp = DynamicProgramming(tim_data, search_cfg)
dyp.initialize()
dyp.execute()

nbins = search_cfg.nbins
bin_cycles = 1.0 / nbins
poly_order = search_cfg.prune_poly_order  # 4 -> [snap, jerk, accel, vel, delay]

print("\n" + "=" * 72)
print("Real dyp-derived scales")
print("=" * 72)
print(f"param_names            : {search_cfg.param_names}")
print(f"param_grid_count       : {dyp.param_grid_count}")
print(f"tseg (post-execute)    : {dyp.tseg} s")
print(f"tobs                   : {tobs} s")
print(f"nbins                  : {nbins}  (1 phase bin = {bin_cycles:.6g} cycles)")
if int(np.asarray(dyp.param_grid_count).ravel()[2]) == 1:
    print(
        "\nWARNING: the accel grid has only ONE point at this configuration, so the\n"
        "         resolved accel cell cannot differ between methods no matter how\n"
        "         large the coefficient difference is. Phase error below is therefore\n"
        "         a statement about approximation quality, NOT about grid selection."
    )

# True Taylor vector, referenced to coord_cur: [snap, jerk, accel, vel, delay]
d_true = np.array([snap_true, jerk_true, accel_true, 0.0, 0.0])
leaves_batch = np.zeros((1, poly_order + 2, 2))
leaves_batch[0, :-1, 0] = d_true
leaves_batch[0, -1, 0] = freq_true


def resolve(d_vec, freq, tseg, *, economized):
    """The resolve step's accel/freq output, with economization on or off.

    Mirrors poly_taylor_resolve_batch (src/pyloki/core/taylor.py:222-227) with the
    same coord layout the pruning loop uses: cur and init share an origin, add is one
    segment later, half-width = tseg/2.
    """
    t0_cur = t0_init = 0.5 * tseg
    t0_add, half_width_add = 1.5 * tseg, 0.5 * tseg
    param_vec_batch = np.asarray(d_vec, dtype=np.float64)[np.newaxis, :]
    f0_batch = np.array([freq], dtype=np.float64)

    dvec_t_add = transforms.shift_taylor_params(param_vec_batch, t0_add - t0_cur)
    if economized:
        dvec_t_add = transforms.economize_taylor_params(
            dvec_t_add, half_width_add, n_keep=3
        )
    dvec_t_init = transforms.shift_taylor_params(param_vec_batch, t0_init - t0_cur)
    accel_new = dvec_t_add[:, -3]
    vel_new = dvec_t_add[:, -2] - dvec_t_init[:, -2]
    freq_new = f0_batch * (1 - vel_new / C_VAL)
    return float(accel_new[0]), float(freq_new[0])


def eval_taylor(d_vec, t):
    """Evaluate sum_k d_vec[-(k+1)]/k! * t**k for d_vec ordered [d_kmax, ..., d_0]."""
    n = len(d_vec)
    out = np.zeros_like(t, dtype=np.float64)
    for k in range(n):
        out = out + d_vec[n - 1 - k] / math.factorial(k) * t**k
    return out


def sup_norm_phase_error(d_vec, freq, tseg, accel_new, freq_new, n_grid=2001):
    """max |Phi_model - Phi_exact| in cycles over the added segment.

    Phi_exact uses the full (snap, jerk, accel) polynomial about coord_cur; the model
    uses the resolved constant-accel form about coord_add. A common DC offset is
    removed at the segment midpoint (it is absorbed by the tree's phase bookkeeping
    and carries no sensitivity cost).
    """
    half = 0.5 * tseg
    t_add = np.linspace(-half, half, n_grid)
    t_cur = t_add + tseg  # same absolute times, referenced to coord_cur

    phi_exact = freq * t_cur - freq * eval_taylor(np.asarray(d_vec), t_cur) / C_VAL
    phi_model = freq_new * t_add - freq_new * (0.5 * accel_new * t_add**2) / C_VAL

    mid = n_grid // 2
    phi_model = phi_model - (phi_model[mid] - phi_exact[mid])
    return float(np.max(np.abs(phi_model - phi_exact)))


# --- Native scale -------------------------------------------------------------
tseg_native = float(dyp.tseg)
a_n, f_n = resolve(d_true, freq_true, tseg_native, economized=False)
a_e, f_e = resolve(d_true, freq_true, tseg_native, economized=True)

print("\n" + "=" * 72)
print(f"Resolved coefficients at native tseg = {tseg_native:g} s")
print("=" * 72)
print(f"TRUE  accel={accel_true:<14g} freq={freq_true:.10f}")
print(f"naive accel={a_n:<14.6f} freq={f_n:.10f}")
print(f"econ  accel={a_e:<14.6f} freq={f_e:.10f}")

e_n = sup_norm_phase_error(d_true, freq_true, tseg_native, a_n, f_n)
e_e = sup_norm_phase_error(d_true, freq_true, tseg_native, a_e, f_e)
print(f"\nsup-norm phase error: naive={e_n:.6e}  econ={e_e:.6e} cycles")
print(f"                      ratio econ/naive = {e_e / e_n:.4f}")
print(f"                      naive = {e_n / bin_cycles:.3e} phase bins  <- negligible")

# --- Sweep over segment half-width -------------------------------------------
# The true (snap, jerk, accel) vector is held FIXED while tseg is stretched. This is a
# mathematical sensitivity sweep: it locates the segment duration at which the dropped
# high-order content first costs real sensitivity. It is NOT a claim that this dataset
# (tobs = 2.1 s) or these parameter values are physical at the longer durations -- see
# the closing note about what snap amplitude the crossing implies.
print("\n" + "=" * 72)
print("Sup-norm phase error vs segment duration (true d_vec held fixed)")
print("=" * 72)
print(
    f"{'tseg [s]':>12} {'naive [cyc]':>14} {'econ [cyc]':>14} "
    f"{'econ/naive':>11} {'naive [bins]':>13}"
)
print("-" * 72)

tsegs = tseg_native * 2.0 ** np.arange(0, 15)
cross_bin = cross_tenth = None
valid_ratios = []
for tseg in tsegs:
    a_n, f_n = resolve(d_true, freq_true, tseg, economized=False)
    a_e, f_e = resolve(d_true, freq_true, tseg, economized=True)
    e_n = sup_norm_phase_error(d_true, freq_true, tseg, a_n, f_n)
    e_e = sup_norm_phase_error(d_true, freq_true, tseg, a_e, f_e)
    ratio = e_e / e_n if e_n > 0 else float("nan")
    flag = ""
    if e_n < 1.0:
        valid_ratios.append(ratio)
    else:
        flag = "  (model invalid: err >> 1 cyc, ratio meaningless)"
    if cross_bin is None and e_n >= bin_cycles:
        cross_bin = tseg
        flag = flag or "  <- naive crosses 1 phase bin"
    if cross_tenth is None and e_n >= 0.1:
        cross_tenth = tseg
        flag = flag or "  <- naive crosses 0.1 cycles"
    print(
        f"{tseg:>12.4g} {e_n:>14.4e} {e_e:>14.4e} {ratio:>11.4f} "
        f"{e_n / bin_cycles:>13.3e}{flag}"
    )

print("\n" + "=" * 72)
print("Summary")
print("=" * 72)
print(
    f"economization/naive sup-norm ratio = {np.mean(valid_ratios):.4f} "
    f"(spread {min(valid_ratios):.4f}-{max(valid_ratios):.4f})\n"
    f"  over the {len(valid_ratios)} rows where the constant-accel model is still\n"
    f"  valid (naive error < 1 cycle). Rows above that are numerical nonsense\n"
    f"  and are excluded from this average."
)

if cross_bin is not None:
    print(
        f"\nnaive error reaches 1 phase bin (1/{nbins} cyc) at tseg ~ {cross_bin:.4g} s,"
        f"\n  for the injected (snap, jerk, accel) = "
        f"({snap_true:g}, {jerk_true:g}, {accel_true:g})."
    )

    def naive_err(snap, tseg):
        d_vec = np.array([snap, jerk_true, accel_true, 0.0, 0.0])
        a, f = resolve(d_vec, freq_true, tseg, economized=False)
        return sup_norm_phase_error(d_vec, freq_true, tseg, a, f)

    # Solved numerically rather than by inverting an assumed error ~ snap * t_s**4:
    # the measured scaling is ~t_s**3.9 because the jerk term enters at t_s**3 and
    # jerk is held fixed here, so a pure t_s**-4 inversion would be wrong.
    print("\n  snap amplitude needed to cost 1 phase bin at a realistic tseg")
    print("  (jerk and accel held at their injected values, solved by bisection):")
    for tseg_real in (10.0, 100.0, 1000.0):
        if naive_err(0.0, tseg_real) >= bin_cycles:
            print(f"    tseg = {tseg_real:>7g} s : >= 1 bin already at snap = 0 "
                  f"(the fixed jerk alone suffices)")
            continue
        lo, hi = 0.0, 1e-12
        while naive_err(hi, tseg_real) < bin_cycles and hi < 1e12:
            lo, hi = hi, hi * 4.0
        if hi >= 1e12:
            print(f"    tseg = {tseg_real:>7g} s : no solution below 1e12")
            continue
        for _ in range(100):
            mid = 0.5 * (lo + hi)
            if naive_err(mid, tseg_real) >= bin_cycles:
                hi = mid
            else:
                lo = mid
        print(f"    tseg = {tseg_real:>7g} s : snap ~ {hi:.4g} m/s^4")
else:
    print("naive error never reaches 1 phase bin over the swept range")
