"""Shared resolve-step metrics: naive truncation vs Chebyshev economization.

Mirrors poly_taylor_resolve_batch (src/pyloki/core/taylor.py) so a candidate search
config can be screened for the preconditions an end-to-end SNR comparison needs,
without paying for the pruning stage.
"""

import math

import numpy as np

from pyloki.utils import transforms
from pyloki.utils.misc import C_VAL


def resolve(d_vec, freq, tseg, *, economized):
    """The resolve step's (accel, freq) output, with economization on or off.

    Coord layout matches the pruning loop: cur and init share an origin, the added
    segment is one segment later, half-width = tseg/2 (snail.py:249).
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

    A common DC offset is removed at the segment midpoint: it is absorbed by the
    tree's phase bookkeeping and carries no sensitivity cost.
    """
    half = 0.5 * tseg
    t_add = np.linspace(-half, half, n_grid)
    t_cur = t_add + tseg  # same absolute times, referenced to coord_cur

    phi_exact = freq * t_cur - freq * eval_taylor(np.asarray(d_vec), t_cur) / C_VAL
    phi_model = freq_new * t_add - freq_new * (0.5 * accel_new * t_add**2) / C_VAL

    mid = n_grid // 2
    phi_model = phi_model - (phi_model[mid] - phi_exact[mid])
    return float(np.max(np.abs(phi_model - phi_exact)))


def phase_error_pair(d_vec, freq, tseg):
    """(naive, econ) sup-norm phase error in cycles for one segment scale."""
    a_n, f_n = resolve(d_vec, freq, tseg, economized=False)
    a_e, f_e = resolve(d_vec, freq, tseg, economized=True)
    return (
        sup_norm_phase_error(d_vec, freq, tseg, a_n, f_n),
        sup_norm_phase_error(d_vec, freq, tseg, a_e, f_e),
    )


def snap_for_phase_budget(jerk, accel, freq, tseg, budget_cycles):
    """Smallest snap whose NAIVE resolve error reaches budget_cycles, by bisection.

    Solved numerically rather than by inverting an assumed snap * tseg**4 law: jerk
    enters at tseg**3 and is held fixed, so a pure power-law inversion is wrong.
    """
    def err(snap):
        d_vec = np.array([snap, jerk, accel, 0.0, 0.0])
        a, f = resolve(d_vec, freq, tseg, economized=False)
        return sup_norm_phase_error(d_vec, freq, tseg, a, f)

    if err(0.0) >= budget_cycles:
        return 0.0
    lo, hi = 0.0, 1.0
    while err(hi) < budget_cycles:
        lo, hi = hi, hi * 4.0
        if hi > 1e18:
            return float("inf")
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if err(mid) >= budget_cycles:
            hi = mid
        else:
            lo = mid
    return hi


def max_v_over_c(d_vec, tobs):
    """Peak |velocity|/c over the observation for a Taylor d_vec [snap,jerk,accel,...].

    pyloki models the Doppler shift to first order only (f = f0 * (1 - v/c), see
    poly_taylor_resolve_batch), so a config is only physically meaningful while this
    stays small. Raising snap to make the resolve error reach a phase bin also raises
    v/c as snap * tobs**3, so the two demands pull against each other.
    """
    C = 299792458.0
    d = np.asarray(d_vec, dtype=np.float64)
    # velocity = d/dt of the position polynomial
    t = np.linspace(-0.5 * tobs, 0.5 * tobs, 4001)
    n = len(d)
    v = np.zeros_like(t)
    for k in range(1, n):  # derivative of the t**k term
        coeff = d[n - 1 - k] / math.factorial(k)
        v = v + coeff * k * t ** (k - 1)
    return float(np.max(np.abs(v)) / C)
