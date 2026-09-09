"""Parameter-space metric (mismatch tensor) for the kinematic Taylor basis.

Phase 1 of `metric_PLAN.md`: standalone and pure NumPy. Nothing here is wired into
the search yet, and no existing behaviour is affected. Numba comes in Phase 4.

Conventions (see `docs/metric_gridding/DECISIONS.md`; all observed from the code):

- **Axis ordering (C1).** The branchable axes are the kinematic Taylor coefficients in
  reverse order, `[d_kmax, ..., d_2, d_1]`, i.e. `poly_order` axes carrying derivative
  orders `k = poly_order .. 1`. Index `i` holds order `k = poly_order - i`. `d_0` is
  deliberately absent: it contributes a constant phase offset, which folding absorbs,
  and it is never branched on (`poly_taylor_seed`: "we never branch on d0").
- **Units (C5).** `d_k` in m/s^k, `f0` in Hz, `C_VAL` in m/s.
- **Phase (C6).** In *cycles*, and the code's delay-to-phase relation is
  `phase = f0 * d / C_VAL` with no 2*pi. Here `Phi(t) = f0 * [(t - t_ref) - d(t)/c]`,
  so `dPhi/dd_k = -(f0/c) * (t - t_ref)^k / k!`.
- **Mismatch (D6).** `m = delta^T g delta` is the fractional **amplitude** (S/N) loss,
  `m = 1 - A/A_0`. The `2*pi**2` is folded into `g`. If the power convention is wanted
  instead, this constant becomes `4*pi**2` and every `m_max` halves.
"""

from __future__ import annotations

import math
from typing import NamedTuple

import numpy as np

from pyloki.utils.misc import C_VAL

__all__ = [
    "MMaxBridge",
    "cholesky_factor",
    "ellipsoid_axis_extents",
    "m_max_from_eta",
    "mismatch",
    "poly_phase_metric",
    "shift_matrix",
    "transform_metric",
]

# D6: m is the fractional amplitude (S/N) loss, so g carries 2*pi**2.
# Power convention would be 4*pi**2.
_AMPLITUDE_SCALE = 2.0 * np.pi**2


def _axis_orders(poly_order: int) -> np.ndarray:
    """Derivative order carried by each axis, in leaf order (C1).

    `[poly_order, ..., 2, 1]` — reverse `k`, and `d_0` excluded.
    """
    return np.arange(poly_order, 0, -1, dtype=np.int64)


def _tau_moment(n: int, t_ref: float, t_start: float, t_end: float) -> float:
    """Time average of `(t - t_ref)**n` over `[t_start, t_end]`.

    `<tau^n> = (v^(n+1) - u^(n+1)) / ((n+1) * (v - u))` with `u = t_start - t_ref`
    and `v = t_end - t_ref`. Exact, no quadrature.
    """
    u = t_start - t_ref
    v = t_end - t_ref
    if v == u:
        msg = "t_start and t_end must differ"
        raise ValueError(msg)
    return (v ** (n + 1) - u ** (n + 1)) / ((n + 1) * (v - u))


def poly_phase_metric(
    t_ref: float,
    t_start: float,
    t_end: float,
    poly_order: int,
    f0: float,
    nbins: int | None = None,  # noqa: ARG001  (D7: unused; kept for the plan's API)
) -> np.ndarray:
    """Mismatch tensor `g` for the kinematic Taylor basis about `t_ref`.

    Parameters
    ----------
    t_ref
        Expansion epoch of the Taylor coefficients.
    t_start, t_end
        Absolute interval over which the phase is averaged.
    poly_order
        Number of branchable axes; orders `k = poly_order .. 1` (C1).
    f0
        Spin frequency in Hz.
    nbins
        Unused. `g` does not depend on it; it enters only via `m_max_from_eta`.
        Accepted because the plan's API specifies it (D7).

    Returns
    -------
    np.ndarray
        `(poly_order, poly_order)` symmetric positive-definite matrix, ordered to match
        leaf rows `[:-2]`, such that `delta^T g delta` is the fractional amplitude loss.

    Notes
    -----
    `g_ij = 2*pi**2 * (f0/c)**2 * Cov(tau^ki / ki!, tau^kj / kj!)`, the covariance being
    a time average over the interval. Subtracting the product of means projects out the
    constant-phase mode, which is the unobservable direction.
    """
    if poly_order < 1:
        msg = f"poly_order must be >= 1, got {poly_order}"
        raise ValueError(msg)
    orders = _axis_orders(poly_order)
    inv_fact = np.array([1.0 / math.factorial(int(k)) for k in orders])

    # First moments <tau^k / k!> per axis, and cross moments <tau^(ki+kj)>.
    mean = np.array(
        [_tau_moment(int(k), t_ref, t_start, t_end) for k in orders]
    ) * inv_fact
    cross = np.empty((poly_order, poly_order), dtype=np.float64)
    for i in range(poly_order):
        for j in range(poly_order):
            n = int(orders[i] + orders[j])
            cross[i, j] = _tau_moment(n, t_ref, t_start, t_end)
    cross = cross * inv_fact[:, None] * inv_fact[None, :]

    cov = cross - np.outer(mean, mean)
    g = _AMPLITUDE_SCALE * (f0 / C_VAL) ** 2 * cov
    # Symmetrise against round-off so Cholesky is reliable.
    return 0.5 * (g + g.T)


def mismatch(g: np.ndarray, delta: np.ndarray) -> float | np.ndarray:
    """`m = delta^T g delta`, the fractional amplitude loss (D6).

    `delta` may be `(n_params,)` or a batch `(..., n_params)`; the return shape is the
    leading shape of `delta`.
    """
    delta = np.asarray(delta, dtype=np.float64)
    out = np.einsum("...i,ij,...j->...", delta, g, delta)
    return float(out) if out.ndim == 0 else out


def ellipsoid_axis_extents(g: np.ndarray, m_max: float) -> np.ndarray:
    """Bounding-box half-widths of `{delta : delta^T g delta <= m_max}`.

    `sqrt(m_max * diag(inv(g)))`, directly comparable to a leaf's per-axis `dparam`
    half-width (C3). Note this is the *bounding box* of the ellipsoid, so it is the
    largest excursion along each axis, not a covering box for the ellipsoid's interior
    in any lattice sense.
    """
    if m_max < 0:
        msg = f"m_max must be >= 0, got {m_max}"
        raise ValueError(msg)
    return np.sqrt(m_max * np.diag(np.linalg.inv(g)))


def shift_matrix(delta_t: float, n_params: int) -> np.ndarray:
    """The Taylor re-centring matrix `T(delta_t)` restricted to the branchable axes.

    Built exactly as `transforms.shift_taylor_errors` builds it, then restricted to the
    leading `n_params x n_params` block. That restriction is exact rather than an
    approximation: the full matrix is lower triangular in the reverse-`k` ordering and
    `d_0` is the *last* row, so the branchable axes transform among themselves and never
    pick up a `d_0` contribution.

    Sign convention (C4): `delta_t = t_new - t_old`, and `d_new = T @ d_old`.
    """
    idx = np.arange(n_params)
    powers = np.tril(idx[:, np.newaxis] - idx)
    fact = np.vectorize(math.factorial)(powers).astype(np.float64)
    return np.asarray(
        delta_t**powers / fact * np.tril(np.ones_like(powers)), dtype=np.float64
    )


def transform_metric(g: np.ndarray, t_mat: np.ndarray) -> np.ndarray:
    """Metric in the shifted coordinates: `g' = T^-T g T^-1` for `d' = T d`.

    Invariance of `m` fixes this: with `delta' = T delta`,
    `delta^T g delta = delta'^T (T^-T g T^-1) delta'`.
    """
    t_inv = np.linalg.inv(t_mat)
    g_new = t_inv.T @ g @ t_inv
    return 0.5 * (g_new + g_new.T)


def cholesky_factor(g: np.ndarray) -> np.ndarray:
    """Lower-triangular `L` with `g = L @ L.T`.

    Whitening coordinates are `w = L.T @ delta`, in which `m = w @ w` and the
    mismatch ellipsoid is the unit ball — the property Phase 2's lattice relies on.
    """
    return np.linalg.cholesky(g)


class MMaxBridge(NamedTuple):
    """Ways of tying `m_max` to the existing `eta`-box (D4; resolves O2).

    The old criterion is sup-norm and the new one mean-square, so no single bridge is
    canonical. All of them are returned and the caller states which it used.

    Attributes
    ----------
    per_axis
        `m_max` matching the box exactly on each axis, in leaf axis order.
    axis_tight
        `min(per_axis)`: ellipsoid inscribed in the box, never wider on any axis.
    axis_loose
        `max(per_axis)`: box inscribed in the ellipsoid.
    volume
        `m_max` giving the ellipsoid the same volume as the box.
    box_half_widths
        The box half-widths used, for reporting.
    """

    per_axis: np.ndarray
    axis_tight: float
    axis_loose: float
    volume: float
    box_half_widths: np.ndarray


def m_max_from_eta(
    eta: float,
    nbins: int,
    poly_order: int,
    t_ref: float,
    t_start: float,
    t_end: float,
    f0: float,
    *,
    use_cheby: bool = False,
) -> MMaxBridge:
    """Bridge the `eta`-box criterion to an `m_max` — for comparison only.

    Not a design rule: it exists so Phase 3 can compare like with like. Per D5 the
    default is `use_cheby=False`, i.e. the **un-coarsened** box, because `"metric"` does
    not inherit the `2**k` coarsening. Pass `use_cheby=True` to bridge against the box
    the existing strategies actually use today.
    """
    from pyloki.utils import psr_utils  # local import: avoids an import cycle

    g = poly_phase_metric(t_ref, t_start, t_end, poly_order, f0, nbins)
    # poly_taylor_step_d_vec wants the span it is sizing for, and returns full step
    # sizes in leaf axis order; a leaf's column 1 is the HALF-width (C3).
    steps = psr_utils.poly_taylor_step_d_vec(
        poly_order,
        t_end - t_start,
        nbins,
        eta,
        np.array([f0], dtype=np.float64),
        t_ref=0,
        use_cheby=use_cheby,
    )[0]
    half = np.asarray(steps, dtype=np.float64) / 2.0

    g_inv_diag = np.diag(np.linalg.inv(g))
    per_axis = half**2 / g_inv_diag

    n = poly_order
    # Ellipsoid volume = V_n * sqrt(det(inv g)) * m_max**(n/2);
    # box volume = prod(2 * half).
    log_unit_ball = (n / 2.0) * math.log(math.pi) - math.lgamma(n / 2.0 + 1.0)
    sign, log_det_inv = np.linalg.slogdet(np.linalg.inv(g))
    if sign <= 0:
        msg = "metric is not positive definite; cannot match volume"
        raise ValueError(msg)
    log_box = float(np.sum(np.log(2.0 * half)))
    log_m = (2.0 / n) * (log_box - log_unit_ball - 0.5 * log_det_inv)

    return MMaxBridge(
        per_axis=per_axis,
        axis_tight=float(np.min(per_axis)),
        axis_loose=float(np.max(per_axis)),
        volume=float(math.exp(log_m)),
        box_half_widths=half,
    )
