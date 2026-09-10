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

from pyloki.utils import psr_utils
from pyloki.utils.misc import C_VAL

__all__ = [
    "MMaxBridge",
    "cholesky_factor",
    "ellipsoid_axis_extents",
    "harmonic_weight",
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
    """Return the derivative order carried by each axis, in leaf order (C1).

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


def harmonic_weight(
    nbins: int,
    ducy: float,
    *,
    weighting: str = "power",
    filter_width: int | None = None,
) -> float:
    """Harmonic weighting factor `<n**2>` to fold into `g` (D10; resolves O5(a)).

    A fundamental-only metric mis-predicts the real S/N loss: folding with wrong
    parameters multiplies harmonic `n` of the profile by
    `kappa_n = <exp(2*pi*i*n*dPhi)>`, so the penalty grows as `n**2`. To second order,

        A / A_0 ~ 1 - 2*pi**2 * Var(dPhi) * <n**2>,

    so the correction is a single scalar `<n**2>`, the weighted mean squared harmonic
    number. Two weightings, because it depends on the filter:

    - ``"power"`` (default): weight by `|p_n|**2`, the ideal matched filter.
      **Validated** against a direct `kappa_n`: `loss / (m * <n**2>)` is 0.96-1.00
      for `ducy` in 0.05-0.5 (see `tests/test_metric.py`).
    - ``"cross"``: weight by `Re(b_n^* p_n)`, the cross-spectrum with the boxcar filter
      the search actually scores with. Closer to the pipeline, and 2-3x *smaller*. It
      cannot be validated in isolation: a boxcar has sinc sidelobes, so `Re(b_n^* p_n)`
      is negative for some `n`, and attenuating those harmonics can *increase* the
      score: boxcar loss is not monotonic in smearing. Measured directly it even
      returns a small negative loss at `ducy=0.1`.

    ``"power"`` is the default for two reasons: it is the one that is validated, and
    it is the conservative choice for a *covering* criterion, since over-predicting the
    loss yields smaller leaves and so a safer covering. Phase 3's recalibration is the
    place to decide whether the 2-3x of ``"cross"`` is worth claiming.

    Parameters
    ----------
    nbins
        Number of phase bins in the folded profile.
    ducy
        Duty cycle; sets the profile width (and the boxcar width for ``"cross"``).
    weighting
        ``"power"`` or ``"cross"``, as above.
    filter_width
        Boxcar width in bins for ``"cross"``; defaults to `round(ducy * nbins)`.

    Returns
    -------
    float
        The weighting factor. Grows roughly as `ducy**-2`.
    """
    # Local import: pyloki.simulation.pulse imports pyloki.core, and core/__init__
    # reaches this module via taylor, so a top-level import here is a cycle
    # whenever pulse is the entry point.
    from pyloki.simulation.pulse import generate_folded_profile  # noqa: PLC0415

    if not 0.0 < ducy < 1.0:
        msg = f"ducy must be in (0, 1), got {ducy}"
        raise ValueError(msg)
    if weighting not in {"power", "cross"}:
        msg = f"weighting must be 'power' or 'cross', got {weighting!r}"
        raise ValueError(msg)

    profile = np.asarray(generate_folded_profile(nbins=nbins, ducy=ducy))
    p_spec = np.fft.rfft(profile)[1:]  # drop DC: a constant baseline carries no phase
    n = np.arange(1, len(p_spec) + 1)

    if weighting == "power":
        w = np.abs(p_spec) ** 2
    else:
        width = filter_width or max(1, round(ducy * nbins))
        box = np.zeros(nbins)
        box[:width] = 1.0 / math.sqrt(width)
        box = np.roll(box, int(np.argmax(profile)) - width // 2)
        w = np.real(np.conj(np.fft.rfft(box)[1:]) * p_spec)

    denom = float(np.sum(w))
    if denom <= 0:
        msg = f"degenerate weighting for nbins={nbins}, ducy={ducy}, {weighting!r}"
        raise ValueError(msg)
    return float(np.sum(w * n**2) / denom)


def poly_phase_metric(
    t_ref: float,
    t_start: float,
    t_end: float,
    poly_order: int,
    f0: float,
    nbins: int | None = None,
    ducy: float | None = None,
    weighting: str = "power",
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
        Phase bins. Needed only when `ducy` is given, to compute the harmonic weight.
    ducy
        Duty cycle. When given, `g` is scaled by `harmonic_weight(nbins, ducy)` so that
        `m` predicts the loss under the search's boxcar scoring (D10). When `None`, `g`
        is the single-harmonic (fundamental-only) metric, which **under-predicts the
        real loss** — pass `ducy` for anything that sizes a leaf.
    weighting
        Passed to `harmonic_weight`; see there. Default `"power"`.

    Returns
    -------
    np.ndarray
        `(poly_order, poly_order)` symmetric positive-definite matrix, ordered to match
        leaf rows `[:-2]`, such that `delta^T g delta` is the fractional amplitude loss.

    Notes
    -----
    `g_ij = W * 2*pi**2 * (f0/c)**2 * Cov(tau^ki / ki!, tau^kj / kj!)`, the covariance
    being a time average over the interval and `W` the harmonic weight (1 if `ducy` is
    None). Subtracting the product of means projects out the constant-phase mode, which
    is the unobservable direction.
    """
    if poly_order < 1:
        msg = f"poly_order must be >= 1, got {poly_order}"
        raise ValueError(msg)
    if ducy is not None and nbins is None:
        msg = "nbins is required when ducy is given"
        raise ValueError(msg)
    orders = _axis_orders(poly_order)
    inv_fact = np.array([1.0 / math.factorial(int(k)) for k in orders])

    # First moments <tau^k / k!> per axis, and cross moments <tau^(ki+kj)>.
    mean = np.array(
        [_tau_moment(int(k), t_ref, t_start, t_end) for k in orders],
    ) * inv_fact
    cross = np.empty((poly_order, poly_order), dtype=np.float64)
    for i in range(poly_order):
        for j in range(poly_order):
            n = int(orders[i] + orders[j])
            cross[i, j] = _tau_moment(n, t_ref, t_start, t_end)
    cross = cross * inv_fact[:, None] * inv_fact[None, :]

    cov = cross - np.outer(mean, mean)
    weight = (
        1.0 if ducy is None else harmonic_weight(int(nbins), ducy, weighting=weighting)
    )
    g = weight * _AMPLITUDE_SCALE * (f0 / C_VAL) ** 2 * cov
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
    """Build the Taylor re-centring matrix `T(delta_t)` for the branchable axes.

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
        delta_t**powers / fact * np.tril(np.ones_like(powers)), dtype=np.float64,
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
    ducy: float | None = None,
    use_cheby: bool = False,
) -> MMaxBridge:
    """Bridge the `eta`-box criterion to an `m_max` — for comparison only.

    Not a design rule: it exists so Phase 3 can compare like with like. Per D5 the
    default is `use_cheby=False`, i.e. the **un-coarsened** box, because `"metric"` does
    not inherit the `2**k` coarsening. Pass `use_cheby=True` to bridge against the box
    the existing strategies actually use today.
    """
    g = poly_phase_metric(t_ref, t_start, t_end, poly_order, f0, nbins, ducy)
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


# --------------------------------------------------------------------------------
# Phase 2: covering a parent region with metric-sized children.
# --------------------------------------------------------------------------------


def cubic_lattice_spacing(n_dim: int, m_max: float) -> float:
    """Spacing of a hypercubic lattice whose covering radius is `sqrt(m_max)`.

    In whitened coordinates a child's region is the ball of radius `r = sqrt(m_max)`.
    A cubic lattice of spacing `a` in `n` dimensions has covering radius
    `a * sqrt(n) / 2` (the worst case is a cell corner), so `a = 2 * r / sqrt(n)`
    makes the balls centred on lattice points cover the whole space.

    This is the deliberately unoptimised choice of `metric_PLAN.md`: A_n* would need
    fewer points for the same covering radius, but the plan defers that until Phase 3
    shows redundancy is the bottleneck.
    """
    if n_dim < 1:
        msg = f"n_dim must be >= 1, got {n_dim}"
        raise ValueError(msg)
    if m_max <= 0:
        msg = f"m_max must be > 0, got {m_max}"
        raise ValueError(msg)
    return 2.0 * math.sqrt(m_max) / math.sqrt(n_dim)


def _retention_form(a_mat: np.ndarray, m_max: float) -> np.ndarray:
    """Quadratic form `q` whose region `w^T q w <= 1` holds every needed lattice point.

    In whitened child coordinates the parent region is the ellipsoid
    `E = {w : w^T a_mat w <= m_max}` and each child covers a ball of radius
    `r = sqrt(m_max)`. A lattice point can only be needed if it lies within `r` of some
    point of `E`, i.e. inside the Minkowski sum `E + B(r)`, so any outer approximation
    of that sum gives a set whose retention keeps coverage exact.

    `E + B(r)` is not itself an ellipsoid. Two candidate outer forms:

    - Dilating uniformly by `1 + sqrt(lambda_max(a_mat))` -- equivalently
      `1 + r / a_min`, since `a_min = sqrt(m_max / lambda_max)`. Sound, and tight for
      the *shortest* semi-axis, but it scales every longer axis by the same factor. With
      observed axis ratios up to 335 (D8) that inflated the enumerated set 167x over the
      volume bound.
    - Dilating each semi-axis as `a_i -> a_i + r`. Tempting and much tighter, but
      UNSOUND: comparing support functions needs
      `|v| * sqrt(sum a_i^2 v_i^2) <= sum a_i v_i^2`, which is Cauchy-Schwarz the wrong
      way round, and sampling `E + B(r)` puts points at 1.40 in that form.

    Used here instead is the standard S-procedure external ellipsoid of a Minkowski sum:
    for any `t` in `(0, 1)`,

        E + B(r)  subset  {w : sum_i (w . q_i)^2 / c_i^2 <= 1},
        c_i^2 = a_i^2 / t + r^2 / (1 - t),

    with `q_i` the eigenvectors of `a_mat`. This is sound for every `t`, and for equal
    semi-axes it reduces to `c = a + r` exactly, so nothing is given away in the
    isotropic case. `t` is chosen to minimise the volume `prod c_i`.
    """
    lam, evec = np.linalg.eigh(a_mat)
    if np.any(lam <= 0.0):
        msg = "parent metric must be positive definite in whitened coordinates"
        raise ValueError(msg)
    a_sq = m_max / lam
    # Minimise sum(log c_i^2) over t; smooth and unimodal, so a fine grid then a local
    # refinement is ample and avoids a solver dependency.
    grid = np.linspace(1e-4, 1.0 - 1e-4, 2001)
    c_sq = a_sq[None, :] / grid[:, None] + m_max / (1.0 - grid)[:, None]
    best = int(np.argmin(np.log(c_sq).sum(axis=1)))
    lo = grid[max(best - 1, 0)]
    hi = grid[min(best + 1, grid.size - 1)]
    fine = np.linspace(lo, hi, 2001)
    c_sq = a_sq[None, :] / fine[:, None] + m_max / (1.0 - fine)[:, None]
    c_sq_best = c_sq[int(np.argmin(np.log(c_sq).sum(axis=1)))]
    return (evec * (1.0 / c_sq_best)) @ evec.T


def _enumerate_lattice_in_ellipsoid(
    q_mat: np.ndarray,
    bound: float,
    *,
    max_points: int | None = None,
) -> np.ndarray:
    """Integer vectors `z` with `z^T q_mat z <= bound`, by Fincke-Pohst enumeration.

    Writing `q_mat = U^T U` with `U` upper triangular,
    `z^T q_mat z = sum_i (sum_{j>=i} U[i,j] z_j)**2`, so the coordinates can be fixed
    from the last to the first with an exact interval at each level and the partial sum
    subtracted from the budget. Only points that can still satisfy the bound are
    visited, which is what makes an elongated ellipsoid tractable -- enumerating its
    axis-aligned bounding box instead is hopeless once the metric is ill-conditioned.

    Raises
    ------
    ValueError
        If more than `max_points` points satisfy the bound; the caller should widen
        `m_max` or accept a coarser covering rather than silently truncate.
    """
    n_dim = q_mat.shape[0]
    upper = np.linalg.cholesky(q_mat).T  # q_mat = upper^T @ upper
    out: list[np.ndarray] = []
    z = np.zeros(n_dim, dtype=np.int64)

    def recurse(level: int, budget: float) -> None:
        if budget < -1e-12:
            return
        if level < 0:
            out.append(z.copy())
            if max_points is not None and len(out) > max_points:
                msg = (
                    f"metric branching kept more than {max_points} lattice points; "
                    f"the parent region is too large for the child spacing"
                )
                raise ValueError(msg)
            return
        # Partial sum from the already-fixed higher coordinates.
        tail = float(np.dot(upper[level, level + 1 :], z[level + 1 :]))
        diag = float(upper[level, level])
        span = math.sqrt(max(budget, 0.0))
        lo = math.ceil((-span - tail) / diag - 1e-12)
        hi = math.floor((span - tail) / diag + 1e-12)
        for value in range(lo, hi + 1):
            z[level] = value
            term = diag * value + tail
            recurse(level - 1, budget - term * term)
        z[level] = 0

    recurse(n_dim - 1, bound)
    if not out:
        return np.zeros((1, n_dim), dtype=np.int64)
    return np.asarray(out, dtype=np.int64)


def lattice_children(
    g_parent: np.ndarray,
    g_child: np.ndarray,
    m_max: float,
    *,
    max_children: int | None = None,
) -> np.ndarray:
    """Offsets of child centres covering a parent region, in Taylor coordinates.

    Parameters
    ----------
    g_parent
        Metric defining the parent's region, `{d : d^T g_parent d <= m_max}`, i.e. the
        metric of the previous stage's interval.
    g_child
        Metric defining each child's region, i.e. the current stage's. Children are
        spaced so their regions cover the parent's.
    m_max
        Mismatch budget, shared by parent and children.
    max_children
        Optional cap; exceeding it raises rather than silently under-covering.

    Returns
    -------
    np.ndarray
        `(n_children, n_params)` offsets relative to the parent centre, in the same
        Taylor coordinates and axis order as the leaf array (C1). Always contains at
        least one point.

    Notes
    -----
    Coverage is guaranteed by construction, not by sampling: the lattice spacing gives
    covering radius `sqrt(m_max)` everywhere, and `_retention_form` keeps every
    lattice point within that radius of the parent region.
    """
    n_dim = g_child.shape[0]
    chol = cholesky_factor(g_child)
    # Whitened coordinates w = L^T d: a child's region is the ball of radius
    # sqrt(m_max), the parent's is the ellipsoid w^T A w <= m_max.
    a_mat = np.linalg.solve(chol, np.linalg.solve(chol, g_parent).T).T
    a_mat = 0.5 * (a_mat + a_mat.T)

    spacing = cubic_lattice_spacing(n_dim, m_max)
    form = _retention_form(a_mat, m_max)

    # Lattice points are w = spacing * z, so w^T form w <= 1 becomes
    # z^T (spacing**2 * form) z <= 1.
    z_int = _enumerate_lattice_in_ellipsoid(
        form * spacing**2, 1.0, max_points=max_children,
    )
    kept = z_int.astype(np.float64) * spacing
    # Back to Taylor coordinates: d = L^-T w.
    return np.linalg.solve(chol.T, kept.T).T
