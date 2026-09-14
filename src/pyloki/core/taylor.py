from __future__ import annotations

import math

import numpy as np
from numba import njit

from pyloki.core import metric
from pyloki.core.common import get_leaves
from pyloki.utils import np_utils, psr_utils, transforms
from pyloki.utils.misc import C_VAL, FLOAT_EPSILON
from pyloki.utils.snail import MiddleOutScheme


@njit(cache=True, fastmath=True)
def poly_taylor_seed(
    param_arr: list[np.ndarray],
    dparams: np.ndarray,
    poly_order: int,
    coord_init: tuple[float, float],
) -> np.ndarray:
    """Generate the seed leaves for Taylor polynomial search.

    Parameters
    ----------
    param_arr : list[np.ndarray]
        Parameter array for each dimension; only (acceleration, frequency).
    dparams : np.ndarray
        Parameter step (grid) sizes for each dimension. Shape is (poly_order,).
        Order is reversed [..., acc, freq].
    poly_order : int
        The order of the Taylor polynomial.
    coord_init : tuple[float, float]
        The coordinate of the starting segment (level 0).
        - coord_init[0] -> t0 (reference time) measured from t=0
        - coord_init[1] -> scale (half duration of the segment)

    Returns
    -------
    np.ndarray
        The seed leaves. Shape is (n_leaves, poly_order + 2, 2).

    Notes
    -----
    Conventions for each seed leaf:
    leaf[:-1, 0] -> Taylor polynomial coefficients,
                    order is [d_poly_order, ..., d_1, d_0]
    leaf[:-1, 1] -> Grid size (error) on each coefficient,
    leaf[-1, 0]  -> Frequency at t_init (f0), assuming f=f0 at t_init
    leaf[-1, 1]  -> Flag to indicate basis change (0: Polynomial, 1: Physical)
    """
    _, _ = coord_init
    leaves_taylor = get_leaves(param_arr, dparams)
    f0_batch = leaves_taylor[:, -1, 0]
    df_batch = leaves_taylor[:, -1, 1]
    leaves = np.zeros((len(leaves_taylor), poly_order + 2, 2), dtype=np.float64)
    # Copy till accel
    leaves[:, :-3] = leaves_taylor[:, :-1]
    # f = f0(1 - v / C) => dv = -(C/f0) * df
    leaves[:, -3, 0] = 0
    leaves[:, -3, 1] = df_batch * (C_VAL / f0_batch)
    # intialize d0 (measure from t=t_init)
    leaves[:, -2, 0] = 0  # we never branch on d0
    leaves[:, -1, 0] = f0_batch
    leaves[:, -1, 1] = 0  # Polynomial basis
    return leaves


@njit(cache=True, fastmath=True)
def poly_taylor_branch_batch(
    leaves_batch: np.ndarray,
    coord_cur: tuple[float, float],
    nbins: int,
    eta: float,
    poly_order: int,
    branch_max: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Branch a batch of tree parameter nodes to leaves.

    Parameters
    ----------
    leaves_batch : np.ndarray
        Leaf parameter sets. Shape: (n_leaves, poly_order + 2, 2).
    coord_cur : tuple[float, float]
        Coordinates for the accumulated segment in the current stage.
    nbins : int
        Number of bins in the folded profile.
    eta : float
        Tolerance for the parameter step size in bins.
    poly_order : int
        The order of the Taylor polynomial.
    param_limits : np.ndarray
        The limits for each parameter in Taylor basis (reverse order).
    branch_max : int
        Maximum number of branches that can be generated.

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        - leaves_branch_batch: Array of leaf centers.
          Shape: (n_branch, poly_order + 2, 2).
        - batch_origins: Array of original indices.
          Shape: (n_branch,).
    """
    n_batch, _, _ = leaves_batch.shape
    _, t_obs_minus_t_ref = coord_cur
    n_params = poly_order

    param_cur_batch = leaves_batch[:, :-2, 0]
    dparam_cur_batch = leaves_batch[:, :-2, 1]
    d0_cur_batch = leaves_batch[:, -2, 0]
    f0_batch = leaves_batch[:, -1, 0]
    basis_flag_batch = leaves_batch[:, -1, 1]

    dparam_new_batch = psr_utils.poly_taylor_step_d_vec(
        n_params,
        t_obs_minus_t_ref,
        nbins,
        eta,
        f0_batch,
        t_ref=0,
    )
    shift_bins_batch = psr_utils.poly_taylor_shift_d_vec(
        dparam_cur_batch,
        dparam_new_batch,
        t_obs_minus_t_ref,
        nbins,
        f0_batch,
        t_ref=0,
    )

    # Vectorized Padded Branching
    pad_branched_params = np.empty((n_batch, n_params, branch_max), dtype=np.float64)
    branched_dparams = np.empty((n_batch, n_params), dtype=np.float64)
    branched_counts = np.empty((n_batch, n_params), dtype=np.int64)
    for i in range(n_batch):
        for j in range(n_params):
            dparam_act, count = psr_utils.branch_param_padded(
                pad_branched_params[i, j],
                param_cur_batch[i, j],
                dparam_cur_batch[i, j],
                dparam_new_batch[i, j],
            )
            branched_dparams[i, j] = dparam_act
            branched_counts[i, j] = count

    # Vectorized Selection
    for i in range(n_batch):
        for j in range(n_params):
            if shift_bins_batch[i, j] < (eta - FLOAT_EPSILON):
                pad_branched_params[i, j, :] = 0
                pad_branched_params[i, j, 0] = param_cur_batch[i, j]
                branched_dparams[i, j] = dparam_cur_batch[i, j]
                branched_counts[i, j] = 1

    # Optimized Padded Cartesian Product
    leaf_params_branch_cart, batch_origins = np_utils.cartesian_prod_padded(
        pad_branched_params,
        branched_counts,
        n_batch,
        n_params,
    )
    n_branch = len(batch_origins)
    leaves_branch_batch = np.zeros((n_branch, n_params + 2, 2), dtype=np.float64)
    leaves_branch_batch[:, :-2, 0] = leaf_params_branch_cart
    leaves_branch_batch[:, :-2, 1] = branched_dparams[batch_origins]
    leaves_branch_batch[:, -2, 0] = d0_cur_batch[batch_origins]
    leaves_branch_batch[:, -1, 0] = f0_batch[batch_origins]
    leaves_branch_batch[:, -1, 1] = basis_flag_batch[batch_origins]
    return leaves_branch_batch, batch_origins


def poly_taylor_branch_metric_batch(
    leaves_batch: np.ndarray,
    coord_cur: tuple[float, float],
    coord_prev: tuple[float, float],
    coord_next: tuple[float, float],
    nbins: int,
    ducy: float,
    poly_order: int,
    m_max: float,
    branch_max: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Branch a batch of leaves onto a metric-sized covering of each parent region.

    Same contract as `poly_taylor_branch_batch` -- identical leaf layout in and out,
    and the same `(leaves_branch_batch, batch_origins)` return -- but the children are
    a lattice covering of the parent's mismatch ellipsoid rather than an axis-aligned
    Cartesian product of per-axis steps.

    Parameters
    ----------
    leaves_batch : np.ndarray
        Leaf parameter sets. Shape: (n_leaves, poly_order + 2, 2).
    coord_cur : tuple[float, float]
        Coordinates for the accumulated segment in the current (child) stage.
    coord_prev : tuple[float, float]
        Coordinates for the previous (parent) stage. The box strategy does not need
        this -- a leaf's `dparam` column already encodes its own spacing -- but a metric
        does: per-axis half-widths cannot represent the parent ellipsoid's orientation,
        so the parent's own metric has to be rebuilt from its interval.
    coord_next : tuple[float, float]
        Coordinates the leaves will be re-centred on at the end of the stage. Only its
        reference time is used, to place the child's averaging window relative to the
        leaf's current epoch (D20).
    nbins, ducy : int, float
        Folded-profile resolution and duty cycle, setting the harmonic weighting.
    poly_order : int
        The order of the Taylor polynomial.
    m_max : float
        Mismatch budget, shared by parent and children.
    branch_max : int
        Cap on the number of children *per parent*. Note this differs from the box
        strategy's use of `branch_max` as a per-axis padding width; a metric covering
        is not separable per axis, so a total is the only meaningful cap.

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        - leaves_branch_batch: Array of leaf centers.
          Shape: (n_branch, poly_order + 2, 2).
        - batch_origins: Array of original indices. Shape: (n_branch,).

    Notes
    -----
    Not `@njit`: the covering needs `eigh`, Cholesky solves and a recursive
    Fincke-Pohst enumeration. `dyn_poly_taylor.branch_func` *is* `@njit` and so cannot
    dispatch here yet -- see `DECISIONS.md` D13.

    The metric is exactly proportional to `f0**2` (verified to machine precision), so
    the offsets scale as `1 / f0` and the enumeration -- much the most expensive part --
    runs once per stage rather than once per distinct `f0` in the batch.
    """
    ref_cur, t_half_cur = coord_cur
    _, t_half_prev = coord_prev
    ref_next, _ = coord_next
    offsets_unit, extents_unit = metric_branch_tables(
        t_half_prev,
        t_half_cur,
        ref_next - ref_cur,
        nbins,
        ducy,
        poly_order,
        m_max,
        branch_max,
    )
    return poly_taylor_branch_metric_apply(leaves_batch, offsets_unit, extents_unit)


def metric_branch_tables(
    t_half_prev: float,
    t_half_cur: float,
    delta_t: float,
    nbins: int,
    ducy: float,
    poly_order: int,
    m_max: float,
    branch_max: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Per-stage covering tables: child offsets at `f0 = 1`, and child axis extents.

    Everything expensive lives here -- `eigh`, Cholesky solves and the Fincke-Pohst
    enumeration -- and none of it depends on the leaves. The tables are a function of
    the *stage* alone, because the metric is exactly proportional to `f0**2` (D11), so
    a leaf's own `f0` enters only as a `1 / f0` rescaling. That is what lets the
    numba-hostile part run once per stage in Python while the per-batch work stays in
    `poly_taylor_branch_metric_apply`, which is `@njit`.

    Parameters
    ----------
    t_half_prev, t_half_cur
        **Half-widths** of the parent and child accumulated windows, i.e. `coord[1]`.
    delta_t
        `coord_next[0] - coord_cur[0]`: how far ahead of the leaf's expansion epoch the
        child window is centred.

    Notes
    -----
    The averaging windows, not the half-widths, are what `poly_phase_metric` needs, and
    getting that wrong was a real bug (D20). Both leaves are expanded about the
    *previous* centre. The parent's window is symmetric about it, `[-t_half_prev,
    +t_half_prev]`; the child's window has grown and moved on, so it is centred
    `delta_t` ahead: `[delta_t - t_half_cur, delta_t + t_half_cur]`. Passing
    `[0, coord[1]]` -- the epoch at the window's edge, over half its true length --
    misprices the extents by 4x to 64x.
    """
    g_parent = metric.poly_phase_metric(
        0.0, -t_half_prev, t_half_prev, poly_order, 1.0, nbins, ducy,
    )
    g_child = metric.poly_phase_metric(
        0.0, delta_t - t_half_cur, delta_t + t_half_cur, poly_order, 1.0, nbins, ducy,
    )
    offsets_unit = metric.lattice_children(
        g_parent, g_child, m_max, max_children=branch_max,
    )
    extents_unit = metric.ellipsoid_axis_extents(g_child, m_max)
    return np.ascontiguousarray(offsets_unit), np.ascontiguousarray(extents_unit)


def metric_transform_extents(
    t_half_cur: float,
    nbins: int,
    ducy: float,
    poly_order: int,
    m_max: float,
) -> np.ndarray:
    """Child axis extents once re-centred on the new epoch, at `f0 = 1`.

    A child's region is the ellipsoid `{d : d^T g d <= m_max}`. The transform moves the
    expansion epoch to the centre of the newly accumulated window, where that window is
    symmetric, so the honest per-axis half-width there is just
    `ellipsoid_axis_extents(g, m_max)` for `g` built about the new epoch.

    This is the *exact* transform the plan's Phase 2 step 2 asks for, and it cannot be
    done from leaf state alone: a leaf stores only the ellipsoid's bounding box (D12),
    and a bounding box does not determine the ellipsoid it bounds. It can be done here
    because the window is a property of the *stage*, so -- exactly as for the covering
    tables (D11) -- one evaluation at `f0 = 1` serves every leaf: `g` is proportional
    to `f0**2`, hence `inv(g)` to `f0**-2` and the extents to `1 / f0`.

    Equivalently one could carry the pre-transform metric through
    `metric.transform_metric(g, shift_matrix(delta_t, poly_order))`; Phase 1 test 2 is
    the statement that the two agree, and `TestTransformExtents` checks it here too.
    Rebuilding about the new epoch is cheaper and needs no `delta_t`.

    Notes
    -----
    Returns `poly_order` entries, for the branchable axes in leaf order (C1). `d_0` is
    the constant-phase mode, which `poly_phase_metric` projects out, so it has no
    defined half-width and keeps the zero that the metric branch gives it.
    """
    g_next = metric.poly_phase_metric(
        0.0, -t_half_cur, t_half_cur, poly_order, 1.0, nbins, ducy,
    )
    return np.ascontiguousarray(metric.ellipsoid_axis_extents(g_next, m_max))


def generate_bp_poly_taylor_metric(
    tseg_ffa: float,
    nsegments: int,
    ref_seg: int,
    nbins: int,
    ducy: float,
    poly_order: int,
    m_max: float,
    branch_max: int,
    *,
    use_moving_grid: bool,
) -> np.ndarray:
    """Exact per-level branching factor `B(s)` under `tiling_strategy="metric"`.

    `generate_bp_poly_taylor` cannot be reused: it models branching as a product of
    independent per-axis counts, which is precisely the axis-aligned assumption the
    metric covering replaces. It is also `@njit`, so it could not call the enumeration.

    Under `"metric"` the answer is both simpler and exact rather than averaged: the
    covering depends only on the stage (D11), so every parent at level `s` emits the
    same number of children, with no `f0` dependence to average over.
    """
    scheme = MiddleOutScheme(nsegments, ref_seg, tseg_ffa, stride=1)
    branching_pattern = np.empty(nsegments - 1, dtype=np.float64)
    for prune_level in range(1, nsegments):
        ref_cur, t_half_cur = scheme.get_current_coord(prune_level, use_moving_grid)
        _, t_half_prev = scheme.get_previous_coord(prune_level, use_moving_grid)
        ref_next, _ = scheme.get_coord(prune_level)
        offsets, _ = metric_branch_tables(
            t_half_prev,
            t_half_cur,
            ref_next - ref_cur,
            nbins,
            ducy,
            poly_order,
            m_max,
            branch_max,
        )
        branching_pattern[prune_level - 1] = float(len(offsets))
    return branching_pattern


@njit(cache=True, fastmath=True)
def poly_taylor_branch_metric_apply(
    leaves_batch: np.ndarray,
    offsets_unit: np.ndarray,
    extents_unit: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Apply per-stage covering tables to a batch of leaves.

    The numba-compatible half of the metric branch: a broadcast add and a `1 / f0`
    rescale, with no linear algebra and no enumeration. `offsets_unit` and
    `extents_unit` come from `metric_branch_tables`, which runs once per stage.

    Splitting the branch here is what unblocks dispatch from the `@njit`
    `dyn_poly_taylor.branch_func`: see DECISIONS.md D14.
    """
    n_batch = leaves_batch.shape[0]
    n_child = offsets_unit.shape[0]
    n_params = offsets_unit.shape[1]

    n_branch = n_batch * n_child
    leaves_branch_batch = np.zeros((n_branch, n_params + 2, 2), dtype=np.float64)
    batch_origins = np.empty(n_branch, dtype=np.int64)

    for i in range(n_batch):
        f0 = leaves_batch[i, -1, 0]
        inv_f0 = 1.0 / f0
        for c in range(n_child):
            row = i * n_child + c
            batch_origins[row] = i
            for j in range(n_params):
                leaves_branch_batch[row, j, 0] = (
                    leaves_batch[i, j, 0] + offsets_unit[c, j] * inv_f0
                )
                # The child's per-axis half-widths: the bounding box of its mismatch
                # ellipsoid, so consumers of the dparam column keep a width (C3).
                leaves_branch_batch[row, j, 1] = extents_unit[j] * inv_f0
            leaves_branch_batch[row, -2, 0] = leaves_batch[i, -2, 0]
            leaves_branch_batch[row, -1, 0] = f0
            leaves_branch_batch[row, -1, 1] = leaves_batch[i, -1, 1]
    return leaves_branch_batch, batch_origins


def metric_resolve_mismatch(
    leaves_batch: np.ndarray,
    coord_add: tuple[float, float],
    coord_cur: tuple[float, float],
    coord_init: tuple[float, float],
    param_arr: list[np.ndarray],
    param_grid_count_init: np.ndarray,
    param_limits: np.ndarray,
    nbins: int,
    ducy: float,
) -> np.ndarray:
    """Mismatch each child loses by being resolved onto the base grid `G0`.

    metric_PLAN.md Phase 2 step 4. Under the box strategies a child sits on a
    rectangular refinement of `G0` by construction, so resolving it is exact up to the
    cell it was built from. Under `"metric"` children sit wherever the lattice puts
    them, so `resolve` rounds them to the nearest `G0` cell centre and that rounding is
    an *extra* mismatch the covering never accounted for.

    This measures it. For each child it replays `poly_taylor_resolve_batch`'s forward
    map to the added segment, looks up the cell centre actually loaded, and returns the
    mismatch between the two in the **base-segment** metric -- the metric of the single
    FFA segment being added, which is the interval over which that fold is valid.

    Returns
    -------
    np.ndarray
        Per-child mismatch, directly comparable to `m_max`. Diagnostic only: nothing in
        the search consumes it, and this is not called on the hot path.

    Notes
    -----
    `G0` is two-dimensional -- `get_nearest_indices_2d_batch` grids only acceleration
    and frequency -- so the residual lives entirely in `[d_2, d_1]` and the base metric
    is built at `poly_order=2`. Higher derivatives are not resolved to a grid at all;
    they enter only through the shift to the segment's epoch.
    """
    t0_cur, _ = coord_cur
    t0_init, _ = coord_init
    t0_add, t_half_add = coord_add

    param_vec_batch = np.ascontiguousarray(leaves_batch[:, :-1, 0])
    f0_batch = leaves_batch[:, -1, 0]

    # Exactly poly_taylor_resolve_batch's forward map.
    dvec_t_add = transforms.shift_taylor_params(param_vec_batch, t0_add - t0_cur)
    dvec_t_init = transforms.shift_taylor_params(param_vec_batch, t0_init - t0_cur)
    accel_new = dvec_t_add[:, -3]
    vel_new = dvec_t_add[:, -2] - dvec_t_init[:, -2]
    freq_new = f0_batch * (1 - vel_new / C_VAL)
    param_idx = psr_utils.get_nearest_indices_2d_batch(
        accel_new,
        freq_new,
        param_grid_count_init,
        param_limits,
    )
    # What the loaded fold actually corresponds to: `range_param` lays down cell
    # centres and the index is the containing cell, so this is the nearest grid point.
    accel_grid = np.asarray(param_arr[-2])[param_idx[:, -2]]
    freq_grid = np.asarray(param_arr[-1])[param_idx[:, -1]]
    vel_grid = C_VAL * (1.0 - freq_grid / f0_batch)

    delta = np.empty((len(leaves_batch), 2), dtype=np.float64)
    delta[:, 0] = accel_grid - accel_new
    delta[:, 1] = vel_grid - vel_new

    # g is exactly proportional to f0**2 (D11), so build it once at f0 = 1 and scale.
    # The window is the added segment, symmetric about its own epoch (D20).
    g_unit = metric.poly_phase_metric(
        0.0, -t_half_add, t_half_add, 2, 1.0, nbins, ducy,
    )
    return f0_batch**2 * np.einsum("ni,ij,nj->n", delta, g_unit, delta)


@njit(cache=True, fastmath=True)
def poly_taylor_resolve_batch(
    leaves_batch: np.ndarray,
    coord_add: tuple[float, float],
    coord_cur: tuple[float, float],
    coord_init: tuple[float, float],
    param_grid_count_init: np.ndarray,
    param_limits: np.ndarray,
    nbins: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Resolve a batch of leaf params to find the closest grid index and phase shift.

    Parameters
    ----------
    leaves_batch : np.ndarray
        The leaf parameter set. Shape is (n_leaves, poly_order + 2, 2).
    coord_add : tuple[float, float]
        The coordinates for the added segment (level current).
    coord_cur : tuple[float, float]
        The coordinates for the current pruning suggestion tree.
    coord_init : tuple[float, float]
        The coordinates for the initial pruning suggestion tree.
    param_grid_count_init : np.ndarray
        Number of points in the initial (FFA) grid for the ``coord_add`` segment
        Currently this is simply [n_accel, n_freq].
    param_limits : np.ndarray
        Parameter limits (min, max).
    nbins : int
        Number of bins in the folded profile.

    Returns
    -------
    tuple[np.ndarray, float]
        The resolved parameter index in the ``param_arr`` and the relative phase shift.

    Notes
    -----
    leaf is referenced to coord_cur, so we need to shift it to coord_add to get
    the resolved parameters index and relative phase shift. We also need to correct for
    the tree phase offset from coord_init to coord_cur.

    relative_phase is complete phase shift with fractional part.
    """
    t0_cur, _ = coord_cur
    t0_init, _ = coord_init
    t0_add, _ = coord_add

    param_vec_batch = leaves_batch[:, :-1, 0]
    f0_batch = leaves_batch[:, -1, 0]

    dvec_t_add = transforms.shift_taylor_params(param_vec_batch, t0_add - t0_cur)
    dvec_t_init = transforms.shift_taylor_params(param_vec_batch, t0_init - t0_cur)
    accel_new_batch = dvec_t_add[:, -3]
    vel_new_batch = dvec_t_add[:, -2] - dvec_t_init[:, -2]
    freq_new_batch = f0_batch * (1 - vel_new_batch / C_VAL)
    delay_batch = (dvec_t_add[:, -1] - dvec_t_init[:, -1]) / C_VAL
    relative_phase_batch = psr_utils.get_phase_idx(
        t0_add - t0_init,
        f0_batch,
        nbins,
        delay_batch,
    )
    # Pass the full param_limits to infer correct n_params
    param_idx_batch = psr_utils.get_nearest_indices_2d_batch(
        accel_new_batch,
        freq_new_batch,
        param_grid_count_init,
        param_limits,
    )
    return param_idx_batch, relative_phase_batch


@njit(cache=True, fastmath=True)
def poly_taylor_fixed_resolve_batch(
    leaves_batch: np.ndarray,
    coord_add: tuple[float, float],
    coord_init: tuple[float, float],
    param_grid_count_init: np.ndarray,
    param_limits: np.ndarray,
    nbins: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Resolve a batch of leaf params to find the closest grid index and phase shift."""
    t0_init, _ = coord_init
    t0_add, _ = coord_add

    param_vec_batch = leaves_batch[:, :-1, 0]
    f0_batch = leaves_batch[:, -1, 0]

    dvec_t_add = transforms.shift_taylor_params(param_vec_batch, t0_add - t0_init)
    accel_new_batch = dvec_t_add[:, -3]
    vel_new_batch = dvec_t_add[:, -2]
    freq_new_batch = f0_batch * (1 - vel_new_batch / C_VAL)
    delay_batch = dvec_t_add[:, -1] / C_VAL
    relative_phase_batch = psr_utils.get_phase_idx(
        t0_add - t0_init,
        f0_batch,
        nbins,
        delay_batch,
    )
    param_idx_batch = psr_utils.get_nearest_indices_2d_batch(
        accel_new_batch,
        freq_new_batch,
        param_grid_count_init,
        param_limits,
    )
    return param_idx_batch, relative_phase_batch


@njit(cache=True, fastmath=True)
def poly_taylor_ascend_resolve_batch(
    leaves_batch: np.ndarray,
    coord_segments: np.ndarray,
    coord_cur: tuple[float, float],
    param_grid_count_init: np.ndarray,
    param_limits: np.ndarray,
    nbins: int,
) -> tuple[np.ndarray, np.ndarray]:
    t0_cur, _ = coord_cur
    n_leaves = len(leaves_batch)
    nsegments = len(coord_segments)
    n_params = param_limits.shape[0]

    param_vec_batch = leaves_batch[:, :-1, 0]
    f0_batch = leaves_batch[:, -1, 0]

    param_idx_batch_arr = np.empty((n_leaves, nsegments, n_params), dtype=np.int64)
    relative_phase_batch_arr = np.empty((n_leaves, nsegments), dtype=np.float64)
    for isegment in range(nsegments):
        t0_seg, _ = coord_segments[isegment]
        dvec_t_seg = transforms.shift_taylor_params(param_vec_batch, t0_seg - t0_cur)
        accel_new_batch = dvec_t_seg[:, -3]
        freq_new_batch = f0_batch * (1 - dvec_t_seg[:, -2] / C_VAL)
        delay_batch = dvec_t_seg[:, -1] / C_VAL
        relative_phase_batch = psr_utils.get_phase_idx(
            t0_seg - t0_cur,
            f0_batch,
            nbins,
            delay_batch,
        )
        # Pass the full param_limits to infer correct n_params
        param_idx_batch = psr_utils.get_nearest_indices_2d_batch(
            accel_new_batch,
            freq_new_batch,
            param_grid_count_init,
            param_limits,
        )
        param_idx_batch_arr[:, isegment, :] = param_idx_batch
        relative_phase_batch_arr[:, isegment] = relative_phase_batch
    return param_idx_batch_arr, relative_phase_batch_arr


@njit(cache=True, fastmath=True)
def poly_taylor_transform_batch(
    leaves_batch: np.ndarray,
    coord_next: tuple[float, float],
    coord_cur: tuple[float, float],
    tiling_strategy: str,
    transform_extents: np.ndarray,
) -> np.ndarray:
    """Re-center (in-place) the leaves to the next segment reference time.

    Under `tiling_strategy="metric"` the values shift exactly as they always have, but
    column 1 does not come from the old column 1: it is the per-stage table built by
    `metric_transform_extents` and rescaled by `1 / f0`. The box strategies propagate
    their half-widths through `shift_taylor_full` and ignore `transform_extents`.
    """
    delta_t = coord_next[0] - coord_cur[0]
    leaves_batch_trans = np.zeros_like(leaves_batch)
    if tiling_strategy == "metric":
        if transform_extents.shape[0] == 0:
            msg = (
                "tiling_strategy='metric' requires per-stage transform extents; "
                "call taylor.metric_transform_extents and pass them to transform()"
            )
            raise ValueError(msg)
        leaves_batch_trans[:, :-1, 0] = transforms.shift_taylor_params(
            np.ascontiguousarray(leaves_batch[:, :-1, 0]),
            delta_t,
        )
        n_params = transform_extents.shape[0]
        for i in range(leaves_batch.shape[0]):
            inv_f0 = 1.0 / leaves_batch[i, -1, 0]
            for j in range(n_params):
                leaves_batch_trans[i, j, 1] = transform_extents[j] * inv_f0
        # Row [-2] is d_0, the constant-phase mode the metric projects out. It has no
        # defined half-width, and keeps the zero the metric branch gave it.
    else:
        leaves_batch_trans[:, :-1] = transforms.shift_taylor_full(
            leaves_batch[:, :-1],
            delta_t,
            tiling_strategy,
        )
    leaves_batch_trans[:, -1] = leaves_batch[:, -1]
    return leaves_batch_trans


@njit(cache=True, fastmath=True)
def poly_taylor_report_batch(leaves_batch: np.ndarray) -> np.ndarray:
    param_sets_batch = leaves_batch.copy()
    param_sets_vals = leaves_batch[:, :-3, 0]
    param_sets_sigs = leaves_batch[:, :-3, 1]
    v_final = leaves_batch[:, -3, 0]
    dv_final = leaves_batch[:, -3, 1]
    f0_batch = leaves_batch[:, -1, 0]
    s_factor = 1 - v_final / C_VAL
    # Gauge transform + error propagation
    param_sets_batch[:, :-3, 0] = param_sets_vals / s_factor[:, None]
    param_sets_batch[:, :-3, 1] = np.sqrt(
        (param_sets_sigs / s_factor[:, None]) ** 2
        + ((param_sets_vals / (C_VAL * s_factor[:, None] ** 2)) ** 2)
        * (dv_final[:, None] ** 2),
    )
    param_sets_batch[:, -3, 0] = f0_batch * s_factor
    param_sets_batch[:, -3, 1] = f0_batch * dv_final / C_VAL
    return param_sets_batch


@njit(cache=True, fastmath=True)
def generate_bp_poly_taylor_approx(
    param_arr: list[np.ndarray],
    dparams_act: np.ndarray,
    tseg_ffa: float,
    nsegments: int,
    nbins: int,
    eta: float,
    ref_seg: int,
    use_moving_grid: bool,
    tiling_strategy: str,
    itree: int = 0,
    branch_max: int = 256,
) -> np.ndarray:
    """Generate the approximate branching pattern for the Taylor pruning search."""
    poly_order = len(dparams_act)
    snail_scheme = MiddleOutScheme(nsegments, ref_seg, tseg_ffa, stride=1)
    coord_init = snail_scheme.get_coord(0)
    leaves_init = poly_taylor_seed(param_arr, dparams_act, poly_order, coord_init)
    leaf = leaves_init[itree : itree + 1]  # shape: (1, total_size)
    branching_pattern = np.empty(nsegments - 1, dtype=np.float64)
    for prune_level in range(1, nsegments):
        coord_next = snail_scheme.get_coord(prune_level)
        coord_cur = snail_scheme.get_current_coord(
            prune_level,
            moving_grid=use_moving_grid,
        )
        leaves_arr, _ = poly_taylor_branch_batch(
            leaf,
            coord_cur,
            nbins,
            eta,
            poly_order,
            branch_max,
        )
        branching_pattern[prune_level - 1] = len(leaves_arr)
        if use_moving_grid:
            leaves_arr = poly_taylor_transform_batch(
                leaves_arr,
                coord_next,
                coord_cur,
                tiling_strategy,
                # No metric table: "metric" is refused before reaching the approximate
                # branching pattern, which models branching per axis anyway.
                np.empty(0, dtype=np.float64),
            )
        leaf = leaves_arr[0:1]  # shape: (1, total_size)
    # Check if any branches is truncated due to branch_max
    if np.any(branching_pattern == branch_max):
        msg = "Branching pattern is truncated due to branch_max. Increase branch_max."
        raise ValueError(msg)
    return branching_pattern


@njit(cache=True, fastmath=True)
def generate_bp_poly_taylor(
    param_arr: list[np.ndarray],
    dparams_act: np.ndarray,
    tseg_ffa: float,
    nsegments: int,
    nbins: int,
    eta: float,
    ref_seg: int,
    use_moving_grid: bool,
    tiling_strategy: str,
    use_cheby_coarsening: bool = True,
) -> np.ndarray:
    """Generate the exact branching pattern for the Taylor pruning search."""
    n_params = len(dparams_act)
    f0_batch = param_arr[-1]
    n_freqs = len(f0_batch)
    snail_scheme = MiddleOutScheme(nsegments, ref_seg, tseg_ffa, stride=1)
    weights = np.ones(n_freqs, dtype=np.float64)
    branching_pattern = np.empty(nsegments - 1, dtype=np.float64)

    dparam_cur_batch = np.empty((n_freqs, n_params), dtype=np.float64)
    dparam_cur_next = np.empty((n_freqs, n_params), dtype=np.float64)
    dparam_d_vec = np.empty((n_freqs, n_params + 1), dtype=np.float64)
    for i in range(n_freqs):
        dparam_cur_batch[i, :n_params] = dparams_act
    # f = f0(1 - v / C) => dv = -(C/f0) * df
    dparam_cur_batch[:, n_params - 1] *= C_VAL / f0_batch

    for prune_level in range(1, nsegments):
        coord_next = snail_scheme.get_coord(prune_level)
        coord_cur = snail_scheme.get_current_coord(
            prune_level,
            moving_grid=use_moving_grid,
        )
        _, t_obs_minus_t_ref = coord_cur

        dparam_new_batch = psr_utils.poly_taylor_step_d_vec(
            n_params,
            t_obs_minus_t_ref,
            nbins,
            eta,
            f0_batch,
            t_ref=0,
            use_cheby=use_cheby_coarsening,
        )
        shift_bins_batch = psr_utils.poly_taylor_shift_d_vec(
            dparam_cur_batch,
            dparam_new_batch,
            t_obs_minus_t_ref,
            nbins,
            f0_batch,
            t_ref=0,
            use_cheby=use_cheby_coarsening,
        )
        n_branches = np.ones(n_freqs, dtype=np.int64)

        for i in range(n_freqs):
            for j in range(n_params):  # skip d0
                if shift_bins_batch[i, j] < (eta - FLOAT_EPSILON):
                    dparam_cur_next[i, j] = dparam_cur_batch[i, j]
                    continue
                ratio = dparam_cur_batch[i, j] / dparam_new_batch[i, j]
                num_points = max(1, math.ceil(ratio - FLOAT_EPSILON))
                n_branches[i] *= num_points
                dparam_cur_next[i, j] = dparam_cur_batch[i, j] / num_points
        # Compute average branching factor
        children = np.sum(weights * n_branches)
        parents = np.sum(weights)
        branching_pattern[prune_level - 1] = children / parents
        # Update weights and dparams
        weights *= n_branches

        if use_moving_grid:
            # Transform dparams to the next segment
            delta_t = coord_next[0] - coord_cur[0]
            dparam_d_vec[:, :-1] = dparam_cur_next
            dparam_d_vec_new = transforms.shift_taylor_errors(
                dparam_d_vec,
                delta_t,
                tiling_strategy,
            )
            dparam_cur_batch = dparam_d_vec_new[:, :-1]
        else:
            dparam_cur_batch = dparam_cur_next
    return branching_pattern
