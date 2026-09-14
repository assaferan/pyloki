"""Quantify the sensitivity loss of each tiling strategy.

This is the measurement Kumar & Zackay (2026) sec. 5.2.4 actually asks for --
"Further work is required to quantify the sensitivity loss and compare basis
strategies" -- and it is not the coverage question. D36 established that the union of
recorded tiles covers the space exactly, so no signal is *missed*. What can still go
wrong is that a signal sits further from its template centre than the grid promised.

The promise (sec. 3.4, eq. `grid_criteria`) is a **sup-norm** bound:

    |dPhi(t)| <= eta / N_b   for all t in the accumulated interval.

The grid is built by bounding each Taylor coefficient *independently* (eq.
`grid_params`), so a signal offset in one coefficient alone respects it by construction.
A signal offset in **all** coefficients at once -- a corner of the cell -- need not: the
per-axis contributions add. That is the "corner" loss the paper's Figure `grid_tiling`
is about, and how big it is depends on the cell each strategy transports.

So: walk the stages, and at each one sample offsets uniformly inside the leaf's own cell
and evaluate

    dPhi(tau) = (f0 / c) * sum_k  delta_k * tau**k / k!        [cycles]

over the true window, in units of the tolerance `eta / N_b`. A value of 1 means the
grid delivered exactly what it promised; 10 means a signal there loses ten times the
phase budget.
"""

from __future__ import annotations

import numpy as np

from pyloki.core import metric
from pyloki.utils import maths, psr_utils, transforms
from pyloki.utils.misc import C_VAL
from pyloki.utils.snail import MiddleOutScheme

RNG = np.random.default_rng(20260914)


def phase_excursion(
    delta: np.ndarray, tau: np.ndarray, f0: float, poly_order: int,
) -> np.ndarray:
    """max_t |dPhi(t)| in cycles, for coefficient offsets `delta` (leaf axis order)."""
    # Leaf order is [d_kmax .. d_2, d_1]; d_1 is the k=1 term.
    orders = np.arange(poly_order, 0, -1)
    basis = tau[None, :] ** orders[:, None] / maths.fact(orders)[:, None]
    return np.abs((f0 / C_VAL) * (delta @ basis)).max(axis=1)


def walk_stages(
    cfg,
    strategy: str,
    poly_order: int,
    f0: float,
    n_sample: int = 4000,
    n_time: int = 512,
) -> dict[str, np.ndarray]:
    """Per-stage sensitivity loss, in units of the promised tolerance."""
    nseg = int(np.ceil(cfg.nsamps / cfg.bseg_ffa))
    scheme = MiddleOutScheme(nseg, nseg // 2, cfg.tseg_ffa, stride=1)
    tol = cfg.eta / cfg.nbins  # the promise, in cycles

    dparam = np.array(cfg.get_dparams_actual(cfg.niters_ffa), dtype=float)
    dparam[-1] *= C_VAL / f0
    vec = np.zeros(poly_order + 1)
    vec[:poly_order] = dparam

    median, worst, mismatch = [], [], []
    for lvl in range(1, nseg):
        ref_cur, t_half = scheme.get_current_coord(lvl, moving_grid=True)
        ref_next, _ = scheme.get_coord(lvl)
        delta_t = ref_next - ref_cur

        # Branch exactly as the box strategies do.
        step = psr_utils.poly_taylor_step_d_vec(
            poly_order, t_half, cfg.nbins, cfg.eta, np.array([f0]), t_ref=0)
        shift = psr_utils.poly_taylor_shift_d_vec(
            vec[None, :poly_order], step, t_half, cfg.nbins, np.array([f0]), t_ref=0)[0]
        npts = np.array([
            max(1, int(np.ceil(vec[j] / step[0, j] - 1e-12)))
            if shift[j] >= cfg.eta - 1e-12 else 1
            for j in range(poly_order)], dtype=float)
        vec[:poly_order] /= npts

        # The window the leaf is actually valid over, relative to its own epoch (D20):
        # the leaf is expanded about the previous centre, the data has moved on by
        # delta_t, and spans +/- t_half about the new centre.
        tau = np.linspace(delta_t - t_half, delta_t + t_half, n_time)
        half = vec[:poly_order] / 2.0
        offsets = RNG.uniform(-1.0, 1.0, size=(n_sample, poly_order)) * half

        excursion = phase_excursion(offsets, tau, f0, poly_order) / tol
        median.append(float(np.median(excursion)))
        worst.append(float(excursion.max()))

        g = metric.poly_phase_metric(
            0.0, delta_t - t_half, delta_t + t_half, poly_order, f0,
            cfg.nbins, cfg.ducy_max)
        mismatch.append(float(np.median(np.einsum("ni,ij,nj->n", offsets, g, offsets))))

        vec = transforms.shift_taylor_errors(vec[None, :], delta_t, strategy)[0]

    return {
        "median": np.asarray(median),
        "worst": np.asarray(worst),
        "mismatch": np.asarray(mismatch),
    }


def self_check(cfg, poly_order: int, f0: float) -> None:
    """Zero costs nothing; the naive grid keeps its promise; the shipped one does not.

    The third assertion is the point of D38. `poly_taylor_step_d_vec` defaults to
    `use_cheby=True`, applying the paper's `2**(k-1)` coarsening (sec. 3.4, appendix
    `app:optimal_gridding`). That factor comes from diagonalising the parameter metric
    in the Chebyshev basis, where the bound is imposed on the *orthogonal* coefficients
    -- so the shipped grid deliberately exceeds the naive per-axis sup-norm bound, by
    exactly `2**(k-1)` on the order-`k` axis. Worth pinning, because it means `eta` does
    not mean "phase error stays under eta/N_b per axis".
    """
    tau = np.linspace(-10.0, 10.0, 256)
    zero = phase_excursion(np.zeros((1, poly_order)), tau, f0, poly_order)
    if zero[0] != 0.0:
        msg = f"zero offset gave {zero[0]}"
        raise AssertionError(msg)

    t_half, tol = 10.0, cfg.eta / cfg.nbins
    fine = np.linspace(0.0, t_half, 2048)

    def one_axis(axis: int, *, use_cheby: bool, frac: float = 1.0) -> float:
        step = psr_utils.poly_taylor_step_d_vec(
            poly_order, t_half, cfg.nbins, cfg.eta, np.array([f0]),
            t_ref=0, use_cheby=use_cheby)[0]
        d = np.zeros((1, poly_order))
        d[0, axis] = step[axis] * frac
        return float(phase_excursion(d, fine, f0, poly_order)[0] / tol)

    # The naive grid is calibrated so a FULL step offset costs exactly the tolerance.
    # (Column 1 is a full span, D24, so a cell corner sits at half of this per axis.)
    for axis in range(poly_order):
        got = one_axis(axis, use_cheby=False)
        if not np.isclose(got, 1.0, rtol=1e-6):
            msg = f"naive grid, axis {axis}: full step cost {got:.4f}, expected 1"
            raise AssertionError(msg)

    # The shipped grid is coarser by exactly 2**(k-1) on the order-k axis (D38).
    # Leaf order is [d_kmax .. d_1], so axis `i` carries order k = poly_order - i.
    for axis in range(poly_order):
        k = poly_order - axis
        ratio = one_axis(axis, use_cheby=True) / one_axis(axis, use_cheby=False)
        if not np.isclose(ratio, 2.0 ** (k - 1), rtol=1e-6):
            msg = (
                f"shipped grid, axis {axis} (order {k}): coarsening {ratio:.4f}, "
                f"expected 2**(k-1) = {2.0 ** (k - 1)}"
            )
            raise AssertionError(msg)
