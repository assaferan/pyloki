"""The exact nearest-template search, ported to the Chebyshev basis.

None of the Taylor results transfer. `T(delta_t)` there is lower-triangular with unit
diagonal, which is what made the box strategies' cells tile exactly (D36) and kept
`aggressive`'s diagonal-only tracking nearly harmless (D41). The Chebyshev
interval-change matrix `C` is triangular but its diagonal is `(ts_new/ts_old)^k`, so
under `aggressive` -- which keeps only `|diag(C)|` -- the tracked width is wrong by the
whole off-diagonal part, and that part is large (the first row of `C` at a single
segment doubling runs 2.07, 4.15, 5.64, 7.72, 3.92).

The port is sound for the same structural reason as the Taylor case: the Chebyshev
branch also goes through `branch_param_padded` with `(param_cur, dparam_cur,
dparam_new)`, and `shift_cheby_errors` propagates the error vector without reference to
the values, so every leaf at a stage still shares `dparam_cur` and child offsets are
still parent-independent. The leaf set is therefore again exactly a Minkowski sum.

Two differences that matter for the conventions:

- `poly_cheb_step_vec` is **uniform across axes**, `(eta/N_b) * c/f_max`, since
  `|T_k| <= 1` bounds every coefficient identically. The corner of the nominal cell is
  therefore `k_max/2 * eta/N_b` -- linear in the order, against the Taylor grid's
  `2^(k_max-1) - 1/2` (D46).
- the branch guard uses `poly_cheb_shift_vec`, which is `|dparam_old - dparam_new|`
  rather than the width itself, so an axis branches on the *change* in claimed width.

The phase residual is `dPhi(x) = (f0/c) * sum_k dalpha_k T_k(x)` over `x` in `[-1, 1]`,
the accumulated domain. `alpha_0` is excluded: it is a constant phase offset, which
rotates the folded profile without costing S/N. That is safe here because the last row
of `C` is `e_last` (verified), so `alpha_0` never leaks into a higher coefficient.
"""

from __future__ import annotations

import numpy as np

from pyloki.utils import maths, psr_utils, transforms
from pyloki.utils.misc import C_VAL
from pyloki.utils.snail import MiddleOutScheme

from nearest_template import min_excursion as _min_excursion


def min_excursion(*args, **kwargs):
    """`nearest_template.min_excursion` with the Chebyshev phase basis."""
    kwargs.setdefault("basis_fn", cheby_basis)
    return _min_excursion(*args, **kwargs)


def cheby_basis(x: np.ndarray, poly_order: int) -> np.ndarray:
    """(n_x, poly_order) matrix of T_k(x), columns ordered [alpha_kmax .. alpha_1]."""
    orders = np.arange(poly_order, 0, -1)
    return np.stack([np.polynomial.chebyshev.Chebyshev.basis(int(k))(x)
                     for k in orders], axis=1)


def build_sets_cheby(
    cfg,
    strategy: str,
    poly_order: int,
    f0: float,
    truth_kin: np.ndarray,
    target_stage: int,
    n_seed: int = 1,
) -> tuple[list[np.ndarray], np.ndarray, float]:
    """Offset sets whose Minkowski sum is the Chebyshev leaf set at `target_stage`.

    `truth_kin` is given in the *Taylor* kinematic basis, as elsewhere in these docs, and
    converted here with the same `taylor_to_cheby` the seed uses.
    """
    nseg = int(np.ceil(cfg.nsamps / cfg.bseg_ffa))
    scheme = MiddleOutScheme(nseg, nseg // 2, cfg.tseg_ffa, stride=1)
    coord_init = scheme.get_coord(0)
    scale_init = coord_init[1]

    dparam = np.array(cfg.get_dparams_actual(cfg.niters_ffa), dtype=float)
    dparam[-1] *= C_VAL / f0

    # Seed on the base grid, in Taylor coordinates, then convert exactly as
    # `poly_chebyshev_seed` does. The conversion is linear, so offsets convert too.
    centre = np.round(truth_kin / dparam) * dparam - truth_kin
    axes = [centre[j] + np.arange(-n_seed, n_seed + 1) * dparam[j]
            for j in range(poly_order)]
    seed_t = np.array(np.meshgrid(*axes, indexing="ij")).reshape(poly_order, -1).T
    seed_full = np.zeros((len(seed_t), poly_order + 1))
    seed_full[:, :poly_order] = seed_t
    seed = transforms.taylor_to_cheby(seed_full, scale_init)[:, :poly_order]
    sets = [np.ascontiguousarray(seed)]

    w_full = np.zeros(poly_order + 1)
    w_full[:poly_order] = dparam
    w = transforms.taylor_to_cheby_errors(w_full, scale_init)[:poly_order]

    x = np.linspace(-1.0, 1.0, 256)
    for lvl in range(1, target_stage + 1):
        coord_cur = scheme.get_current_coord(lvl, moving_grid=True)
        coord_prev = scheme.get_previous_coord(lvl, moving_grid=True)

        # Transform into the current domain first -- the Chebyshev branch does this
        # inside itself, unlike the Taylor path which branches then transforms.
        c_mat = maths.poly_chebyshev_transform_matrix(
            poly_order, coord_prev[0], coord_prev[1], coord_cur[0], coord_cur[1], 1)
        c_kin = c_mat[:poly_order, :poly_order]
        sets = [s @ c_kin for s in sets]
        w = transforms.shift_cheby_errors(
            np.concatenate([w, [0.0]])[None, :], coord_cur, coord_prev, strategy,
        )[0][:poly_order]

        step = psr_utils.poly_cheb_step_vec(
            poly_order, cfg.nbins, cfg.eta, np.array([f0]))
        shift = psr_utils.poly_cheb_shift_vec(
            w[None, :], step, cfg.nbins, np.array([f0]))[0]

        axes = []
        for j in range(poly_order):
            if shift[j] >= cfg.eta - 1e-12:
                n = max(1, int(np.ceil(w[j] / step[0, j] - 1e-12)))
                act = w[j] / n
                axes.append(-w[j] / 2.0 + act / 2.0 + np.arange(n) * act)
                w[j] = act
            else:
                axes.append(np.zeros(1))
        if max(len(a) for a in axes) > 1:
            o = np.array(np.meshgrid(*axes, indexing="ij")).reshape(poly_order, -1).T
            sets.append(np.ascontiguousarray(o))

    return sets, x, f0 / C_VAL
