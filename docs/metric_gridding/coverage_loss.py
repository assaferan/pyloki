"""How much of the parameter space does each box strategy actually leave uncovered?

Every other comparison in Phase 3 measures **cost**. This measures the **benefit** side:
the holes a strategy leaves in the searched volume, which is what the whole project is
supposed to fix.

The question has to be asked of the *union* of leaves, not of one leaf. Re-centring the
expansion maps coefficients by `T`, so a signal that leaves its own leaf simply enters a
neighbour -- the grid moves coherently. What matters is whether boxes of the recorded
size, placed on the **sheared** lattice of leaf centres, still tile the space.

`T` is lower-triangular with **unit diagonal** (C4), so:

  aggressive   records `h * |diag(T)| == h`: widths unchanged while the lattice shears
               underneath them. The box no longer contains a fundamental domain of the
               sheared lattice, and the shortfall is the gap.
  quadrature   records `sqrt((h**2) @ (T**2).T)`, wider, still not guaranteed to cover.
  conservative records `h @ |T|.T`, the AABB of the sheared cell: covers by
               construction, at the cost the whole project is trying to avoid.

Method: draw points uniformly from one sheared fundamental cell, and ask whether each is
inside the recorded box of *some* nearby lattice point. The uncovered fraction is the
gap. Sanity checks that must hold, and are asserted in `self_check()`: `conservative`
never leaves a hole, and at zero shear every strategy covers perfectly.
"""

from __future__ import annotations

import itertools

import numpy as np

from pyloki.utils import maths, transforms

RNG = np.random.default_rng(20260914)


def shift_matrix(delta_t: float, n: int) -> np.ndarray:
    """The same `T` that `shift_taylor_errors` builds."""
    powers = np.tril(np.arange(n)[:, None] - np.arange(n))
    return delta_t**powers / maths.fact(powers) * np.tril(np.ones_like(powers))


def uncovered_fraction(
    spacing: np.ndarray,
    delta_t: float,
    strategy: str,
    n_sample: int = 20000,
) -> float:
    """Fraction of a sheared lattice cell left uncovered by the recorded boxes.

    The covering lattice point is **solved for by forward substitution**, not searched.
    Two earlier attempts searched a block of neighbours instead and both reported
    spurious holes: `T` is lower-triangular, so the index that covers a point can be
    arbitrarily far from the nearest one on the lower-order axes -- the offset needed on
    row `i` is driven by `T[i, j] * d_j * s_j / s_i` from every row above it.

    Row `i` of `T @ ((idx - n) * spacing)` involves only `n_0 .. n_i`, so choosing
    `n_0, n_1, ...` in order makes each row's residual exactly the rounding error of one
    integer choice. That is optimal per row and therefore finds a covering point
    whenever one exists.
    """
    n_dim = len(spacing)
    t_mat = shift_matrix(delta_t, n_dim)
    t_inv = np.linalg.inv(t_mat)
    recorded_half = np.abs(
        transforms.shift_taylor_errors(spacing[None, :], delta_t, strategy)[0] / 2.0,
    )

    pts = RNG.uniform(-0.5, 0.5, size=(n_sample, n_dim)) * spacing
    sheared = pts @ t_mat.T
    idx = (sheared @ t_inv.T) / spacing

    # Forward substitution: fix n_i to minimise |row i| given n_0..n_{i-1}.
    n_int = np.zeros_like(idx)
    for i in range(n_dim):
        # residual on row i from the already-chosen indices, in units of spacing[i]
        carry = np.zeros(len(idx))
        for j in range(i):
            carry += t_mat[i, j] * (idx[:, j] - n_int[:, j]) * spacing[j]
        n_int[:, i] = np.round(idx[:, i] + carry / spacing[i])

    centres = (n_int * spacing) @ t_mat.T
    covered = np.all(np.abs(sheared - centres) <= recorded_half + 1e-12, axis=1)
    return float(1.0 - covered.mean())


def self_check() -> None:
    """The two properties that must hold if the measurement means anything."""
    spacing = np.array([1.0, 0.5, 0.25])
    for strategy in ("aggressive", "quadrature", "conservative"):
        zero = uncovered_fraction(spacing, 0.0, strategy, n_sample=4000)
        if zero > 1e-9:
            msg = f"{strategy}: gap {zero} at zero shear, must be 0"
            raise AssertionError(msg)
    cons = uncovered_fraction(spacing, 0.7, "conservative", n_sample=4000)
    if cons > 1e-9:
        msg = f"conservative leaves a hole ({cons}); it is the AABB and cannot"
        raise AssertionError(msg)
