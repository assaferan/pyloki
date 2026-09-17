"""Correctness of the exact nearest-template search (`docs/metric_gridding`).

These pin the three properties whose absence made `pruning_multiplicity.py` unsound
(D47): the structural claim the method rests on, exact agreement with brute force, and
monotonicity of the reported minimum in the search width.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

from pyloki.utils import psr_utils

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "docs" / "metric_gridding"))

import phase3_config as P  # noqa: E402
from nearest_template import _basis, build_sets, min_excursion  # noqa: E402

F0 = 1.0 / P.PERIOD
PO = P.POLY_ORDER
TRUTH = np.array([0.001, 0.05, 1.0, 0.0])


@pytest.mark.parametrize(
    ("c1", "c2", "dcur", "dnew"),
    [(0.0, 7.3, 4.0, 0.9), (-11.2, 250.0, 16.0, 1.1), (3.0, 3.0, 2.0, 5.0)],
)
def test_offsets_are_parent_independent(c1, c2, dcur, dnew):
    """The whole method rests on this: child offsets do not depend on the parent.

    `branch_param_padded` derives them from `(dparam_cur, dparam_new)` alone, and every
    leaf at a stage shares those, which is what makes the leaf set an exact Minkowski
    sum of per-stage offset sets rather than something that has to be materialised.
    """
    o1 = np.zeros(32)
    o2 = np.zeros(32)
    a1, n1 = psr_utils.branch_param_padded(o1, c1, dcur, dnew)
    a2, n2 = psr_utils.branch_param_padded(o2, c2, dcur, dnew)
    assert n1 == n2
    assert a1 == pytest.approx(a2)
    np.testing.assert_allclose(o1[:n1] - c1, o2[:n2] - c2, atol=1e-12)


@pytest.mark.parametrize("strategy", ["aggressive", "quadrature", "conservative"])
@pytest.mark.parametrize(("stage", "n_seed"), [(4, 0), (6, 0), (5, 1)])
def test_branch_and_bound_matches_brute_force(strategy, stage, n_seed):
    """The bound never discards the optimum: B&B equals exhaustive enumeration."""
    cfg = P.make_config(strategy)
    sets, tau, scale = build_sets(cfg, strategy, PO, F0, TRUTH, stage, n_seed=n_seed)
    total = float(np.prod([float(len(s)) for s in sets]))
    if total > 4e6:
        pytest.skip(f"{total:.2g} leaves is too many to enumerate")

    value, exact, _ = min_excursion(sets, tau, scale, PO, cfg.eta / cfg.nbins)
    assert exact

    acc = sets[0]
    for nxt in sets[1:]:
        acc = (acc[:, None, :] + nxt[None, :, :]).reshape(-1, PO)
    basis = _basis(tau, PO) * (scale / (cfg.eta / cfg.nbins))
    brute = float(np.abs(acc @ basis.T).max(axis=1).min())
    assert value == pytest.approx(brute, rel=1e-12, abs=1e-12)


@pytest.mark.parametrize("strategy", ["aggressive", "quadrature"])
def test_minimum_is_monotone_in_seed_width(strategy):
    """A wider seed only adds candidates, so the minimum must not increase.

    `pruning_multiplicity.py` failed exactly here -- its keep-best-N cap made the
    reported minimum rise with the search width, which a minimum cannot do.
    """
    cfg = P.make_config(strategy)
    values = []
    for n_seed in (0, 1, 2):
        sets, tau, scale = build_sets(
            cfg, strategy, PO, F0, TRUTH, 10, n_seed=n_seed)
        value, exact, _ = min_excursion(sets, tau, scale, PO, cfg.eta / cfg.nbins)
        assert exact
        values.append(value)
    assert values[1] <= values[0] + 1e-12
    assert values[2] <= values[1] + 1e-12


def test_seeded_incumbent_gives_sound_verdicts():
    """`init_best` turns the search into a decision procedure without changing it."""
    cfg = P.make_config("quadrature")
    sets, tau, scale = build_sets(cfg, "quadrature", PO, F0, TRUTH, 8, n_seed=0)
    tol = cfg.eta / cfg.nbins
    true_min, exact, _ = min_excursion(sets, tau, scale, PO, tol)
    assert exact

    # Seeded above the optimum: must still find it.
    v_hi, ex_hi, _ = min_excursion(sets, tau, scale, PO, tol,
                                   init_best=true_min * 2.0)
    assert ex_hi
    assert v_hi == pytest.approx(true_min)

    # Seeded below the optimum: must prove no leaf beats it, and not invent one.
    target = true_min * 0.5
    v_lo, ex_lo, _ = min_excursion(sets, tau, scale, PO, tol, init_best=target)
    assert ex_lo
    assert v_lo == pytest.approx(target)


# --- amplitude conversion -------------------------------------------------------

from amplitude_loss import snr_ratio  # noqa: E402

from pyloki.utils.misc import C_VAL  # noqa: E402

TAU = np.linspace(-10.0, 10.0, 2048)


def _ramp(x: float) -> np.ndarray:
    """A pure d_1 offset whose phase residual has sup-norm x * (eta/N_b) cycles."""
    d = np.zeros(PO)
    d[-1] = x * (1.0 / P.NBINS) * C_VAL / F0 / np.abs(TAU).max()
    return d


@pytest.mark.parametrize("filt", ["boxcar", "matched"])
def test_zero_phase_error_costs_nothing(filt):
    assert snr_ratio(np.zeros(PO), TAU, F0, PO, P.NBINS, 0.1, filt=filt) == (
        pytest.approx(1.0, abs=1e-12))


@pytest.mark.parametrize("filt", ["boxcar", "matched"])
@pytest.mark.parametrize("ducy", [0.05, 0.1, 0.2])
def test_loss_is_monotone_in_phase_error(filt, ducy):
    """More smearing cannot raise the recovered S/N.

    Before the phase average was added this failed: a Gaussian smeared towards a flat
    top scores *better* against a boxcar bank than a sharp one, which produced negative
    losses of the same size as the effect under study.
    """
    ratios = [
        snr_ratio(_ramp(x), TAU, F0, PO, P.NBINS, ducy, filt=filt)
        for x in (0.0, 0.5, 1.0, 2.0, 3.0)
    ]
    for lo, hi in zip(ratios[:-1], ratios[1:], strict=True):
        assert hi <= lo + 1e-9


@pytest.mark.parametrize("ducy", [0.1, 0.2])
def test_boxcar_and_matched_filters_agree(ducy):
    """The loss must be a property of the grid, not of the filter used to read it."""
    for x in (1.0, 2.0):
        box = 1.0 - snr_ratio(_ramp(x), TAU, F0, PO, P.NBINS, ducy, filt="boxcar")
        mat = 1.0 - snr_ratio(_ramp(x), TAU, F0, PO, P.NBINS, ducy, filt="matched")
        assert box == pytest.approx(mat, rel=0.35)


def test_narrower_pulses_lose_more():
    """Smearing costs a narrow pulse more than a broad one, at equal phase error."""
    losses = [
        1.0 - snr_ratio(_ramp(1.0), TAU, F0, PO, P.NBINS, d, filt="matched")
        for d in (0.20, 0.10, 0.05)
    ]
    assert losses[0] < losses[1] < losses[2]


# --- Chebyshev port -------------------------------------------------------------

from nearest_template_cheby import (  # noqa: E402
    build_sets_cheby,
    cheby_basis,
    min_excursion as min_excursion_cheby,
)


@pytest.mark.parametrize("strategy", ["aggressive", "quadrature", "conservative"])
@pytest.mark.parametrize(("stage", "n_seed"), [(4, 0), (6, 0)])
def test_cheby_branch_and_bound_matches_brute_force(strategy, stage, n_seed):
    """The port keeps the guarantee: B&B equals exhaustive enumeration."""
    cfg = P.make_config(strategy)
    sets, x, scale = build_sets_cheby(
        cfg, strategy, PO, F0, TRUTH, stage, n_seed=n_seed)
    total = float(np.prod([float(len(s)) for s in sets]))
    if total > 3e6:
        pytest.skip(f"{total:.2g} leaves is too many to enumerate")

    value, exact, _ = min_excursion_cheby(sets, x, scale, PO, cfg.eta / cfg.nbins)
    assert exact

    acc = sets[0]
    for nxt in sets[1:]:
        acc = (acc[:, None, :] + nxt[None, :, :]).reshape(-1, PO)
    basis = cheby_basis(x, PO) * (scale / (cfg.eta / cfg.nbins))
    brute = float(np.abs(acc @ basis.T).max(axis=1).min())
    assert value == pytest.approx(brute, rel=1e-12, abs=1e-12)


def test_cheby_step_is_uniform_and_corner_is_kmax_over_two():
    """|T_k| <= 1 bounds every coefficient alike, so the corner is linear in the order.

    Against the Taylor grid's `2**(k-1) - 1/2` (D46), which is geometric. Both are
    nominal-cell corners and neither is a covering radius (D48).
    """
    for k_max in range(2, 9):
        step = psr_utils.poly_cheb_step_vec(
            k_max, P.NBINS, 1.0, np.array([F0]))[0]
        assert np.allclose(step, step[0])
        corner = (F0 / C_VAL) * np.sum(step / 2.0) * P.NBINS
        assert corner == pytest.approx(k_max / 2.0, rel=1e-12)


def test_cheby_transform_diagonal_is_not_unity():
    """The reason no Taylor tiling result transfers: diag(C) = (ts_new/ts_old)**k."""
    from pyloki.utils import maths
    ts1, ts2 = 10.0, 20.0
    c_mat = maths.poly_chebyshev_transform_matrix(4, 0.0, ts1, 0.0, ts2, 1)
    expected = (ts2 / ts1) ** np.arange(4, -1, -1)
    np.testing.assert_allclose(np.diag(c_mat), expected, rtol=1e-12)
    assert not np.allclose(np.diag(c_mat), 1.0)
    # alpha_0 is a constant phase offset and never leaks into a higher coefficient,
    # which is what lets the search drop it from the excursion basis.
    np.testing.assert_allclose(c_mat[-1], np.eye(5)[-1], atol=1e-12)
