"""Phase 2 tests: the metric covering and the metric branch function.

`metric_PLAN.md`'s Phase 2 exit criterion is "coverage test green; redundancy known;
`aggressive` path untouched", and the three parts are tested separately here:

- `TestCoverage` is the guarantee. It is the reason `_retention_form` exists, and the
  reason it uses a *sound* outer ellipsoid of a Minkowski sum rather than the tighter
  but unsound per-axis dilation (see the docstring there).
- `TestRedundancy` pins the measured overhead against its predicted decomposition, so
  a regression in the retention bound shows up as a number rather than as a slowdown.
- `TestAggressiveUntouched` and `TestBranchContract` guard the shipped path and the
  leaf-array contract.
"""

from __future__ import annotations

import logging
from types import SimpleNamespace

import numpy as np
import pytest
from numba import njit, prange
from scipy.special import gamma

from pyloki.config import ParamLimits, PulsarSearchConfig
from pyloki.core import metric, taylor
from pyloki.dynamic import PrunePolyTaylorDPFuncts
from pyloki.prune import Pruning
from pyloki.utils import psr_utils, transforms
from pyloki.utils.snail import MiddleOutScheme

RNG = np.random.default_rng(20260910)
F0 = 1.0 / 0.007  # 142.857 Hz, the spin period used throughout the plan
NBINS = 64
DUCY = 0.1
M_MAX = 0.2
# One EP stage: the accumulated segment doubles.
T_PARENT = 16.78
T_CHILD = 33.55


def metrics_for(poly_order: int, f0: float = F0) -> tuple[np.ndarray, np.ndarray]:
    """Parent- and child-stage metrics for one segment doubling."""
    g_parent = metric.poly_phase_metric(
        0.0, 0.0, T_PARENT, poly_order, f0, NBINS, DUCY,
    )
    g_child = metric.poly_phase_metric(0.0, 0.0, T_CHILD, poly_order, f0, NBINS, DUCY)
    return g_parent, g_child


def sample_parent_region(
    g_parent: np.ndarray, m_max: float, n_samp: int,
) -> np.ndarray:
    """Uniform samples of `{d : d^T g_parent d <= m_max}`.

    Uniform *in the region*, not in the whitened ball's radius: the radius has to be
    drawn as `u**(1/n)` or the samples pile up at the centre and the corners of the
    region -- exactly where coverage is hardest -- go untested.
    """
    n_dim = g_parent.shape[0]
    raw = RNG.normal(size=(n_samp, n_dim))
    unit = raw / np.linalg.norm(raw, axis=1, keepdims=True)
    radii = RNG.random(n_samp) ** (1.0 / n_dim)
    chol = metric.cholesky_factor(g_parent)
    return np.sqrt(m_max) * np.linalg.solve(chol.T, (unit * radii[:, None]).T).T


def cover_stats(
    samples: np.ndarray, offsets: np.ndarray, g_child: np.ndarray, m_max: float,
) -> tuple[int, float, float]:
    """(uncovered count, mean children covering a sample, worst m / m_max).

    Chunked over samples: the full `(n_samp, n_child, n_dim)` difference is ~1 GB at
    `poly_order=4`, which is an OOM rather than a slow test.
    """
    n_child, n_dim = offsets.shape
    chunk = max(1, int(4e7 / max(n_child * n_dim, 1)))
    uncovered = 0
    covers = []
    worst = 0.0
    for start in range(0, len(samples), chunk):
        blk = samples[start : start + chunk]
        diff = blk[:, None, :] - offsets[None, :, :]
        m_to = np.einsum("sci,ij,scj->sc", diff, g_child, diff)
        n_cov = (m_to <= m_max * (1.0 + 1e-9)).sum(axis=1)
        uncovered += int((n_cov == 0).sum())
        covers.append(n_cov)
        worst = max(worst, float(np.max(np.min(m_to, axis=1)) / m_max))
    return uncovered, float(np.concatenate(covers).mean()), worst


@njit(cache=True, fastmath=True)
def njit_caller(
    leaves_batch: np.ndarray, offsets: np.ndarray, extents: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Stand-in for `dyn_poly_taylor.branch_func`, which is `@njit` too."""
    return taylor.poly_taylor_branch_metric_apply(leaves_batch, offsets, extents)


@njit(cache=True, fastmath=True, parallel=True)
def njit_parallel_caller(
    leaves_batch: np.ndarray,
    offsets: np.ndarray,
    extents: np.ndarray,
    n_rep: int,
) -> float:
    """Pruning branches inside a `prange`; nested regions have bitten before."""
    total = 0.0
    for _ in prange(n_rep):
        out, _origins = taylor.poly_taylor_branch_metric_apply(
            leaves_batch, offsets, extents,
        )
        total += out[0, 0, 0]
    return total


def cubic_thickness(n_dim: int) -> float:
    """Covering thickness of `Z^n` at covering radius `a * sqrt(n) / 2`."""
    v_n = np.pi ** (n_dim / 2) / gamma(n_dim / 2 + 1)
    return float(v_n * (np.sqrt(n_dim) / 2) ** n_dim)


class TestLatticeSpacing:
    def test_covering_radius_is_sqrt_m_max(self) -> None:
        """A cell corner is the worst case, and must land exactly on the budget."""
        for n_dim in (1, 2, 3, 4, 5):
            a = metric.cubic_lattice_spacing(n_dim, M_MAX)
            corner = a * np.sqrt(n_dim) / 2.0
            assert corner == pytest.approx(np.sqrt(M_MAX))

    def test_rejects_bad_arguments(self) -> None:
        with pytest.raises(ValueError, match="n_dim"):
            metric.cubic_lattice_spacing(0, M_MAX)
        with pytest.raises(ValueError, match="m_max"):
            metric.cubic_lattice_spacing(2, 0.0)


class TestRetentionForm:
    """The retention set must be a superset of `E_parent + B(sqrt(m_max))`."""

    @pytest.mark.parametrize("poly_order", [2, 3, 4])
    def test_contains_minkowski_sum(self, poly_order: int) -> None:
        g_parent, g_child = metrics_for(poly_order)
        chol = metric.cholesky_factor(g_child)
        a_mat = np.linalg.solve(chol, np.linalg.solve(chol, g_parent).T).T
        a_mat = 0.5 * (a_mat + a_mat.T)
        form = metric._retention_form(a_mat, M_MAX)

        # Sample the boundary of E + B(r): a parent-boundary point plus a radius-r
        # kick. The boundary is where the containment is tight.
        lam, evec = np.linalg.eigh(a_mat)
        a_semi = np.sqrt(M_MAX / lam)
        r = np.sqrt(M_MAX)
        n = 50_000
        x = RNG.normal(size=(n, poly_order))
        x /= np.linalg.norm(x, axis=1, keepdims=True)
        x = ((evec * a_semi) @ x.T).T
        u = RNG.normal(size=(n, poly_order))
        u = r * u / np.linalg.norm(u, axis=1, keepdims=True)
        w = x + u

        value = np.einsum("si,ij,sj->s", w, form, w)
        assert value.max() <= 1.0 + 1e-9, "retention form does not contain E + B(r)"
        # Non-vacuity: it must also be tight, not an arbitrarily huge ellipsoid.
        assert value.max() > 0.5, f"retention form is loose (max {value.max():.3f})"

    @pytest.mark.parametrize("poly_order", [2, 3, 4])
    def test_per_axis_dilation_would_be_unsound(self, poly_order: int) -> None:
        """Guards the reasoning, not the code: records *why* the cheap form is wrong.

        Dilating each semi-axis `a_i -> a_i + r` is the obvious tightening and is
        unsound, because comparing support functions needs Cauchy-Schwarz in the
        direction it does not hold. If this ever stops failing, the argument in
        `_retention_form` needs revisiting.
        """
        g_parent, g_child = metrics_for(poly_order)
        chol = metric.cholesky_factor(g_child)
        a_mat = np.linalg.solve(chol, np.linalg.solve(chol, g_parent).T).T
        a_mat = 0.5 * (a_mat + a_mat.T)
        lam, evec = np.linalg.eigh(a_mat)
        a_semi = np.sqrt(M_MAX / lam)
        r = np.sqrt(M_MAX)

        n = 50_000
        x = RNG.normal(size=(n, poly_order))
        x /= np.linalg.norm(x, axis=1, keepdims=True)
        x = ((evec * a_semi) @ x.T).T
        u = RNG.normal(size=(n, poly_order))
        u = r * u / np.linalg.norm(u, axis=1, keepdims=True)
        coords = (x + u) @ evec
        value = (coords**2 / (a_semi + r) ** 2).sum(axis=1)
        assert value.max() > 1.0 + 1e-6, (
            "per-axis dilation now looks sound; re-derive _retention_form"
        )


class TestCoverage:
    """The Phase 2 exit criterion: the covering must actually cover."""

    @pytest.mark.parametrize("poly_order", [2, 3, 4])
    def test_every_parent_point_is_covered(self, poly_order: int) -> None:
        g_parent, g_child = metrics_for(poly_order)
        offsets = metric.lattice_children(
            g_parent, g_child, M_MAX, max_children=500_000,
        )
        samples = sample_parent_region(g_parent, M_MAX, 4000)
        # Guard: the samples really are in the parent region, else this proves nothing.
        assert np.all(metric.mismatch(g_parent, samples) <= M_MAX * (1.0 + 1e-9))

        uncovered, mean_cover, worst = cover_stats(
            samples, offsets, g_child, M_MAX,
        )
        assert uncovered == 0, f"{uncovered}/{len(samples)} parent points uncovered"
        assert worst <= 1.0, f"worst mismatch {worst:.4f} exceeds m_max"
        # Non-vacuity: a trivially huge child count would also pass the above.
        assert mean_cover >= 1.0

    @pytest.mark.parametrize("poly_order", [2, 3])
    def test_offsets_are_a_lattice_in_whitened_coordinates(
        self, poly_order: int,
    ) -> None:
        """Offsets must sit on `spacing * Z^n` after whitening by the child metric."""
        g_parent, g_child = metrics_for(poly_order)
        offsets = metric.lattice_children(g_parent, g_child, M_MAX)
        chol = metric.cholesky_factor(g_child)
        whitened = offsets @ chol
        spacing = metric.cubic_lattice_spacing(poly_order, M_MAX)
        resid = whitened / spacing - np.round(whitened / spacing)
        assert np.max(np.abs(resid)) < 1e-9

    def test_always_returns_at_least_one_child(self) -> None:
        """A parent much smaller than a child still has to produce its own centre."""
        g_parent, _ = metrics_for(2)
        offsets = metric.lattice_children(g_parent, g_parent * 1e6, M_MAX)
        assert len(offsets) >= 1

    def test_cap_raises_rather_than_under_covering(self) -> None:
        g_parent, g_child = metrics_for(4)
        with pytest.raises(ValueError, match="lattice points"):
            metric.lattice_children(g_parent, g_child, M_MAX, max_children=10)


class TestRedundancy:
    """Redundancy is *known*, i.e. pinned to a predicted decomposition."""

    @pytest.mark.parametrize("poly_order", [2, 3, 4])
    def test_mean_cover_matches_cubic_thickness(self, poly_order: int) -> None:
        """Mean coverings per point should be the lattice's covering thickness.

        This is what makes the overhead explainable rather than mysterious: the
        redundancy in the interior is a property of `Z^n`, independent of the region.
        """
        g_parent, g_child = metrics_for(poly_order)
        offsets = metric.lattice_children(
            g_parent, g_child, M_MAX, max_children=500_000,
        )
        samples = sample_parent_region(g_parent, M_MAX, 4000)
        _, mean_cover, _ = cover_stats(samples, offsets, g_child, M_MAX)
        assert mean_cover == pytest.approx(cubic_thickness(poly_order), rel=0.15)

    @pytest.mark.parametrize(
        ("poly_order", "max_overhead"),
        [(2, 8.0), (3, 20.0), (4, 60.0)],
    )
    def test_child_count_overhead_over_volume_bound(
        self, poly_order: int, max_overhead: float,
    ) -> None:
        """Child count vs the packing lower bound `sqrt(det g_child / det g_parent)`.

        No covering can use fewer children than the volume ratio. The bounds here are
        the measured overheads (4.9 / 13.6 / 45.0) with headroom; they exist to catch a
        regression in `_retention_form`, which previously cost 6.4 / 35.9 / 166.7.
        """
        g_parent, g_child = metrics_for(poly_order)
        offsets = metric.lattice_children(
            g_parent, g_child, M_MAX, max_children=500_000,
        )
        _, logdet_p = np.linalg.slogdet(g_parent)
        _, logdet_c = np.linalg.slogdet(g_child)
        vol_bound = float(np.exp(0.5 * (logdet_c - logdet_p)))
        overhead = len(offsets) / vol_bound
        assert overhead >= 1.0, "fewer children than the volume bound is impossible"
        assert overhead <= max_overhead, (
            f"poly_order={poly_order}: overhead {overhead:.1f} over the volume bound "
            f"{vol_bound:.1f} exceeds {max_overhead}"
        )

    def test_metric_covering_costs_more_than_the_box_strategy(self) -> None:
        """Records the Phase 3 decision point rather than asserting it is fine.

        A *complete* covering cannot beat the volume bound, and the box strategy
        undercovers by roughly a factor two concentrated on a single axis, so the metric
        strategy is intrinsically more expensive per parent. This test documents the
        gap so the Phase 3 choice is made against a measured number.
        """
        poly_order, eta = 4, 1.0
        dpar_parent = psr_utils.poly_taylor_step_d_vec(
            poly_order, T_PARENT, NBINS, eta, np.array([F0]), t_ref=0,
        )
        leaves = np.zeros((1, poly_order + 2, 2))
        leaves[0, :-2, 1] = dpar_parent[0]
        leaves[0, -1, 0] = F0
        box, _ = taylor.poly_taylor_branch_batch(
            leaves, (0.0, T_CHILD), NBINS, eta, poly_order, 64,
        )
        g_parent, g_child = metrics_for(poly_order)
        met = metric.lattice_children(
            g_parent, g_child, M_MAX, max_children=500_000,
        )
        assert len(box) < len(met), (
            "metric covering is no longer more expensive than the box strategy; "
            "the Phase 3 lattice decision should be revisited"
        )


class TestBranchContract:
    """`poly_taylor_branch_metric_batch` must be drop-in for the box version."""

    def _leaves(self, n_leaves: int, poly_order: int) -> np.ndarray:
        leaves = np.zeros((n_leaves, poly_order + 2, 2))
        leaves[:, :-2, 0] = RNG.normal(scale=1e-3, size=(n_leaves, poly_order))
        leaves[:, -2, 0] = RNG.normal(size=n_leaves)
        leaves[:, -1, 0] = F0
        leaves[:, -1, 1] = 0.0
        return leaves

    @pytest.mark.parametrize("poly_order", [2, 3])
    def test_shapes_and_origins(self, poly_order: int) -> None:
        n_leaves = 5
        leaves = self._leaves(n_leaves, poly_order)
        out, origins = taylor.poly_taylor_branch_metric_batch(
            leaves, (0.0, T_CHILD), (0.0, T_PARENT), NBINS, DUCY, poly_order,
            M_MAX, 500_000,
        )
        assert out.ndim == 3
        assert out.shape[1:] == (poly_order + 2, 2)
        assert out.shape[0] == len(origins)
        assert origins.min() == 0
        assert origins.max() == n_leaves - 1
        # Every parent must be represented, and equally.
        counts = np.bincount(origins, minlength=n_leaves)
        assert np.all(counts == counts[0])
        assert counts[0] > 1, "branching produced a single child; test is vacuous"

    @pytest.mark.parametrize("poly_order", [2, 3])
    def test_passthrough_columns(self, poly_order: int) -> None:
        leaves = self._leaves(4, poly_order)
        out, origins = taylor.poly_taylor_branch_metric_batch(
            leaves, (0.0, T_CHILD), (0.0, T_PARENT), NBINS, DUCY, poly_order,
            M_MAX, 500_000,
        )
        np.testing.assert_allclose(out[:, -2, 0], leaves[origins, -2, 0])
        np.testing.assert_allclose(out[:, -1, 0], leaves[origins, -1, 0])
        np.testing.assert_allclose(out[:, -1, 1], leaves[origins, -1, 1])

    @pytest.mark.parametrize("poly_order", [2, 3])
    def test_children_cover_their_own_parent(self, poly_order: int) -> None:
        """End-to-end: the branch output covers each parent's region.

        Distinct `f0` per leaf, so this also exercises the `1 / f0` rescaling that lets
        one enumeration serve the whole batch.
        """
        leaves = self._leaves(3, poly_order)
        leaves[:, -1, 0] = np.array([F0, 2.0 * F0, 0.5 * F0])
        out, origins = taylor.poly_taylor_branch_metric_batch(
            leaves, (0.0, T_CHILD), (0.0, T_PARENT), NBINS, DUCY, poly_order,
            M_MAX, 500_000,
        )
        for i in range(len(leaves)):
            f0 = leaves[i, -1, 0]
            g_parent, g_child = metrics_for(poly_order, f0=f0)
            centre = leaves[i, :-2, 0]
            offsets = out[origins == i, :-2, 0] - centre
            samples = sample_parent_region(g_parent, M_MAX, 1500)
            uncovered, mean_cover, worst = cover_stats(
                samples, offsets, g_child, M_MAX,
            )
            assert uncovered == 0, f"leaf {i} (f0={f0:.1f}): {uncovered} uncovered"
            assert worst <= 1.0
            assert mean_cover >= 1.0

    def test_dparam_column_is_the_child_ellipsoid_box(self) -> None:
        poly_order = 3
        leaves = self._leaves(2, poly_order)
        out, _ = taylor.poly_taylor_branch_metric_batch(
            leaves, (0.0, T_CHILD), (0.0, T_PARENT), NBINS, DUCY, poly_order,
            M_MAX, 500_000,
        )
        _, g_child = metrics_for(poly_order)
        expected = metric.ellipsoid_axis_extents(g_child, M_MAX)
        np.testing.assert_allclose(out[0, :-2, 1], expected, rtol=1e-10)
        assert np.all(out[:, :-2, 1] > 0.0)

    def test_cap_propagates(self) -> None:
        leaves = self._leaves(2, 3)
        with pytest.raises(ValueError, match="lattice points"):
            taylor.poly_taylor_branch_metric_batch(
                leaves, (0.0, T_CHILD), (0.0, T_PARENT), NBINS, DUCY, 3, M_MAX, 5,
            )


class TestNjitDispatch:
    """The Phase 2 step 6 unblock: an `@njit` caller must be able to branch.

    `dyn_poly_taylor.branch_func` is `@njit(cache=True, fastmath=True)`, so the
    original pure-NumPy branch could not be dispatched from it. The branch is now split
    into `metric_branch_tables` (Python, once per stage: eigh, Cholesky, Fincke-Pohst)
    and `poly_taylor_branch_metric_apply` (`@njit`, once per batch: a broadcast add and
    a 1/f0 rescale). These tests pin that the njit half really is reachable.
    """

    def _tables_and_leaves(
        self, poly_order: int,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        offsets, extents = taylor.metric_branch_tables(
            T_PARENT, T_CHILD, NBINS, DUCY, poly_order, M_MAX, 500_000,
        )
        leaves = np.zeros((4, poly_order + 2, 2))
        leaves[:, :-2, 0] = RNG.normal(scale=1e-3, size=(4, poly_order))
        leaves[:, -2, 0] = np.arange(4.0)
        leaves[:, -1, 0] = F0
        leaves[1, -1, 0] = 2.0 * F0  # distinct f0, exercising the rescale
        return offsets, extents, leaves

    def test_callable_from_an_njit_function(self) -> None:
        offsets, extents, leaves = self._tables_and_leaves(3)
        out, origins = njit_caller(leaves, offsets, extents)
        assert out.shape == (len(leaves) * len(offsets), 3 + 2, 2)
        assert len(origins) == len(out)
        # Guard: compiled in nopython mode, not an object-mode fallback.
        assert njit_caller.nopython_signatures

    def test_callable_inside_a_prange(self) -> None:
        offsets, extents, leaves = self._tables_and_leaves(2)
        value = njit_parallel_caller(leaves, offsets, extents, 4)
        assert np.isfinite(value)

    def test_agrees_with_the_pure_python_implementation(self) -> None:
        """Only to `fastmath` tolerance: FMA contraction costs about one ulp."""
        offsets, extents, leaves = self._tables_and_leaves(3)
        out_njit, org_njit = taylor.poly_taylor_branch_metric_apply(
            leaves, offsets, extents,
        )
        out_py, org_py = taylor.poly_taylor_branch_metric_apply.py_func(
            leaves, offsets, extents,
        )
        np.testing.assert_array_equal(org_njit, org_py)
        np.testing.assert_allclose(out_njit, out_py, rtol=1e-14, atol=0)

    def test_tables_depend_only_on_the_stage(self) -> None:
        """D11: no leaf state enters the tables, so one call per stage suffices."""
        first = taylor.metric_branch_tables(
            T_PARENT, T_CHILD, NBINS, DUCY, 3, M_MAX, 500_000,
        )
        second = taylor.metric_branch_tables(
            T_PARENT, T_CHILD, NBINS, DUCY, 3, M_MAX, 500_000,
        )
        np.testing.assert_array_equal(first[0], second[0])
        np.testing.assert_array_equal(first[1], second[1])
        assert len(first[0]) > 1, "single child would make this vacuous"


class TestAggressiveUntouched:
    """Phase 2 must not perturb the shipped path."""

    @pytest.mark.parametrize("poly_order", [2, 3, 4])
    def test_box_branch_counts_unchanged(self, poly_order: int) -> None:
        """Per-axis counts for one doubling, measured before Phase 2 landed."""
        expected = {2: [4, 1], 3: [8, 4, 1], 4: [16, 8, 4, 1]}[poly_order]
        eta = 1.0
        dpar_parent = psr_utils.poly_taylor_step_d_vec(
            poly_order, T_PARENT, NBINS, eta, np.array([F0]), t_ref=0,
        )
        leaves = np.zeros((1, poly_order + 2, 2))
        leaves[0, :-2, 1] = dpar_parent[0]
        leaves[0, -1, 0] = F0
        out, origins = taylor.poly_taylor_branch_batch(
            leaves, (0.0, T_CHILD), NBINS, eta, poly_order, 64,
        )
        counts = [len(np.unique(out[:, j, 0])) for j in range(poly_order)]
        assert counts == expected
        assert len(out) == int(np.prod(expected))
        assert len(origins) == len(out)

    def test_metric_is_not_the_default_strategy(self) -> None:
        fields = {f.name: f for f in PulsarSearchConfig.__attrs_attrs__}
        assert fields["tiling_strategy"].default == "aggressive"
        assert fields["m_max"].default == M_MAX
        assert fields["metric_lattice"].default == "cubic"


def _search_config(strategy: str, poly_order: int = 3) -> PulsarSearchConfig:
    """Build a minimal config; only the fields the branch path reads matter."""
    nsamps, tsamp = 2**14, 64e-6
    tobs = nsamps * tsamp
    limits = ParamLimits.from_upper(
        (F0 - 1, F0 + 1), [6.0, 500.0], (-8.0, 8.0), tobs,
    )
    return PulsarSearchConfig(
        nsamps=nsamps, tsamp=tsamp, nbins=NBINS, eta=1,
        param_limits=limits.limits, bseg_brute=nsamps // 8, bseg_ffa=nsamps // 2,
        prune_poly_order=poly_order, ducy_max=0.5, wtsp=1.2, use_fourier=False,
        tiling_strategy=strategy, branch_max=16,
        # prune_poly_order=5 is how config.py selects the circular-orbit search.
        p_orb_min=tobs if poly_order == 5 else 0.0,
    )


def _dp_functs(cfg: PulsarSearchConfig) -> PrunePolyTaylorDPFuncts:
    """`branch` reads only nbins/eta/poly_order/branch_max/tiling_strategy off self."""
    n = cfg.prune_poly_order
    return PrunePolyTaylorDPFuncts(
        [np.array([0.0])] * (n - 1) + [np.linspace(F0 - 1, F0 + 1, 8)],
        np.ones(n),
        np.array([1] * (n - 1) + [8], dtype=np.int64),
        T_PARENT,
        cfg,
    )


class TestStrategyDispatch:
    """Phase 2 step 6: the `@njit` `branch_func` must route on `tiling_strategy`."""

    def test_metric_strategy_reaches_the_metric_covering(self) -> None:
        cfg = _search_config("metric")
        tables = taylor.metric_branch_tables(
            T_PARENT, T_CHILD, NBINS, cfg.metric_ducy, 3, cfg.m_max,
            cfg.metric_branch_max,
        )
        leaves = np.zeros((2, 5, 2))
        leaves[:, -1, 0] = F0
        out, origins = _dp_functs(cfg).branch(
            leaves, (0.0, T_CHILD), (0.0, T_PARENT), *tables,
        )
        expected, _ = taylor.poly_taylor_branch_metric_apply(leaves, *tables)
        np.testing.assert_array_equal(out, expected)
        assert len(out) == 2 * len(tables[0])
        assert len(origins) == len(out)

    def test_box_strategy_ignores_the_tables(self) -> None:
        """The shipped path must not read them, so junk tables must change nothing."""
        cfg = _search_config("aggressive")
        dpar = psr_utils.poly_taylor_step_d_vec(
            3, T_PARENT, NBINS, 1.0, np.array([F0]), t_ref=0,
        )
        leaves = np.zeros((1, 5, 2))
        leaves[0, :-2, 1] = dpar[0]
        leaves[0, -1, 0] = F0
        funcs = _dp_functs(cfg)
        empty = (np.empty((0, 0)), np.empty(0))
        junk = (RNG.normal(size=(7, 3)), RNG.normal(size=3))
        out_a, org_a = funcs.branch(leaves, (0.0, T_CHILD), (0.0, T_PARENT), *empty)
        out_b, org_b = funcs.branch(leaves, (0.0, T_CHILD), (0.0, T_PARENT), *junk)
        np.testing.assert_array_equal(out_a, out_b)
        np.testing.assert_array_equal(org_a, org_b)
        # ...and still equals the untouched box branch.
        out_ref, _ = taylor.poly_taylor_branch_batch(
            leaves, (0.0, T_CHILD), NBINS, 1.0, 3, 16,
        )
        np.testing.assert_array_equal(out_a, out_ref)

    def test_missing_tables_raise_rather_than_under_cover(self) -> None:
        """An unset table would otherwise branch every parent into nothing."""
        cfg = _search_config("metric")
        leaves = np.zeros((2, 5, 2))
        leaves[:, -1, 0] = F0
        with pytest.raises(ValueError, match="per-stage branch tables"):
            _dp_functs(cfg).branch(
                leaves, (0.0, T_CHILD), (0.0, T_PARENT),
                np.empty((0, 0)), np.empty(0),
            )


class TestPruneWiring:
    """Phase 2 step 6b: the per-stage precompute hook in `prune.py`."""

    @staticmethod
    def _pruning(cfg: PulsarSearchConfig) -> Pruning:
        prn = Pruning.__new__(Pruning)
        prn._stage_table_cache = {}
        prn._empty_stage_tables = (np.empty((0, 0)), np.empty(0), np.empty(0))
        prn._dyp = SimpleNamespace(cfg=cfg)
        prn._logger = logging.getLogger("test_prune_wiring")
        return prn

    def test_tables_are_empty_for_the_box_strategies(self) -> None:
        prn = self._pruning(_search_config("aggressive"))
        offsets, extents, trans = prn._metric_stage_tables(
            (0.0, T_CHILD), (0.0, T_PARENT), (1.0, 2 * T_CHILD),
        )
        assert offsets.shape == (0, 0)
        assert extents.shape == (0,)
        assert trans.shape == (0,)

    def test_tables_are_built_once_per_stage(self) -> None:
        """D11/D14: the enumeration is the expensive half, so it must be memoised."""
        cfg = _search_config("metric")
        prn = self._pruning(cfg)
        cur, prev, nxt = (0.0, T_CHILD), (0.0, T_PARENT), (1.0, 2 * T_CHILD)

        first = prn._metric_stage_tables(cur, prev, nxt)
        second = prn._metric_stage_tables(cur, prev, nxt)
        assert first[0] is second[0], "same stage must hit the cache"
        assert len(prn._stage_table_cache) == 1

        prn._metric_stage_tables((0.0, 2 * T_CHILD), cur, (2.0, 4 * T_CHILD))
        assert len(prn._stage_table_cache) == 2, "a new stage must recompute"
        # delta_t alone is enough to make it a different stage for the transform.
        prn._metric_stage_tables(cur, prev, (5.0, 2 * T_CHILD))
        assert len(prn._stage_table_cache) == 3, "delta_t must be part of the key"

        expected = taylor.metric_branch_tables(
            T_PARENT, T_CHILD, NBINS, cfg.metric_ducy, 3, cfg.m_max,
            cfg.metric_branch_max,
        )
        np.testing.assert_array_equal(first[0], expected[0])
        np.testing.assert_array_equal(first[1], expected[1])
        np.testing.assert_array_equal(
            first[2],
            taylor.metric_transform_extents(
                T_CHILD, 1.0, NBINS, cfg.metric_ducy, 3, cfg.m_max,
            ),
        )

    @pytest.mark.parametrize(
        ("poly_basis", "poly_order"),
        [("chebyshev", 3), ("taylor", 5)],  # poly_order=5 selects the circular basis
    )
    def test_metric_is_rejected_on_the_bases_that_cannot_honour_it(
        self, poly_basis: str, poly_order: int,
    ) -> None:
        """Silently ignoring the strategy would be worse than refusing to run."""
        cfg = _search_config("metric", poly_order)
        prn = Pruning.__new__(Pruning)
        prn._dyp = SimpleNamespace(
            cfg=cfg, fold=np.zeros((2,) * 6), param_arr=[], dparams_actual=None,
        )
        with pytest.raises(ValueError, match="only implemented for the Taylor basis"):
            prn._setup_pruning(poly_basis, use_moving_grid=True)


class TestMetricConfigDefaults:
    """The knobs step 6b needed, and the assumption baked into one of them."""

    def test_metric_ducy_follows_ducy_max_unless_set(self) -> None:
        assert _search_config("metric").metric_ducy == 0.5
        cfg = PulsarSearchConfig(
            nsamps=2**14, tsamp=64e-6, nbins=NBINS, eta=1,
            param_limits=ParamLimits.from_upper(
                (F0 - 1, F0 + 1), [6.0, 500.0], (-8.0, 8.0), 2**14 * 64e-6,
            ).limits,
            bseg_brute=2**11, bseg_ffa=2**13, prune_poly_order=3,
            ducy_max=0.5, wtsp=1.2, use_fourier=False,
            tiling_strategy="metric", branch_max=16, metric_ducy=0.05,
        )
        assert cfg.metric_ducy == 0.05

    def test_metric_branch_max_is_not_the_box_branch_max(self) -> None:
        """`branch_max` is a per-axis width; the metric cap is a total (D14 note)."""
        cfg = _search_config("metric")
        assert cfg.branch_max == 16
        assert cfg.metric_branch_max >= 500_000
        # 16 as a total would truncate even the poly_order=3 covering.
        n_children = len(
            taylor.metric_branch_tables(
                T_PARENT, T_CHILD, NBINS, cfg.metric_ducy, 3, cfg.m_max,
                cfg.metric_branch_max,
            )[0],
        )
        assert n_children > cfg.branch_max


class TestTransformExtents:
    """Phase 2 step 2: the exact re-centred extents (DECISIONS.md D18)."""

    @pytest.mark.parametrize("poly_order", [2, 3, 4])
    def test_matches_the_support_function_of_the_mapped_ellipsoid(
        self, poly_order: int,
    ) -> None:
        """The honest half-width is the ellipsoid's own extent, not a box's shear.

        For axis `j` the exact half-width after `d -> T d` is
        `max{ (T d)_j : d^T g d <= m_max } = sqrt(m_max * (T g^-1 T^T)_jj)`.
        Checked here against a direct maximisation over sampled boundary points, which
        must approach it from below.
        """
        delta_t = 7.3
        extents = taylor.metric_transform_extents(
            T_CHILD, delta_t, NBINS, DUCY, poly_order, M_MAX,
        )
        g = metric.poly_phase_metric(0.0, 0.0, T_CHILD, poly_order, 1.0, NBINS, DUCY)
        t_mat = metric.shift_matrix(delta_t, poly_order)

        chol = metric.cholesky_factor(g)
        raw = RNG.normal(size=(200_000, poly_order))
        unit = raw / np.linalg.norm(raw, axis=1, keepdims=True)
        boundary = np.sqrt(M_MAX) * np.linalg.solve(chol.T, unit.T).T
        sampled = np.abs(boundary @ t_mat.T).max(axis=0)

        assert np.all(sampled <= extents * (1 + 1e-9)), "extents must bound the image"
        np.testing.assert_allclose(sampled, extents, rtol=0.05)

    def test_reduces_to_the_branch_extents_at_zero_shift(self) -> None:
        """T(0) = I, so the transform must not move the child's bounding box."""
        _, extents = taylor.metric_branch_tables(
            T_PARENT, T_CHILD, NBINS, DUCY, 3, M_MAX, 500_000,
        )
        at_zero = taylor.metric_transform_extents(
            T_CHILD, 0.0, NBINS, DUCY, 3, M_MAX,
        )
        np.testing.assert_allclose(at_zero, extents, rtol=1e-12)

    def test_scales_as_one_over_f0(self) -> None:
        """D11's f0**2 law, which is what lets this be a per-stage table."""
        unit = taylor.metric_transform_extents(T_CHILD, 3.1, NBINS, DUCY, 4, M_MAX)
        g = metric.poly_phase_metric(0.0, 0.0, T_CHILD, 4, F0, NBINS, DUCY)
        at_f0 = metric.ellipsoid_axis_extents(
            metric.transform_metric(g, metric.shift_matrix(3.1, 4)), M_MAX,
        )
        np.testing.assert_allclose(at_f0 * F0, unit, rtol=1e-10)


class TestTransformDispatch:
    """`poly_taylor_transform_batch` must route on strategy like `branch` does."""

    @staticmethod
    def _leaves(poly_order: int = 3) -> np.ndarray:
        leaves = np.zeros((3, poly_order + 2, 2))
        leaves[:, :-2, 0] = RNG.normal(scale=1e-3, size=(3, poly_order))
        leaves[:, :-2, 1] = RNG.random((3, poly_order)) * 1e-4
        leaves[:, -2, 0] = np.arange(3.0)
        leaves[:, -1, 0] = F0
        leaves[1, -1, 0] = 2.0 * F0
        return leaves

    def test_metric_takes_column_1_from_the_table_and_values_from_the_shift(
        self,
    ) -> None:
        leaves = self._leaves()
        coord_cur, coord_next = (0.0, T_CHILD), (4.5, T_CHILD)
        extents = taylor.metric_transform_extents(
            T_CHILD, 4.5, NBINS, DUCY, 3, M_MAX,
        )
        out = taylor.poly_taylor_transform_batch(
            leaves, coord_next, coord_cur, "metric", extents,
        )
        # Values: the same shift the box strategies apply.
        expected_vals = transforms.shift_taylor_params(
            np.ascontiguousarray(leaves[:, :-1, 0]), 4.5,
        )
        np.testing.assert_allclose(out[:, :-1, 0], expected_vals, rtol=1e-12)
        # Column 1: the per-stage table, rescaled per leaf, and nothing from the input.
        for i in range(len(leaves)):
            np.testing.assert_allclose(
                out[i, :-2, 1], extents / leaves[i, -1, 0], rtol=1e-12,
            )
        assert np.all(out[:, -2, 1] == 0.0), "d_0 has no defined half-width"
        np.testing.assert_array_equal(out[:, -1], leaves[:, -1])

    def test_box_strategies_ignore_the_table(self) -> None:
        leaves = self._leaves()
        coord_cur, coord_next = (0.0, T_CHILD), (4.5, T_CHILD)
        out_empty = taylor.poly_taylor_transform_batch(
            leaves, coord_next, coord_cur, "aggressive", np.empty(0),
        )
        out_junk = taylor.poly_taylor_transform_batch(
            leaves, coord_next, coord_cur, "aggressive", RNG.normal(size=3),
        )
        np.testing.assert_array_equal(out_empty, out_junk)

    def test_missing_table_raises(self) -> None:
        with pytest.raises(ValueError, match="per-stage transform extents"):
            taylor.poly_taylor_transform_batch(
                self._leaves(), (4.5, T_CHILD), (0.0, T_CHILD), "metric", np.empty(0),
            )


class TestMetricBranchingPattern:
    """Phase 2 step 2, second half: B(s) for the threshold scheme."""

    def test_matches_the_actual_per_level_child_counts(self) -> None:
        cfg = _search_config("metric")
        nsegments, ref_seg = 4, 1
        pattern = taylor.generate_bp_poly_taylor_metric(
            cfg.tseg_ffa, nsegments, ref_seg, NBINS, cfg.metric_ducy,
            cfg.prune_poly_order, cfg.m_max, cfg.metric_branch_max,
            use_moving_grid=True,
        )
        scheme = MiddleOutScheme(nsegments, ref_seg, cfg.tseg_ffa, stride=1)
        expected = []
        for lvl in range(1, nsegments):
            _, t_cur = scheme.get_current_coord(lvl, moving_grid=True)
            _, t_prev = scheme.get_previous_coord(lvl, moving_grid=True)
            offsets, _ = taylor.metric_branch_tables(
                t_prev, t_cur, NBINS, cfg.metric_ducy, cfg.prune_poly_order,
                cfg.m_max, cfg.metric_branch_max,
            )
            expected.append(float(len(offsets)))
        assert len(pattern) == nsegments - 1
        np.testing.assert_array_equal(pattern, expected)
        assert np.all(pattern >= 1)

    def test_config_dispatches_to_it(self) -> None:
        cfg = _search_config("metric")
        pattern = cfg.generate_branching_pattern(
            kind="poly_taylor_moving", ref_seg=1,
        )
        nsegments = int(np.ceil(cfg.nsamps / cfg.bseg_ffa))
        assert len(pattern) == nsegments - 1
        np.testing.assert_array_equal(
            pattern,
            taylor.generate_bp_poly_taylor_metric(
                cfg.tseg_ffa, nsegments, 1, cfg.nbins, cfg.metric_ducy,
                cfg.prune_poly_order, cfg.m_max, cfg.metric_branch_max,
                use_moving_grid=True,
            ),
        )

    def test_non_taylor_kinds_are_refused(self) -> None:
        cfg = _search_config("metric")
        with pytest.raises(ValueError, match="only implemented for the Taylor basis"):
            cfg.generate_branching_pattern(kind="poly_chebyshev_moving", ref_seg=1)

    def test_the_approximate_pattern_is_refused(self) -> None:
        """A worst-case per-axis factor is meaningless for a covering."""
        cfg = _search_config("metric")
        with pytest.raises(ValueError, match="no approximate branching pattern"):
            cfg.generate_branching_pattern_approx(kind="poly_taylor_moving", ref_seg=1)
