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

import attrs
import numpy as np
import pytest
from numba import njit, prange
from scipy.special import gamma

from pyloki.config import ParamLimits, PulsarSearchConfig
from pyloki.core import metric, taylor
from pyloki.dynamic import PrunePolyTaylorDPFuncts
from pyloki.prune import Pruning
from pyloki.utils import psr_utils, transforms
from pyloki.utils.misc import C_VAL
from pyloki.utils.snail import MiddleOutScheme

RNG = np.random.default_rng(20260910)
F0 = 1.0 / 0.007  # 142.857 Hz, the spin period used throughout the plan
NBINS = 64
DUCY = 0.1
M_MAX = 0.2
# One EP stage: the accumulated segment doubles.
T_PARENT = 16.78
T_CHILD = 33.55
# How far ahead of the leaf's epoch the child window is centred (D20). The leaf is
# expanded about the previous centre, so the newly added segment moves the window on.
DELTA_T = T_CHILD - T_PARENT


def metrics_for(poly_order: int, f0: float = F0) -> tuple[np.ndarray, np.ndarray]:
    """Parent- and child-stage metrics for one segment doubling.

    Both are about the leaf's epoch, which is the parent window's centre. The parent
    window is therefore symmetric; the child window has grown and moved on, so it sits
    `DELTA_T` ahead. `T_PARENT` and `T_CHILD` are half-widths, not endpoints (D20).
    """
    g_parent = metric.poly_phase_metric(
        0.0, -T_PARENT, T_PARENT, poly_order, f0, NBINS, DUCY,
    )
    g_child = metric.poly_phase_metric(
        0.0, DELTA_T - T_CHILD, DELTA_T + T_CHILD, poly_order, f0, NBINS, DUCY,
    )
    return g_parent, g_child



def parent_region(poly_order: int, f0: float = 1.0) -> np.ndarray:
    """Parent region as an explicit form (D23): the stage-(s-1) mismatch ellipsoid.

    This is the region the old two-metric API assumed implicitly. It stays useful for
    testing the covering itself, where a *self-consistent* parent is exactly right --
    D22 was about the search feeding in a region no leaf occupied, not about the
    covering being wrong on the region it is given.
    """
    g_parent, _ = metrics_for(poly_order, f0)
    return metric.region_from_metric(g_parent, M_MAX)


def branch_tables(
    poly_order: int, max_children: int = 500_000,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """`metric_branch_tables` on the self-consistent parent region."""
    return taylor.metric_branch_tables(
        parent_region(poly_order), T_CHILD, DELTA_T, NBINS, DUCY,
        poly_order, M_MAX, max_children,
    )

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
    """`poly_taylor_branch_metric_batch` must be drop-in for the box version.

    Since D23 the wrapper reads the parent region off the leaves' column 1, exactly as
    the search does for a seed, so these leaves carry a real cell there.
    """

    def _leaves(self, n_leaves: int, poly_order: int) -> np.ndarray:
        leaves = np.zeros((n_leaves, poly_order + 2, 2))
        leaves[:, :-2, 0] = RNG.normal(scale=1e-3, size=(n_leaves, poly_order))
        leaves[:, -2, 0] = RNG.normal(size=n_leaves)
        leaves[:, -1, 0] = F0
        leaves[:, -1, 1] = 0.0
        # Column 1 is a full span (D24); make the cell big enough that the stage can
        # actually refine it, or the guard fires and the contract tests go vacuous.
        span = 2.0 * metric.region_axis_extents(parent_region(poly_order)) / F0
        leaves[:, :-2, 1] = span
        return leaves

    @staticmethod
    def _branch(
        leaves: np.ndarray, poly_order: int, cap: int = 500_000,
    ) -> tuple[np.ndarray, np.ndarray]:
        return taylor.poly_taylor_branch_metric_batch(
            leaves, (0.0, T_CHILD), (0.0, T_PARENT), (DELTA_T, T_CHILD), NBINS, DUCY,
            poly_order, M_MAX, cap,
        )

    @pytest.mark.parametrize("poly_order", [2, 3])
    def test_shapes_and_origins(self, poly_order: int) -> None:
        n_leaves = 5
        leaves = self._leaves(n_leaves, poly_order)
        out, origins = self._branch(leaves, poly_order)
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
        out, origins = self._branch(leaves, poly_order)
        np.testing.assert_allclose(out[:, -2, 0], leaves[origins, -2, 0])
        np.testing.assert_allclose(out[:, -1, 0], leaves[origins, -1, 0])
        np.testing.assert_allclose(out[:, -1, 1], leaves[origins, -1, 1])

    @pytest.mark.parametrize("poly_order", [2, 3])
    def test_children_cover_their_own_parent(self, poly_order: int) -> None:
        """End-to-end: the branch output covers the region it was actually given.

        Distinct `f0` per leaf, so this also exercises the `1 / f0` rescaling that lets
        one enumeration serve the whole batch. The region is the one the branch derived
        from column 1, not an assumed mismatch ellipsoid -- that substitution is what
        D22 was about.
        """
        leaves = self._leaves(3, poly_order)
        leaves[:, -1, 0] = np.array([F0, 2.0 * F0, 0.5 * F0])
        region_unit = taylor.metric_seed_region(leaves)
        out, origins = self._branch(leaves, poly_order)
        _, g_child = metrics_for(poly_order, f0=1.0)
        for i in range(len(leaves)):
            f0 = leaves[i, -1, 0]
            centre = leaves[i, :-2, 0]
            offsets = out[origins == i, :-2, 0] - centre
            # Region and metric both scale as f0**2, so evaluate at this leaf's f0.
            samples = sample_parent_region(region_unit * f0**2, 1.0, 1500)
            uncovered, mean_cover, worst = cover_stats(
                samples, offsets, g_child * f0**2, M_MAX,
            )
            assert uncovered == 0, f"leaf {i} (f0={f0:.1f}): {uncovered} uncovered"
            assert worst <= 1.0
            assert mean_cover >= 1.0

    def test_dparam_column_is_the_child_region_span(self) -> None:
        """Column 1 is the region's per-axis FULL span (D24), not a half-width."""
        poly_order = 3
        leaves = self._leaves(2, poly_order)
        out, _ = self._branch(leaves, poly_order)
        _, _, region_new = taylor.metric_branch_tables(
            taylor.metric_seed_region(leaves), T_CHILD, DELTA_T, NBINS, DUCY,
            poly_order, M_MAX, 500_000,
        )
        expected = 2.0 * metric.region_axis_extents(region_new) / F0
        np.testing.assert_allclose(out[0, :-2, 1], expected, rtol=1e-10)
        assert np.all(out[:, :-2, 1] > 0.0)

    def test_cap_propagates(self) -> None:
        leaves = self._leaves(2, 3)
        with pytest.raises(ValueError, match="lattice points"):
            self._branch(leaves, 3, cap=5)


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
        offsets, extents, _ = branch_tables(poly_order)
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
        first = branch_tables(3)
        second = branch_tables(3)
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
        offsets, dparams, _ = taylor.metric_branch_tables(
            parent_region(3), T_CHILD, DELTA_T, NBINS, cfg.metric_ducy, 3, cfg.m_max,
            cfg.metric_branch_max,
        )
        tables = (offsets, dparams)
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
        prn._empty_stage_tables = (np.empty((0, 0)), np.empty(0), np.empty(0))
        prn._region_form = parent_region(cfg.prune_poly_order)
        prn._dyp = SimpleNamespace(cfg=cfg)
        prn._logger = logging.getLogger("test_prune_wiring")
        return prn

    def test_tables_are_empty_for_the_box_strategies(self) -> None:
        prn = self._pruning(_search_config("aggressive"))
        offsets, dparams, trans = prn._metric_stage_tables(
            (0.0, T_CHILD), (0.0, T_PARENT), (DELTA_T, T_CHILD),
        )
        assert offsets.shape == (0, 0)
        assert dparams.shape == (0,)
        assert trans.shape == (0,)

    def test_tables_match_the_standalone_functions(self) -> None:
        cfg = _search_config("metric")
        prn = self._pruning(cfg)
        region_in = prn._region_form.copy()
        offsets, dparams, trans = prn._metric_stage_tables(
            (0.0, T_CHILD), (0.0, T_PARENT), (DELTA_T, T_CHILD),
        )
        exp_off, exp_dp, region_new = taylor.metric_branch_tables(
            region_in, T_CHILD, DELTA_T, NBINS, cfg.metric_ducy, 3, cfg.m_max,
            cfg.metric_branch_max,
        )
        exp_region, exp_trans = taylor.metric_transform_region(region_new, DELTA_T, 3)
        np.testing.assert_array_equal(offsets, exp_off)
        np.testing.assert_array_equal(dparams, exp_dp)
        np.testing.assert_array_equal(trans, exp_trans)
        np.testing.assert_array_equal(prn._region_form, exp_region)

    def test_the_region_is_carried_across_levels(self) -> None:
        """D23: the parent region is run state, not a function of the baseline."""
        cfg = _search_config("metric")
        prn = self._pruning(cfg)
        # Compare extents, not the forms: the forms are ~1e-18, so `allclose`'s
        # default atol of 1e-8 would call any two of them equal.
        def extents(form: np.ndarray) -> np.ndarray:
            return metric.region_axis_extents(form)

        before = extents(prn._region_form)
        prn._metric_stage_tables(
            (0.0, T_CHILD), (0.0, T_PARENT), (DELTA_T, T_CHILD),
        )
        after = extents(prn._region_form)
        assert np.any(after / before < 0.99), "the region must shrink as it refines"
        # ...and the next level starts from where the last one left off.
        prn._metric_stage_tables(
            (0.0, 2 * T_CHILD), (0.0, T_CHILD), (DELTA_T, 2 * T_CHILD),
        )
        assert np.any(extents(prn._region_form) / after < 0.99)

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
                parent_region(3), T_CHILD, DELTA_T, NBINS, cfg.metric_ducy, 3,
                cfg.m_max, cfg.metric_branch_max,
            )[0],
        )
        assert n_children > cfg.branch_max


class TestTransformExtents:
    """Phase 2 step 2: carrying a region to the new epoch (D18), on the right window."""

    @pytest.mark.parametrize("poly_order", [2, 3, 4])
    def test_agrees_with_transporting_the_metric(self, poly_order: int) -> None:
        """A region form transforms exactly as a metric does.

        `{d : d^T A d <= 1}` under `d' = T d` is `{d' : d'^T T^-T A T^-1 d' <= 1}`, the
        same map `transform_metric` applies. Pinned against the mismatch-ellipsoid case,
        where the answer is independently known: transporting `g / m_max` must equal
        rebuilding `g` about the new epoch over the same absolute window (Phase 1
        test 2), in the exact configuration the pruning loop uses.
        """
        g_old = metric.poly_phase_metric(
            0.0, DELTA_T - T_CHILD, DELTA_T + T_CHILD, poly_order, 1.0, NBINS, DUCY,
        )
        region_next, span = taylor.metric_transform_region(
            metric.region_from_metric(g_old, M_MAX), DELTA_T, poly_order,
        )
        g_new = metric.poly_phase_metric(
            0.0, -T_CHILD, T_CHILD, poly_order, 1.0, NBINS, DUCY,
        )
        expected = metric.region_from_metric(g_new, M_MAX)
        # A symmetric window kills the odd moments exactly, so `expected` has hard
        # zeros off the diagonal while the transported form carries rounding noise
        # there. Scale the absolute tolerance to the form rather than compare to 0.
        np.testing.assert_allclose(
            region_next, expected, rtol=1e-8, atol=1e-14 * np.abs(expected).max(),
        )
        # And the span is the full width of that region, not a half-width (D24).
        np.testing.assert_allclose(
            span, 2.0 * metric.ellipsoid_axis_extents(g_new, M_MAX), rtol=1e-8,
        )

    @pytest.mark.parametrize("poly_order", [2, 3, 4])
    def test_span_matches_the_support_function(self, poly_order: int) -> None:
        """Half the span is `max{d_j : d^T A d <= 1} = sqrt((A^-1)_jj)`.

        Checked against a direct maximisation over sampled boundary points, which must
        approach it from below.
        """
        region = parent_region(poly_order)
        half = metric.region_axis_extents(region)
        chol = metric.cholesky_factor(region)
        raw = RNG.normal(size=(200_000, poly_order))
        unit = raw / np.linalg.norm(raw, axis=1, keepdims=True)
        boundary = np.linalg.solve(chol.T, unit.T).T
        sampled = np.abs(boundary).max(axis=0)
        assert np.all(sampled <= half * (1 + 1e-9)), "extents must bound the region"
        np.testing.assert_allclose(sampled, half, rtol=0.05)

    def test_uses_the_symmetric_window_not_the_half_width_as_an_endpoint(self) -> None:
        """The D20 regression guard.

        `coord[1]` is a half-width, so the window is `[-T, +T]`. Passing `[0, T]` puts
        the epoch at the window's edge over half its length and misprices the region by
        a large factor -- 4x to 64x, growing with `poly_order`.
        """
        poly_order = 3
        symmetric = metric.ellipsoid_axis_extents(
            metric.poly_phase_metric(
                0.0, -T_CHILD, T_CHILD, poly_order, 1.0, NBINS, DUCY,
            ),
            M_MAX,
        )
        one_sided = metric.ellipsoid_axis_extents(
            metric.poly_phase_metric(
                0.0, 0.0, T_CHILD, poly_order, 1.0, NBINS, DUCY,
            ),
            M_MAX,
        )
        _, span = taylor.metric_transform_region(
            metric.region_from_metric(
                metric.poly_phase_metric(
                    0.0, DELTA_T - T_CHILD, DELTA_T + T_CHILD, poly_order, 1.0,
                    NBINS, DUCY,
                ),
                M_MAX,
            ),
            DELTA_T,
            poly_order,
        )
        np.testing.assert_allclose(span, 2.0 * symmetric, rtol=1e-8)
        assert np.all(one_sided / symmetric > 3.0), "the bug was not benign"

    def test_scales_as_one_over_f0(self) -> None:
        """D11's f0**2 law, which is what lets this be a per-stage table."""
        _, span_unit = taylor.metric_transform_region(parent_region(4, 1.0), 3.1, 4)
        _, span_f0 = taylor.metric_transform_region(parent_region(4, F0), 3.1, 4)
        np.testing.assert_allclose(span_f0 * F0, span_unit, rtol=1e-10)


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
        _, extents = taylor.metric_transform_region(parent_region(3), 4.5, 3)
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

    @staticmethod
    def _seed_region(cfg: PulsarSearchConfig) -> np.ndarray:
        return cfg.metric_seed_region(
            cfg.get_dparams_actual(cfg.niters_ffa, use_cheby_coarsening=False),
            cfg.get_param_arr(
                cfg.get_dparams(cfg.niters_ffa, use_cheby_coarsening=False),
            ),
        )

    def test_replays_the_search_region_recursion(self) -> None:
        """The pattern must be what the branch will really emit, guard included."""
        cfg = _search_config("metric")
        nsegments, ref_seg = 4, 1
        seed = self._seed_region(cfg)
        pattern = taylor.generate_bp_poly_taylor_metric(
            seed, cfg.tseg_ffa, nsegments, ref_seg, NBINS, cfg.metric_ducy,
            cfg.prune_poly_order, cfg.m_max, cfg.metric_branch_max,
            use_moving_grid=True,
        )
        scheme = MiddleOutScheme(nsegments, ref_seg, cfg.tseg_ffa, stride=1)
        region, expected = seed, []
        for lvl in range(1, nsegments):
            ref_cur, t_cur = scheme.get_current_coord(lvl, moving_grid=True)
            ref_next, _ = scheme.get_coord(lvl)
            delta_t = ref_next - ref_cur
            offsets, _, region = taylor.metric_branch_tables(
                region, t_cur, delta_t, NBINS, cfg.metric_ducy,
                cfg.prune_poly_order, cfg.m_max, cfg.metric_branch_max,
            )
            expected.append(float(len(offsets)))
            region, _ = taylor.metric_transform_region(
                region, delta_t, cfg.prune_poly_order,
            )
        assert len(pattern) == nsegments - 1
        np.testing.assert_array_equal(pattern, expected)
        assert np.all(pattern >= 1)

    def test_the_guard_shows_up_as_a_branching_factor_of_one(self) -> None:
        """On a short baseline the metric cannot refine, so B(s) must be exactly 1.

        The pre-D23 code reported hundreds of children here, for a region no leaf
        occupied. Getting 1 is the whole point of the fix.
        """
        cfg = _search_config("metric")
        pattern = taylor.generate_bp_poly_taylor_metric(
            self._seed_region(cfg), cfg.tseg_ffa, 4, 1, NBINS, cfg.metric_ducy,
            cfg.prune_poly_order, cfg.m_max, cfg.metric_branch_max,
            use_moving_grid=True,
        )
        np.testing.assert_array_equal(pattern, np.ones(3))

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
                self._seed_region(cfg), cfg.tseg_ffa, nsegments, 1, cfg.nbins,
                cfg.metric_ducy, cfg.prune_poly_order, cfg.m_max,
                cfg.metric_branch_max, use_moving_grid=True,
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


class TestResolveMismatch:
    """Phase 2 step 4: what `resolve` costs by snapping children onto `G0`."""

    T_HALF = 0.131
    T_ADD = 0.393

    @staticmethod
    def _grid() -> tuple[list[np.ndarray], np.ndarray, np.ndarray]:
        """Build a 2D base grid in (accel, freq), the axes `resolve` snaps to."""
        limits = np.array([[-8.0, 8.0], [-4.0, 4.0], [F0 - 1.0, F0 + 1.0]])
        counts = np.array([1, 8, 16], dtype=np.int64)
        param_arr = [
            psr_utils.range_param(
                limits[i, 0], limits[i, 1], (limits[i, 1] - limits[i, 0]) / counts[i],
            )
            for i in range(3)
        ]
        return param_arr, counts, limits

    def _leaf_resolving_to(self, accel: float, freq: float) -> np.ndarray:
        """Build a leaf whose forward map lands exactly on `(accel, freq)`.

        With `coord_init` at the origin and `coord_add` at `T_ADD`, resolve computes
        `accel_new = d_2 + d_3 * t_a` and `vel_new = d_2 * t_a + d_3 * t_a**2 / 2`.
        Two equations, two unknowns.
        """
        t_a = self.T_ADD
        vel = C_VAL * (1.0 - freq / F0)
        d_3 = 2.0 * (accel * t_a - vel) / t_a**2
        d_2 = accel - d_3 * t_a
        leaf = np.zeros((1, 5, 2))
        leaf[0, 0, 0] = d_3
        leaf[0, 1, 0] = d_2
        leaf[0, -1, 0] = F0
        return leaf

    def _mismatch(self, leaves: np.ndarray) -> np.ndarray:
        param_arr, counts, limits = self._grid()
        return taylor.metric_resolve_mismatch(
            leaves,
            (self.T_ADD, self.T_HALF),
            (0.0, self.T_HALF),
            (0.0, self.T_HALF),
            param_arr,
            counts,
            limits,
            NBINS,
            DUCY,
        )

    def test_children_on_grid_centres_cost_nothing(self) -> None:
        """The floor of the diagnostic: a child already on `G0` loses no mismatch."""
        param_arr, _, _ = self._grid()
        leaves = np.concatenate([
            self._leaf_resolving_to(float(param_arr[-2][i]), float(param_arr[-1][j]))
            for i, j in ((0, 0), (3, 5), (7, 15))
        ])
        m = self._mismatch(leaves)
        assert np.all(m < 1e-18), f"on-grid children should be free, got {m}"

    def test_half_cell_offset_matches_a_hand_computation(self) -> None:
        """Half a frequency cell, priced by hand against the base-segment metric."""
        param_arr, _, _ = self._grid()
        step = float(param_arr[-1][1] - param_arr[-1][0])
        # A hair inside the next cell up, so the nearest centre is half a cell away.
        freq = float(param_arr[-1][5]) + 0.5 * step * (1 - 1e-9)
        leaf = self._leaf_resolving_to(float(param_arr[-2][3]), freq)
        got = self._mismatch(leaf)[0]

        d_vel = C_VAL * 0.5 * step / F0
        g = metric.poly_phase_metric(
            0.0, -self.T_HALF, self.T_HALF, 2, F0, NBINS, DUCY,
        )
        assert got == pytest.approx(g[1, 1] * d_vel**2, rel=1e-5)

    def test_the_metric_region_exceeds_the_whole_search_space(self) -> None:
        """The step 4 finding, pinned (D22).

        `metric_branch_tables` takes the parent's region to be the `m_max` ellipsoid of
        the previous stage's metric, computed from the accumulated baseline alone. Over
        the short baselines of the early pruning levels that ellipsoid is larger than
        the entire search space -- a 0.13 s baseline constrains almost nothing -- so the
        covering scatters children far outside `param_limits`, where `resolve` can only
        clamp them to the grid edge.

        Note what is *not* wrong: on the same baseline the ellipsoid tracks the eta-box
        to within a factor of a few (see the ratio below). The metric is fine. What is
        missing is any reference to the region the parent actually occupies.
        """
        g_parent = metric.poly_phase_metric(
            0.0, -self.T_HALF, self.T_HALF, 3, F0, NBINS, DUCY,
        )
        ellipsoid = metric.ellipsoid_axis_extents(g_parent, M_MAX)
        box = psr_utils.poly_taylor_step_d_vec(
            3, self.T_HALF, NBINS, 1.0, np.array([F0]), t_ref=0,
        )[0] / 2.0
        # The metric is not the problem: it agrees with the box to a factor of a few.
        assert np.all(ellipsoid / box < 100.0)

        # The search space is: half-spans of jerk, accel, and the frequency band.
        search_half = np.array([8.0, 4.0, C_VAL * 1.0 / F0])
        ratio = ellipsoid / search_half
        # A 0.13 s baseline cannot constrain the higher derivatives at all, so the
        # ellipsoid runs past the whole search range on those axes by 6-8 orders.
        assert np.all(ratio[:-1] > 1e6), f"jerk/accel should blow out, got {ratio}"
        # Frequency is the exception: there the ellipsoid is comparable to the band,
        # which is why the damage shows up as jerk/accel scatter, not a frequency sweep.
        assert 0.1 < ratio[-1] < 10.0, f"d_1 should be comparable, got {ratio[-1]}"


class TestRegionIntersection:
    """D25: a child's region is parent AND ellipsoid, never just the ellipsoid."""

    @staticmethod
    def _pair(poly_order: int) -> tuple[np.ndarray, np.ndarray]:
        a_form = parent_region(poly_order)
        _, g_child = metrics_for(poly_order, 1.0)
        return a_form, metric.region_from_metric(g_child, M_MAX)

    @pytest.mark.parametrize("poly_order", [2, 3, 4])
    def test_is_a_sound_outer_bound(self, poly_order: int) -> None:
        """Every point in both regions must land inside the returned one."""
        a_form, b_form = self._pair(poly_order)
        out = metric.region_intersection(a_form, b_form)

        # Sample inside A, then keep those also in B -- a box around A misses the
        # intersection almost surely once poly_order > 2.
        pts = sample_parent_region(a_form, 1.0, 200_000)
        in_b = np.einsum("ni,ij,nj->n", pts, b_form, pts) <= 1.0
        assert in_b.sum() > 100, "test is vacuous; sampling missed the intersection"
        q = np.einsum("ni,ij,nj->n", pts[in_b], out, pts[in_b])
        assert q.max() <= 1.0 + 1e-9, f"not an outer bound: max {q.max()}"

    @pytest.mark.parametrize("poly_order", [2, 3, 4])
    def test_never_grows_the_volume(self, poly_order: int) -> None:
        """The point of D25: the region must not expand as it is carried forward.

        Volume is the right invariant, not the per-axis extent: `t = 0` and `t = 1` are
        both in the family, so maximising `det` can only beat whichever input is
        smaller. An *ellipsoidal* outer bound of a lens can still be wider than the
        narrower input on an individual axis -- unavoidable, and harmless, since the
        covering consumes the form and not its bounding box.

        Taking the child ellipsoid alone grew the region on every unrefined axis, and
        because the ellipsoid then shrinks only slowly the guard never fired again: the
        branch re-tiled at every level, compounding to `prod B(s) = 4.3e33` against the
        box strategy's 27.
        """
        a_form, b_form = self._pair(poly_order)
        out = metric.region_intersection(a_form, b_form)
        # det of the form is inverse to volume, so it must not fall below either input.
        sign, logdet = np.linalg.slogdet(out)
        assert sign > 0
        for other in (a_form, b_form):
            s_o, ld_o = np.linalg.slogdet(other)
            assert sign * logdet >= s_o * ld_o - 1e-9

    def test_reduces_to_the_smaller_region_when_nested(self) -> None:
        """If one region contains the other, the intersection is the inner one."""
        inner = metric.region_form_from_box(np.array([1.0, 2.0, 3.0]))
        outer = metric.region_form_from_box(np.array([10.0, 20.0, 30.0]))
        np.testing.assert_allclose(
            metric.region_intersection(inner, outer), inner, rtol=1e-9,
        )


class TestDeferFactor:
    """D26: `metric_defer_factor` — how far a region may overhang before re-covering."""

    def test_overhang_is_one_for_an_exactly_matching_region(self) -> None:
        _, g_child = metrics_for(3, 1.0)
        region = metric.region_from_metric(g_child, M_MAX)
        assert metric.region_overhang(region, g_child, M_MAX) == pytest.approx(1.0)

    @pytest.mark.parametrize("scale", [0.5, 1.0, 2.0, 5.0])
    def test_overhang_tracks_a_uniform_dilation(self, scale: float) -> None:
        """Scaling a region by `s` must move the overhang to `s`."""
        _, g_child = metrics_for(3, 1.0)
        region = metric.region_from_metric(g_child, M_MAX) / scale**2
        assert metric.region_overhang(region, g_child, M_MAX) == pytest.approx(scale)

    def test_defer_factor_one_is_exact_containment(self) -> None:
        """The default must keep the `m_max` guarantee for every leaf."""
        _, g_child = metrics_for(3, 1.0)
        just_outside = metric.region_from_metric(g_child, M_MAX) / 1.01**2
        assert not metric.region_fits_in_one_child(just_outside, g_child, M_MAX)
        just_inside = metric.region_from_metric(g_child, M_MAX) * 1.01**2
        assert metric.region_fits_in_one_child(just_inside, g_child, M_MAX)

    def test_a_larger_factor_defers_branching(self) -> None:
        _, g_child = metrics_for(3, 1.0)
        region = metric.region_from_metric(g_child, M_MAX) / 2.5**2  # overhang 2.5
        assert not metric.region_fits_in_one_child(region, g_child, M_MAX, 1.0)
        assert not metric.region_fits_in_one_child(region, g_child, M_MAX, 2.0)
        assert metric.region_fits_in_one_child(region, g_child, M_MAX, 3.0)

    def test_deferring_emits_one_child_and_keeps_the_region(self) -> None:
        """Branching deferred means the parent rides on unchanged, as the box does."""
        region = parent_region(3)
        g_child = metric.poly_phase_metric(
            0.0, DELTA_T - T_CHILD, DELTA_T + T_CHILD, 3, 1.0, NBINS, DUCY,
        )
        overhang = metric.region_overhang(region, g_child, M_MAX)
        assert overhang > 1.0, "test is vacuous unless the region really overhangs"

        # Just under the overhang: still branches.
        offsets, _, _ = taylor.metric_branch_tables(
            region, T_CHILD, DELTA_T, NBINS, DUCY, 3, M_MAX, 500_000,
            overhang * 0.99,
        )
        assert len(offsets) > 1

        # Just over it: one child, and the region rides on untouched.
        offsets, _, region_new = taylor.metric_branch_tables(
            region, T_CHILD, DELTA_T, NBINS, DUCY, 3, M_MAX, 500_000,
            overhang * 1.01,
        )
        assert len(offsets) == 1
        np.testing.assert_array_equal(offsets, np.zeros((1, 3)))
        np.testing.assert_array_equal(region_new, region)

    def test_below_one_is_rejected(self) -> None:
        _, g_child = metrics_for(3, 1.0)
        with pytest.raises(ValueError, match="defer_factor"):
            metric.region_fits_in_one_child(parent_region(3), g_child, M_MAX, 0.5)

    def test_config_default_keeps_the_guarantee(self) -> None:
        """Shipping a relaxed default would silently weaken every metric run."""
        fields = {f.name: f for f in PulsarSearchConfig.__attrs_attrs__}
        assert fields["metric_defer_factor"].default == 1.0
        with pytest.raises(ValueError, match="metric_defer_factor"):
            _search_config("metric").__class__(
                **{**attrs.asdict(_search_config("metric")),
                   "metric_defer_factor": 0.5},
            )
