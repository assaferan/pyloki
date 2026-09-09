"""Phase 1 tests for `pyloki.core.metric` (the five in `metric_PLAN.md`).

Tests 1-4 are the exit criterion; test 5 reports ratios that must be understood rather
than asserted tightly.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from pyloki.core import metric
from pyloki.utils.misc import C_VAL

RNG = np.random.default_rng(20260909)
F0 = 1.0 / 0.007  # 142.857 Hz, the spin period used throughout the plan
NBINS = 64


def time_grid(t_start: float, t_end: float, n: int) -> np.ndarray:
    """Midpoint sampling of `[t_start, t_end]`.

    A plain `linspace` includes both endpoints, so an equal-weight mean of it is a
    rectangle rule with `O(1/n)` endpoint error — which shows up as a ~1e-5 relative
    discrepancy against the metric's exact moments and looks like a bug in the metric.
    The midpoint rule is `O(1/n**2)` and is the right discretisation of a continuous
    time average.
    """
    h = (t_end - t_start) / n
    return t_start + (np.arange(n) + 0.5) * h


def phase_cycles(delta: np.ndarray, t: np.ndarray, t_ref: float, poly_order: int):
    """Phase difference in cycles from a coefficient offset, evaluated directly.

    `Phi = f0 * [(t - t_ref) - d(t)/c]`, so an offset `delta` in the coefficients gives
    `dPhi(t) = -(f0/c) * sum_i delta_i * tau^k_i / k_i!`. Written out independently of
    `metric.py` so the comparison is a real check.
    """
    tau = t - t_ref
    orders = np.arange(poly_order, 0, -1)
    out = np.zeros_like(tau)
    for d_i, k in zip(delta, orders, strict=True):
        out = out + d_i * tau**k / math.factorial(int(k))
    return -(F0 / C_VAL) * out


class TestMetricBasics:
    def test_shape_symmetry_and_positive_definite(self) -> None:
        for poly_order in (1, 2, 3, 4, 5):
            g = metric.poly_phase_metric(50.0, 0.0, 100.0, poly_order, F0, NBINS)
            assert g.shape == (poly_order, poly_order)
            np.testing.assert_allclose(g, g.T, rtol=0, atol=0)
            # positive definite <=> Cholesky succeeds and eigenvalues > 0
            assert np.all(np.linalg.eigvalsh(g) > 0)
            metric.cholesky_factor(g)

    def test_nbins_does_not_affect_g(self) -> None:
        # D7: g is independent of nbins.
        a = metric.poly_phase_metric(50.0, 0.0, 100.0, 4, F0, 16)
        b = metric.poly_phase_metric(50.0, 0.0, 100.0, 4, F0, 1024)
        np.testing.assert_allclose(a, b, rtol=0, atol=0)

    def test_cholesky_whitens(self) -> None:
        g = metric.poly_phase_metric(50.0, 0.0, 100.0, 4, F0, NBINS)
        chol = metric.cholesky_factor(g)
        np.testing.assert_allclose(chol @ chol.T, g, rtol=1e-10)
        # In whitened coordinates w = L^T delta, m == w . w
        delta = RNG.normal(size=4) * np.array([1e-3, 1e-2, 1e-1, 1e0])
        w = chol.T @ delta
        np.testing.assert_allclose(
            metric.mismatch(g, delta), float(w @ w), rtol=1e-10
        )


class TestBruteForceMismatch:
    """Plan test 1 — metric vs direct evaluation of the phase difference."""

    @pytest.mark.parametrize("poly_order", [2, 3, 4, 5])
    def test_metric_equals_phase_variance(self, poly_order: int) -> None:
        t_ref, t_start, t_end = 33.5, 0.0, 67.1
        g = metric.poly_phase_metric(t_ref, t_start, t_end, poly_order, F0, NBINS)

        # Scale each axis by the metric's own extent so delta is a sensible size.
        scale = metric.ellipsoid_axis_extents(g, 1e-4)
        delta = RNG.normal(size=poly_order) * scale

        m_metric = metric.mismatch(g, delta)

        # Direct: 2*pi**2 * Var_t(dPhi) on a dense uniform grid (D6).
        t = time_grid(t_start, t_end, 400_000)
        d_phi = phase_cycles(delta, t, t_ref, poly_order)
        m_direct = 2.0 * np.pi**2 * float(np.var(d_phi))

        # The metric IS the exact covariance of a polynomial, so this is limited only
        # by quadrature, not by any small-delta expansion.
        np.testing.assert_allclose(m_metric, m_direct, rtol=1e-6)

    def test_scales_quadratically_in_delta(self) -> None:
        poly_order = 4
        g = metric.poly_phase_metric(33.5, 0.0, 67.1, poly_order, F0, NBINS)
        delta = RNG.normal(size=poly_order) * metric.ellipsoid_axis_extents(g, 1e-4)
        m1 = metric.mismatch(g, delta)
        m2 = metric.mismatch(g, 2.0 * delta)
        np.testing.assert_allclose(m2, 4.0 * m1, rtol=1e-12)

    def test_batched_delta(self) -> None:
        g = metric.poly_phase_metric(33.5, 0.0, 67.1, 3, F0, NBINS)
        batch = RNG.normal(size=(7, 3)) * metric.ellipsoid_axis_extents(g, 1e-4)
        out = metric.mismatch(g, batch)
        assert out.shape == (7,)
        for i in range(7):
            np.testing.assert_allclose(out[i], metric.mismatch(g, batch[i]), rtol=1e-12)


class TestEpochShiftCovariance:
    """Plan test 2 — `transform_metric` vs recomputing about the new epoch."""

    @pytest.mark.parametrize("poly_order", [2, 3, 4, 5])
    @pytest.mark.parametrize("delta_t", [-17.0, -3.5, 2.0, 11.25])
    def test_transform_matches_recompute(self, poly_order: int, delta_t: float) -> None:
        t_ref, t_start, t_end = 33.5, 0.0, 67.1
        g_old = metric.poly_phase_metric(t_ref, t_start, t_end, poly_order, F0, NBINS)

        t_mat = metric.shift_matrix(delta_t, poly_order)
        g_shifted = metric.transform_metric(g_old, t_mat)

        # Same absolute interval, expansion moved to t_ref + delta_t.
        g_direct = metric.poly_phase_metric(
            t_ref + delta_t, t_start, t_end, poly_order, F0, NBINS
        )
        np.testing.assert_allclose(g_shifted, g_direct, rtol=1e-8)

    def test_shift_matrix_matches_transforms_module(self) -> None:
        """The metric's `T` must be the same matrix the codebase already shifts with."""
        from pyloki.utils import transforms

        poly_order, delta_t = 4, 7.25
        # Full vector including d_0, shifted by the production routine.
        d_full = RNG.normal(size=poly_order + 1)
        shifted_full = transforms.shift_taylor_params(d_full, delta_t)

        # Our block should reproduce it on the branchable axes (all but the last entry),
        # which is only true because the full matrix is lower triangular with d_0 last.
        t_block = metric.shift_matrix(delta_t, poly_order)
        np.testing.assert_allclose(
            t_block @ d_full[:poly_order], shifted_full[:poly_order], rtol=1e-10
        )

    def test_mismatch_is_invariant_under_shift(self) -> None:
        poly_order, delta_t = 4, -9.0
        t_ref, t_start, t_end = 33.5, 0.0, 67.1
        g_old = metric.poly_phase_metric(t_ref, t_start, t_end, poly_order, F0, NBINS)
        t_mat = metric.shift_matrix(delta_t, poly_order)
        g_new = metric.transform_metric(g_old, t_mat)

        delta = RNG.normal(size=poly_order) * metric.ellipsoid_axis_extents(g_old, 1e-4)
        np.testing.assert_allclose(
            metric.mismatch(g_old, delta),
            metric.mismatch(g_new, t_mat @ delta),
            rtol=1e-8,
        )


class TestEllipsoidInvariance:
    """Plan test 3 — the property that dissolves the tiling dilemma."""

    @pytest.mark.parametrize("poly_order", [2, 3, 4])
    def test_boundary_maps_to_boundary(self, poly_order: int) -> None:
        t_ref, t_start, t_end = 33.5, 0.0, 67.1
        m_max, delta_t = 1e-3, 12.5
        g = metric.poly_phase_metric(t_ref, t_start, t_end, poly_order, F0, NBINS)
        chol = metric.cholesky_factor(g)
        t_mat = metric.shift_matrix(delta_t, poly_order)
        g_new = metric.transform_metric(g, t_mat)

        # Sample the boundary exactly: unit sphere in whitened coords, un-whitened.
        raw = RNG.normal(size=(200, poly_order))
        unit = raw / np.linalg.norm(raw, axis=1, keepdims=True)
        boundary = np.sqrt(m_max) * np.linalg.solve(chol.T, unit.T).T

        np.testing.assert_allclose(
            metric.mismatch(g, boundary), np.full(200, m_max), rtol=1e-8
        )
        mapped = boundary @ t_mat.T
        np.testing.assert_allclose(
            metric.mismatch(g_new, mapped), np.full(200, m_max), rtol=1e-8
        )

    def test_ellipsoid_axis_extents_touch_the_boundary(self) -> None:
        g = metric.poly_phase_metric(33.5, 0.0, 67.1, 4, F0, NBINS)
        m_max = 1e-3
        ext = metric.ellipsoid_axis_extents(g, m_max)
        # The bounding-box corner is outside, but the point where the ellipsoid touches
        # the plane delta_i = ext_i is exactly on the boundary.
        g_inv = np.linalg.inv(g)
        for i in range(4):
            touch = m_max * g_inv[:, i] / ext[i]
            np.testing.assert_allclose(metric.mismatch(g, touch), m_max, rtol=1e-8)


class TestEmpiricalSnrLoss:
    """Plan test 4 — does `m` actually predict the S/N loss?"""

    @pytest.mark.parametrize("m_target", [1e-4, 1e-3, 1e-2])
    def test_single_harmonic_amplitude_loss(self, m_target: float) -> None:
        """`A/A_0 = |<exp(2*pi*i*dPhi)>|` should be `1 - m` to second order (D6)."""
        poly_order = 4
        t_ref, t_start, t_end = 33.5, 0.0, 67.1
        g = metric.poly_phase_metric(t_ref, t_start, t_end, poly_order, F0, NBINS)

        raw = RNG.normal(size=poly_order)
        direction = raw / np.linalg.norm(raw)
        # Scale the direction to land exactly on mismatch m_target.
        unit_m = metric.mismatch(g, direction)
        delta = direction * math.sqrt(m_target / unit_m)
        m = metric.mismatch(g, delta)
        np.testing.assert_allclose(m, m_target, rtol=1e-10)

        t = time_grid(t_start, t_end, 200_000)
        d_phi = phase_cycles(delta, t, t_ref, poly_order)
        amp_ratio = abs(np.mean(np.exp(2j * np.pi * d_phi)))

        # 1 - m is the second-order prediction; the residual is O(m**2).
        predicted = 1.0 - m
        assert abs(amp_ratio - predicted) < 12.0 * m_target**2, (
            f"m={m:.3e} predicted={predicted:.8f} measured={amp_ratio:.8f}"
        )

    def test_power_convention_is_twice_amplitude(self) -> None:
        """Makes the plan's factor-of-2 ambiguity explicit rather than hiding it (D6).

        `poly_phase_metric`'s docstring says fractional S/N loss (amplitude), while plan
        test 4 compares `(S/N ratio)**2` to `1 - m` (power). They differ by exactly 2.
        """
        poly_order = 4
        t_ref, t_start, t_end = 33.5, 0.0, 67.1
        g = metric.poly_phase_metric(t_ref, t_start, t_end, poly_order, F0, NBINS)
        m_target = 1e-4
        raw = RNG.normal(size=poly_order)
        direction = raw / np.linalg.norm(raw)
        delta = direction * math.sqrt(m_target / metric.mismatch(g, direction))

        t = time_grid(t_start, t_end, 200_000)
        d_phi = phase_cycles(delta, t, t_ref, poly_order)
        amp_ratio = abs(np.mean(np.exp(2j * np.pi * d_phi)))

        m_amp = 1.0 - amp_ratio
        m_pow = 1.0 - amp_ratio**2
        np.testing.assert_allclose(m_pow / m_amp, 2.0, rtol=1e-3)
        # And our g is the amplitude one.
        np.testing.assert_allclose(m_amp, metric.mismatch(g, delta), rtol=2e-3)

    def test_duty_cycle_dependence_is_reported(self, capsys) -> None:
        """Boxcar-scored folded profile vs the single-harmonic prediction (O4).

        Narrow pulses carry power at higher harmonics, so they lose more S/N than the
        fundamental-only estimate. This measures the effective harmonic weighting rather
        than assuming it; the numbers feed O4 and DECISIONS.md.
        """
        from pyloki.detection.thresholding import generate_folded_profile

        poly_order = 4
        t_ref, t_start, t_end = 33.5, 0.0, 67.1
        g = metric.poly_phase_metric(t_ref, t_start, t_end, poly_order, F0, NBINS)
        m_target = 4e-3
        raw = RNG.normal(size=poly_order)
        direction = raw / np.linalg.norm(raw)
        delta = direction * math.sqrt(m_target / metric.mismatch(g, direction))

        t = time_grid(t_start, t_end, 200_000)
        d_phi = phase_cycles(delta, t, t_ref, poly_order)

        rows = []
        for ducy in (0.5, 0.2, 0.1, 0.05):
            profile = np.asarray(generate_folded_profile(nbins=NBINS, ducy=ducy))
            # Folding with the wrong parameters smears the profile by the distribution
            # of d_phi over the observation: circular convolution with that histogram.
            hist, _ = np.histogram(
                np.mod(d_phi, 1.0), bins=NBINS, range=(0.0, 1.0), density=False
            )
            kernel = hist / hist.sum()
            smeared = np.real(
                np.fft.ifft(np.fft.fft(profile) * np.fft.fft(kernel))
            )
            # Boxcar-matched amplitude, the same statistic the search maximises.
            width = max(1, int(round(ducy * NBINS)))
            box = np.ones(width) / math.sqrt(width)
            best = max(
                float(np.dot(np.roll(smeared, -s)[:width], box)) for s in range(NBINS)
            )
            best_true = max(
                float(np.dot(np.roll(profile, -s)[:width], box)) for s in range(NBINS)
            )
            rows.append((ducy, best / best_true, (1.0 - best / best_true) / m_target))

        with capsys.disabled():
            print(f"\n  boxcar S/N loss vs metric prediction (m_target={m_target:g}):")
            print(f"  {'ducy':>6} {'A/A0':>10} {'loss / m':>10}")
            for ducy, ratio, factor in rows:
                print(f"  {ducy:>6.2f} {ratio:>10.6f} {factor:>10.3f}")

        # Direction, not magnitude: narrower pulses must not lose *less*.
        factors = [f for _, _, f in rows]
        assert factors[-1] >= factors[0] - 0.05, (
            f"narrow pulses should lose at least as much as broad ones: {factors}"
        )
        assert all(r <= 1.0 + 1e-9 for _, r, _ in rows)


class TestSanityAgainstCurrentSpacing:
    """Plan test 5 — ratios against the existing box. Reported, loosely asserted."""

    @pytest.mark.parametrize("poly_order", [2, 3, 4, 5])
    def test_extents_within_orders_of_magnitude(self, poly_order, capsys) -> None:
        eta, t_start, t_end = 1.0, 0.0, 67.108864
        t_ref = 0.0  # poly_taylor_step_d_vec is called with t_ref=0 in branch
        bridge_raw = metric.m_max_from_eta(
            eta, NBINS, poly_order, t_ref, t_start, t_end, F0, use_cheby=False
        )
        bridge_cheby = metric.m_max_from_eta(
            eta, NBINS, poly_order, t_ref, t_start, t_end, F0, use_cheby=True
        )
        g = metric.poly_phase_metric(t_ref, t_start, t_end, poly_order, F0, NBINS)

        with capsys.disabled():
            print(f"\n  poly_order={poly_order}  (D5: metric does NOT inherit 2**k)")
            print(f"    m_max un-coarsened box: axis_tight={bridge_raw.axis_tight:.4g}"
                  f"  axis_loose={bridge_raw.axis_loose:.4g}"
                  f"  volume={bridge_raw.volume:.4g}")
            print(f"    m_max coarsened   box: axis_tight={bridge_cheby.axis_tight:.4g}"
                  f"  axis_loose={bridge_cheby.axis_loose:.4g}"
                  f"  volume={bridge_cheby.volume:.4g}")
            pairs = (("un-coarsened", bridge_raw), ("coarsened", bridge_cheby))
            for name, bridge in pairs:
                for label in ("axis_tight", "volume"):
                    ext = metric.ellipsoid_axis_extents(g, getattr(bridge, label))
                    ratio = ext / bridge.box_half_widths
                    print(f"    {name:>12} / {label:<10} extent:box per axis = "
                          f"{np.array2string(ratio, precision=3)}")

        # A units or ordering bug would show up as many orders of magnitude, not a
        # factor of a few. The tight bridge is inscribed by construction.
        ext_tight = metric.ellipsoid_axis_extents(g, bridge_raw.axis_tight)
        ratio_tight = ext_tight / bridge_raw.box_half_widths
        assert np.all(ratio_tight <= 1.0 + 1e-9)
        assert np.all(ratio_tight > 1e-3), f"suspiciously small: {ratio_tight}"

    def test_coarsening_ratio_is_two_to_the_k(self) -> None:
        """The two bridges must differ exactly by the documented 2**k per axis (D5)."""
        eta, poly_order = 1.0, 4
        raw = metric.m_max_from_eta(
            eta, NBINS, poly_order, 0.0, 0.0, 67.108864, F0, use_cheby=False
        )
        chb = metric.m_max_from_eta(
            eta, NBINS, poly_order, 0.0, 0.0, 67.108864, F0, use_cheby=True
        )
        ratio = chb.box_half_widths / raw.box_half_widths
        # Axis i carries order k = poly_order - i, and the coarsening is 2**(k-1)
        # because poly_taylor_step_f applies 2**k over k = 0.. for [f_k..f_0] reversed.
        np.testing.assert_allclose(ratio, 2.0 ** np.arange(poly_order - 1, -1, -1))
