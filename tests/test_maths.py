import numpy as np
import pytest
from numpy import polynomial
from scipy import special, stats

from pyloki.utils import maths, transforms

# ``maths.norm_isf_func`` linearly interpolates a table of ``norm.isf(exp(-x))``
# sampled every ``maths.minus_logsf_res`` (= 0.1). The tabulated function is
# steepest as x -> 0 (it diverges to -inf at x = 0), so the interpolation error
# peaks in the first table cells and decays quickly outwards. The tolerances below
# are the accuracy the function actually delivers, measured on a 10**6-point grid
# over [0, 10]:
#     [0.1, 0.2)  max |err| 2.7e-2   <- exceeds the usual decimal=2 bound (1.5e-2)
#     [0.2, 0.3)  max |err| 1.0e-2
#     [0.3, 1.0)  max |err| 5.4e-3
#     [1.0, 10]   max |err| 7.5e-4
# Below 0.1 the function is not usable at all, see
# ``TestMaths.test_norm_isf_func_below_first_table_node``.
# The points are deliberately off the 0.1-spaced table nodes, where the error is
# identically zero and the test would prove nothing.
NORM_ISF_POINTS = [
    (0.125, 3e-2),
    (0.145, 3e-2),
    (0.175, 3e-2),
    (0.230, 1.5e-2),
    (0.270, 1.5e-2),
    (0.350, 6e-3),
    (0.470, 6e-3),
    (0.660, 6e-3),
    (1.050, 1e-3),
    (1.630, 1e-3),
    (2.440, 1e-3),
    (3.870, 1e-3),
    (5.550, 1e-3),
    (7.310, 1e-3),
    (9.420, 1e-3),
]

# Same idea for ``maths.chi_sq_minus_logsf_func``, whose table is sampled every
# ``maths.chi_sq_res`` (= 0.5): points off the table nodes, spanning [0, 10].
CHI_SQ_POINTS = [0.25, 0.75, 1.6, 3.3, 6.2, 9.8]


class TestMaths:
    @pytest.mark.parametrize(
        ("n", "k"),
        [(5, 2), (6, 7), (2, 3), (20, 12), (5, 0), (5, 5), (0, 0)],
    )
    def test_nbinom(self, n: int, k: int) -> None:
        np.testing.assert_almost_equal(maths.nbinom(n, k), special.binom(n, k))

    @pytest.mark.parametrize("n", [0, 1, 5, 20, np.array([0, 1, 5, 20])])
    def test_fact(self, n: int | np.ndarray) -> None:
        np.testing.assert_almost_equal(maths.fact(n), special.factorial(n))

    @pytest.mark.parametrize(("minus_logsf", "atol"), NORM_ISF_POINTS)
    def test_norm_isf_func(self, minus_logsf: float, atol: float) -> None:
        expected = stats.norm.isf(np.exp(-minus_logsf))
        np.testing.assert_allclose(
            maths.norm_isf_func(minus_logsf),
            expected,
            rtol=0,
            atol=atol,
        )

    @pytest.mark.parametrize("minus_logsf", [0.0, 1e-9, 0.05, 0.099])
    def test_norm_isf_func_below_first_table_node(self, minus_logsf: float) -> None:
        """Pin the current behaviour below the first table node.

        ``gen_norm_isf_table`` starts at x = 0, where ``norm.isf(exp(-0)) == -inf``,
        so every interpolation inside the first cell [0, ``minus_logsf_res``) mixes
        in that infinity and comes back non-finite. The true values there are
        perfectly finite (e.g. ~-1.66 at x = 0.05), so this is a defect of the
        table's lower edge rather than an interpolation limit; callers must keep
        ``minus_logsf >= maths.minus_logsf_res``. If the table gains a usable first
        entry, this test is the one to delete.
        """
        result = maths.norm_isf_func(minus_logsf)
        assert not np.isfinite(result), f"expected non-finite, got {result}"

    def test_norm_isf_func_first_cell_accuracy_limit(self) -> None:
        """Characterise the error spike in the first *finite* table cell.

        With ``minus_logsf_res = 0.1`` the table cannot follow the curvature of
        ``norm.isf(exp(-x))`` just above zero: the error peaks at ~2.7e-2 near
        x = 0.145, i.e. the function is not accurate to two decimals there. This
        is why ``NORM_ISF_POINTS`` carries a looser tolerance below x = 0.2.
        """
        x = np.linspace(maths.minus_logsf_res, 2 * maths.minus_logsf_res, 201)
        got = np.array([maths.norm_isf_func(val) for val in x])
        err = np.abs(got - stats.norm.isf(np.exp(-x)))
        np.testing.assert_array_less(
            1.5e-2,
            err.max(),
            err_msg="first cell now meets decimal=2; tighten NORM_ISF_POINTS",
        )
        np.testing.assert_array_less(err.max(), 3e-2)

    @pytest.mark.parametrize("df", [2, 3, 5, 10, 32])
    @pytest.mark.parametrize("chi_sq", CHI_SQ_POINTS)
    def test_chi_sq_minus_logsf_func(self, chi_sq: float, df: int) -> None:
        expected = -stats.chi2.logsf(chi_sq, df)
        np.testing.assert_allclose(
            maths.chi_sq_minus_logsf_func(chi_sq, df),
            expected,
            rtol=0,
            atol=1.5e-2,
        )

    def test_chi_sq_minus_logsf_func_df1_accuracy_limit(self) -> None:
        """Characterise the one df where the chi2 table misses two decimals.

        ``-chi2.logsf(x, 1)`` has an infinite slope at x = 0, and with
        ``chi_sq_res = 0.5`` the first cell is far too coarse for it: the error
        reaches ~1.4e-1 near x = 0.12 and stays above 1.5e-2 out to x ~ 0.47.
        Every other df tested is within 1.1e-2 over [0, 10] (worst: df = 3).
        """
        x = np.linspace(0, maths.chi_sq_res, 201)
        got = np.array([maths.chi_sq_minus_logsf_func(val, 1) for val in x])
        err = np.abs(got + stats.chi2.logsf(x, 1))
        np.testing.assert_array_less(
            1.5e-2,
            err.max(),
            err_msg="df=1 now meets decimal=2; fold df=1 into CHI_SQ_POINTS",
        )
        np.testing.assert_array_less(err.max(), 1.5e-1)

    @pytest.mark.parametrize(("order_max", "n_derivs"), [(3, 1), (5, 3), (10, 10)])
    def test_gen_chebyshev_polys_table(self, order_max: int, n_derivs: int) -> None:
        expected = maths.gen_chebyshev_polys_table_np(order_max, n_derivs)
        np.testing.assert_equal(
            expected.shape,
            (n_derivs + 1, order_max + 1, order_max + 1),
        )
        np.testing.assert_almost_equal(
            maths.gen_chebyshev_polys_table(order_max, n_derivs),
            expected,
            decimal=2,
        )


class TestChebyshevTransform:
    def test_connection_coefficients_s(self) -> None:
        # S_{0,0} = 1
        np.testing.assert_almost_equal(
            maths.compute_connection_coefficient_s(0, 0),
            1.0,
        )
        # S_{2,0} = 1/2, S_{2,2} = 1/2
        np.testing.assert_almost_equal(
            maths.compute_connection_coefficient_s(2, 0),
            0.5,
        )
        np.testing.assert_almost_equal(
            maths.compute_connection_coefficient_s(2, 2),
            0.5,
        )
        # S_{3,1} = 3/4, S_{3,3} = 1/4
        np.testing.assert_almost_equal(
            maths.compute_connection_coefficient_s(3, 1),
            0.75,
        )
        np.testing.assert_almost_equal(
            maths.compute_connection_coefficient_s(3, 3),
            0.25,
        )

    def test_connection_coefficients_r(self) -> None:
        # R_{2,0} = 1, R_{2,2} = 2
        np.testing.assert_almost_equal(
            maths.compute_connection_coefficient_r(2, 0),
            -1.0,
        )
        np.testing.assert_almost_equal(
            maths.compute_connection_coefficient_r(2, 2),
            2.0,
        )

    def test_taylor_to_cheby_manual(self) -> None:
        # Snap parameters
        d_vec = np.array([0.5, 2.3, 1500.0, 1e6, 1e4])
        t_s = 4.2
        # Expected coefficients
        alpha_4 = d_vec[0] * t_s**4 / (8 * maths.fact(4))
        alpha_3 = d_vec[1] * t_s**3 / (4 * maths.fact(3))
        alpha_2 = 0.5 * (
            (d_vec[2] * t_s**2 / maths.fact(2)) + (d_vec[0] * t_s**4 / (maths.fact(4)))
        )
        alpha_1 = d_vec[3] * t_s + (0.75 * d_vec[1] * t_s**3 / maths.fact(3))
        alpha_0 = (
            d_vec[4]
            + d_vec[2] * t_s**2 / (2 * maths.fact(2))
            + 3 * d_vec[0] * t_s**4 / (8 * maths.fact(4))
        )
        alpha_expected = np.array([alpha_4, alpha_3, alpha_2, alpha_1, alpha_0])
        alpha = transforms.taylor_to_cheby(d_vec, t_s)
        np.testing.assert_almost_equal(alpha, alpha_expected, decimal=12)

    def test_cheby_to_taylor_manual(self) -> None:
        alpha_vec = np.array([0.5, 2.3, 1500.0, 1e6, 1e4])
        t_s = 4.2
        # Expected coefficients
        d_4 = 192 * alpha_vec[0] / t_s**4
        d_3 = 24 * alpha_vec[1] / t_s**3
        d_2 = 4 * (alpha_vec[2] - 4 * alpha_vec[0]) / t_s**2
        d_1 = 1 * (alpha_vec[3] - 3 * alpha_vec[1]) / t_s
        d_0 = alpha_vec[4] - alpha_vec[2] + alpha_vec[0]
        d_expected = np.array([d_4, d_3, d_2, d_1, d_0])
        d = transforms.cheby_to_taylor(alpha_vec, t_s)
        np.testing.assert_almost_equal(d, d_expected, decimal=12)

    @pytest.mark.parametrize("k_max", [2, 4, 6])
    def test_roundtrip_identity(self, k_max: int) -> None:
        rng = np.random.default_rng(42)
        d_vec = rng.random(k_max + 1)
        t_s = 1.5
        alpha = transforms.taylor_to_cheby(d_vec, t_s)
        d_reconstructed = transforms.cheby_to_taylor(alpha, t_s)
        np.testing.assert_almost_equal(d_vec, d_reconstructed, decimal=12)

    def test_polynomial_evaluation(self) -> None:
        d_vec = np.array([0.5, 2.3, 1500.0, 1e6, 1e4])
        k_max = len(d_vec) - 1
        t_c, t_s = 4.2, 2.6
        alpha_vec = transforms.taylor_to_cheby(d_vec, t_s)
        t_test = np.linspace(t_c - t_s, t_c + t_s, 11)
        x = (t_test - t_c) / t_s
        k_range = np.arange(k_max + 1)
        c_power = d_vec[::-1] / maths.fact(k_range)
        taylor_poly = polynomial.Polynomial(c_power)
        val_taylor = taylor_poly(t_test - t_c)
        cheby_poly = polynomial.Chebyshev(alpha_vec[::-1], domain=[-1, 1])
        val_cheby = cheby_poly(x)
        np.testing.assert_almost_equal(val_taylor, val_cheby, decimal=8)
