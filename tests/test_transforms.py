import math

import numpy as np
import pytest

from pyloki.utils import transforms


def evaluate_taylor_poly(d_vec: np.ndarray, t: np.ndarray) -> np.ndarray:
    """Evaluate sum_k d_vec[-(k+1)] / k! * t**k for d_vec ordered [d_kmax,...,d_0]."""
    n_params = len(d_vec)
    result = np.zeros_like(t, dtype=np.float64)
    for k in range(n_params):
        d_k = d_vec[n_params - 1 - k]
        result = result + d_k / math.factorial(k) * t**k
    return result


class TestEconomizeTaylorParams:
    rng = np.random.default_rng(42)

    @pytest.mark.parametrize("n_keep", [1, 2, 3, 5])
    def test_no_op_when_n_keep_geq_n_params(self, n_keep: int) -> None:
        d_vec = self.rng.random(n_keep)
        t_s = 1.5
        result = transforms.economize_taylor_params(d_vec, t_s, n_keep)
        np.testing.assert_array_equal(result, d_vec)

    def test_agrees_with_naive_when_no_high_order_content(self) -> None:
        # Only accel, vel, delay populated; jerk/snap exactly zero.
        d_vec = np.array([0.0, 0.0, 1500.0, 2.3, 0.5])
        t_s = 4.2
        econ = transforms.economize_taylor_params(d_vec, t_s, n_keep=3)
        np.testing.assert_allclose(econ[-3:], d_vec[-3:], atol=1e-8)
        np.testing.assert_allclose(econ[:-3], 0.0, atol=1e-8)

    def test_economization_reduces_worst_case_reconstruction_error(self) -> None:
        # Snap-order (k_max=4) vector with non-negligible jerk and snap terms.
        d_vec = np.array([1e6, 2.5e4, 1500.0, 2.3, 0.5])
        t_s = 4.2
        t_samples = np.linspace(-t_s, t_s, 201)

        phi_exact = evaluate_taylor_poly(d_vec, t_samples)
        phi_naive = evaluate_taylor_poly(d_vec[-3:], t_samples)

        d_econ = transforms.economize_taylor_params(d_vec, t_s, n_keep=3)
        phi_econ = evaluate_taylor_poly(d_econ[-3:], t_samples)

        naive_max_err = np.max(np.abs(phi_naive - phi_exact))
        econ_max_err = np.max(np.abs(phi_econ - phi_exact))

        assert econ_max_err < naive_max_err

    def test_roundtrip_identity_full_order(self) -> None:
        # Economizing to the full order should be an exact roundtrip (no info lost).
        d_vec = self.rng.random(5)
        t_s = 1.5
        result = transforms.economize_taylor_params(d_vec, t_s, n_keep=5)
        np.testing.assert_array_equal(result, d_vec)
