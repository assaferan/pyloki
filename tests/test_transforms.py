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

    def test_parity_limitation_when_one_order_dropped(self) -> None:
        # Dropping a single order perturbs only the retained coefficients sharing its
        # parity; opposite-parity coefficients are bit-identical to naive truncation.
        t_s = 3.7
        for k_max in range(2, 7):
            n_params = k_max + 1
            d_vec = self.rng.random(n_params) * 10.0
            econ = transforms.economize_taylor_params(d_vec, t_s, n_keep=k_max)
            for order in range(k_max):
                idx = n_params - 1 - order
                if (order - k_max) % 2:  # opposite parity to the dropped order
                    assert econ[idx] == pytest.approx(d_vec[idx], abs=1e-12), (
                        f"k_max={k_max}, order={order} should be untouched"
                    )
                else:
                    assert econ[idx] != d_vec[idx], (
                        f"k_max={k_max}, order={order} should have moved"
                    )

    def test_poly_order_3_leaves_accel_and_delay_untouched(self) -> None:
        # The resolve step's poly_order=3 case: d_vec is [jerk, accel, vel, delay] and
        # n_keep=3 drops jerk alone (order 3, odd). Accel (order 2) and delay (order 0)
        # are even, so economization cannot move them -- the resolved accel grid cell is
        # exactly what naive truncation picks. Only vel (order 1) changes.
        d_vec = np.array([6.0, 500.0, 12.0, 3.0])  # jerk, accel, vel, delay
        econ = transforms.economize_taylor_params(d_vec, 4.2, n_keep=3)

        assert econ[1] == d_vec[1]  # accel: bit-identical
        assert econ[3] == d_vec[3]  # delay: bit-identical
        assert econ[2] != d_vec[2]  # vel: moved
