"""Deterministic test of the section 11 prediction: how well does each basis' grid
actually cover parameter space?

The recovered-SNR comparison (run_basis.py) is noise-dominated: "score of the closest
surviving candidate" swings by +-3.5 between realizations because which candidates clear
thresholding depends on the noise. That metric cannot resolve a grid property.

This measures the grid property directly and without noise. For a true signal, find the
nearest grid point in each basis using the code's own step rules, then compute the
sup-norm phase error over the segment between the true polynomial and that grid point.
Averaged over random grid offsets, since where the signal falls inside a cell is
arbitrary. No pruning, no folding, no noise.
"""

import math

import numpy as np

from pyloki.utils import psr_utils, transforms
from pyloki.utils.misc import C_VAL

rng = np.random.default_rng(0)
NBINS, ETA = 64, 1.0
TOL = ETA / NBINS
F0 = np.array([1.0 / 0.007])
N_OFFSETS = 20000


def sup_phase_err(d_true, d_approx, t_s, n_grid=513):
    """max |phase(d_true) - phase(d_approx)| in cycles over [-t_s, t_s]."""
    t = np.linspace(-t_s, t_s, n_grid)
    n = len(d_true)
    diff = np.zeros_like(t)
    for k in range(n):
        c = (d_true[n - 1 - k] - d_approx[n - 1 - k]) / math.factorial(k)
        diff = diff + c * t**k
    return float(np.max(np.abs(F0[0] * diff / C_VAL)))


print(f"nominal tolerance eta/nbins = {TOL:.6g} cycles")
print(f"averaging over {N_OFFSETS} random grid offsets per cell\n")
print(f"{'poly_order':>10} {'t_s [s]':>9} | {'taylor mean':>12} {'taylor p95':>11} "
      f"| {'cheby mean':>11} {'cheby p95':>10} | {'mean ratio':>10}")
print("-" * 92)

for poly_order in (3, 4, 5):
    for t_s in (16.78,):
        n_par = poly_order + 1  # [d_k ... d_1, d_0]

        # Taylor steps, exactly as poly_taylor_branch_batch computes them.
        dd = psr_utils.poly_taylor_step_d_vec(
            poly_order, t_s, NBINS, ETA, F0, t_ref=0, use_cheby=True
        )[0]
        # dd covers [d_kmax ... d_1]; d_0 is never branched on.
        step_taylor = np.concatenate([dd, [0.0]])

        # Chebyshev steps: uniform per coefficient.
        da = psr_utils.poly_cheb_step_vec(poly_order, NBINS, ETA, F0)[0]
        step_cheby = np.concatenate([da, [0.0]])

        t_errs, c_errs = [], []
        for _ in range(N_OFFSETS):
            # A signal sitting at a uniformly random position inside one cell:
            # the offset from the nearest grid point is uniform in [-step/2, step/2].
            off_t = (rng.random(n_par) - 0.5) * step_taylor
            off_c = (rng.random(n_par) - 0.5) * step_cheby

            # Taylor: the offset IS the coefficient error, directly in Taylor space.
            t_errs.append(sup_phase_err(off_t, np.zeros(n_par), t_s))

            # Chebyshev: the offset is in alpha space; convert to Taylor to evaluate
            # the same phase functional.
            d_off = transforms.cheby_to_taylor(off_c, t_s)
            c_errs.append(sup_phase_err(d_off, np.zeros(n_par), t_s))

        t_errs, c_errs = np.array(t_errs), np.array(c_errs)
        print(f"{poly_order:>10} {t_s:>9.2f} | {t_errs.mean():>12.5f} "
              f"{np.percentile(t_errs, 95):>11.5f} | {c_errs.mean():>11.5f} "
              f"{np.percentile(c_errs, 95):>10.5f} | "
              f"{t_errs.mean() / c_errs.mean():>10.2f}")

print("\nsame, expressed in units of the nominal tolerance:")
print(f"{'poly_order':>10} | {'taylor mean':>12} {'cheby mean':>11} | "
      f"{'taylor p95':>11} {'cheby p95':>10}")
print("-" * 64)
for poly_order in (3, 4, 5):
    t_s = 16.78
    n_par = poly_order + 1
    dd = psr_utils.poly_taylor_step_d_vec(poly_order, t_s, NBINS, ETA, F0, 0, True)[0]
    step_taylor = np.concatenate([dd, [0.0]])
    da = psr_utils.poly_cheb_step_vec(poly_order, NBINS, ETA, F0)[0]
    step_cheby = np.concatenate([da, [0.0]])
    t_errs, c_errs = [], []
    for _ in range(N_OFFSETS):
        off_t = (rng.random(n_par) - 0.5) * step_taylor
        off_c = (rng.random(n_par) - 0.5) * step_cheby
        t_errs.append(sup_phase_err(off_t, np.zeros(n_par), t_s))
        c_errs.append(sup_phase_err(transforms.cheby_to_taylor(off_c, t_s),
                                    np.zeros(n_par), t_s))
    t_errs, c_errs = np.array(t_errs), np.array(c_errs)
    print(f"{poly_order:>10} | {t_errs.mean() / TOL:>12.3f} {c_errs.mean() / TOL:>11.3f} | "
          f"{np.percentile(t_errs, 95) / TOL:>11.3f} {np.percentile(c_errs, 95) / TOL:>10.3f}")
