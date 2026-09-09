"""Do the Taylor and Chebyshev grids deliver the same sup-norm phase guarantee?

Both are quoted at the same tolerance eta/nbins. If they don't actually deliver it
equally, then comparing leaf counts at equal eta compares different things.
"""
import math
import numpy as np
from pyloki.utils import psr_utils
from pyloki.utils.misc import C_VAL

nbins, eta = 64, 1.0
tol = eta / nbins                     # nominal tolerance, in phase (cycles)
f0 = np.array([1.0 / 0.007])          # 142.86 Hz

print(f"nominal tolerance eta/nbins = {tol:.6g} cycles\n")
print(f"{'poly_order':>10} {'T [s]':>9} {'taylor 2^k':>12} {'taylor raw':>12} {'chebyshev':>11}")
print("-" * 60)

for poly_order in (2, 3, 4, 5):
    for T in (2.0, 16.78, 67.1):
        # --- Taylor cell: coefficients ordered [d_kmax ... d_1]; worst case is all
        # half-step errors aligned, evaluated at the far end of the span.
        def taylor_err(use_cheby):
            dd = psr_utils.poly_taylor_step_d_vec(
                poly_order, T, nbins, eta, f0, t_ref=0, use_cheby=use_cheby
            )[0]
            # dd[::-1] is ascending order k = 1..poly_order (d_1 .. d_kmax)
            asc = dd[::-1]
            phase = 0.0
            for i, dparam in enumerate(asc):
                k = i + 1
                phase += (dparam / 2.0) * T**k / math.factorial(k)
            return float(f0[0] * phase / C_VAL)

        # --- Chebyshev cell: uniform half-step per coefficient, |T_k| <= 1 so the
        # worst case is the plain sum.
        da = psr_utils.poly_cheb_step_vec(poly_order, nbins, eta, f0)[0]
        cheb = float(f0[0] * np.sum(da / 2.0) / C_VAL)

        print(f"{poly_order:>10} {T:>9.2f} {taylor_err(True):>12.4g} "
              f"{taylor_err(False):>12.4g} {cheb:>11.4g}")

print("\nratio to the nominal tolerance (values > 1 mean the cell UNDER-covers):")
print(f"{'poly_order':>10} {'T [s]':>9} {'taylor 2^k':>12} {'chebyshev':>11}")
print("-" * 46)
for poly_order in (2, 3, 4, 5):
    for T in (16.78,):
        dd = psr_utils.poly_taylor_step_d_vec(poly_order, T, nbins, eta, f0, 0, True)[0]
        asc = dd[::-1]
        ph = sum((d / 2.0) * T ** (i + 1) / math.factorial(i + 1) for i, d in enumerate(asc))
        t_err = f0[0] * ph / C_VAL
        da = psr_utils.poly_cheb_step_vec(poly_order, nbins, eta, f0)[0]
        c_err = f0[0] * np.sum(da / 2.0) / C_VAL
        print(f"{poly_order:>10} {T:>9.2f} {t_err / tol:>12.3f} {c_err / tol:>11.3f}")
