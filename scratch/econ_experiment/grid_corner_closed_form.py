"""Closed forms for the grid corner error, and what that error is actually worth.

Successor to `grid_guarantee.py`. Two separate questions:

A. CLOSED FORM for the corner of a cell -- the sup-norm phase distance between a
   cell centre and the cell corner, every branched coefficient displaced by half a
   grid step. `grid_guarantee.py` tabulated this; here it is in closed form and
   verified against the shipped step functions for poly_order 2..8.

       Taylor, use_cheby=True   : (2**(k_max - 1) - 1/2) * eta/N_b
       Taylor, use_cheby=False  :  (k_max / 2)           * eta/N_b
       Chebyshev                :  (k_max / 2)           * eta/N_b

   Per branched order j (1-based ORDER, not array index): Taylor contributes
   2**(j-2) * eta/N_b, Chebyshev contributes exactly 1/2 * eta/N_b. The geometric
   factor is the `2**k` coarsening in `psr_utils.poly_taylor_step_f`, where `k` is
   the 0-based array index, i.e. `2**(j-1)` in order terms. Remove it and the two
   bases give the *same* corner error, so the corner gap is a property of the
   coarsening factor, not of the basis.

B. COVERING RADIUS -- max over points in a cell of the distance to the NEAREST
   grid point (not to the containing cell's centre), optionally modulo an additive
   constant, since the constant term (d_0 / alpha_0) is never branched and a
   constant phase offset only rotates the folded profile. This is the quantity a
   grid search actually achieves. It tells a very different story from (A).

Run:  PYTHONPATH=src python scratch/econ_experiment/grid_corner_closed_form.py
"""

from __future__ import annotations

import itertools
import math

import numpy as np
from numpy.polynomial import chebyshev as npcheb

from pyloki.utils import psr_utils, transforms
from pyloki.utils.misc import C_VAL

NBINS, ETA = 64, 1.0
TOL = ETA / NBINS
F0 = 1.0 / 0.007  # 142.86 Hz
TOBS = 16.78


def cheb_corner(k_max: int, nbins: int, eta: float, f0: float, t_s: float) -> float:
    """Corner error via the codebase's own cheby_to_taylor, in units of eta/nbins."""
    da = psr_utils.poly_cheb_step_vec(k_max, nbins, eta, np.array([f0]))[0]
    alpha_dev = np.zeros(k_max + 1)
    alpha_dev[:k_max] = da / 2.0  # [alpha_kmax .. alpha_1]; alpha_0 is not branched
    d_dev = transforms.cheby_to_taylor(alpha_dev, t_s)[::-1]  # -> [d_0, d_1, ...]
    tt = np.linspace(-t_s, t_s, 200001)
    delay = sum(dk * tt**k / math.factorial(k) for k, dk in enumerate(d_dev))
    return f0 * np.max(np.abs(delay)) / C_VAL / (eta / nbins)


def taylor_corner(k_max: int, tobs: float, nbins: int, eta: float, f0: float,
                  use_cheby: bool) -> float:
    dd = psr_utils.poly_taylor_step_d_vec(
        k_max, tobs, nbins, eta, np.array([f0]), 0, use_cheby
    )[0][::-1]
    ph = sum((d / 2.0) * tobs ** (j + 1) / math.factorial(j + 1)
             for j, d in enumerate(dd))
    return f0 * ph / C_VAL / (eta / nbins)


def _design_taylor(k_max: int, tobs: float) -> np.ndarray:
    t = np.linspace(0.0, tobs, 1501)
    return np.stack([F0 * t**j / math.factorial(j) / C_VAL / TOL
                     for j in range(1, k_max + 1)])


def _design_cheb(k_max: int) -> np.ndarray:
    x = np.linspace(-1.0, 1.0, 1501)
    return np.stack([F0 * npcheb.chebval(x, np.eye(k_max + 1)[j]) / C_VAL / TOL
                     for j in range(1, k_max + 1)])


def covering_radius(step, design, nsamp, mrange, modconst, seed=0):
    rng = np.random.default_rng(seed)
    k = len(step)
    lattice = np.array(list(itertools.product(mrange, repeat=k))) @ (
        np.diag(step) @ design
    )
    pts = (rng.random((nsamp, k)) - 0.5) * step
    worst = 0.0
    for q in pts @ design:
        d = q[None, :] - lattice
        val = 0.5 * (d.max(1) - d.min(1)) if modconst else np.abs(d).max(1)
        worst = max(worst, float(val.min()))
    return worst


def main() -> None:
    print("All numbers are multiples of the nominal tolerance eta/nbins.")
    print("k_max == poly_order == number of BRANCHED coefficients (orders 1..k_max);")
    print("the constant term is never branched, in either basis.\n")

    print("== A1. corner error vs closed form, poly_order 2..8 ==")
    print(f"{'k_max':>6} {'taylor 2^k':>12} {'2^(k-1)-1/2':>12} "
          f"{'taylor raw':>11} {'chebyshev':>10} {'k_max/2':>9}")
    for k_max in range(2, 9):
        print(f"{k_max:>6} {taylor_corner(k_max, TOBS, NBINS, ETA, F0, True):>12.6f} "
              f"{2.0**(k_max - 1) - 0.5:>12.6f} "
              f"{taylor_corner(k_max, TOBS, NBINS, ETA, F0, False):>11.6f} "
              f"{cheb_corner(k_max, NBINS, ETA, F0, TOBS / 2):>10.6f} "
              f"{k_max / 2:>9.6f}")

    print("\n== A2. Chebyshev corner is k_max/2 independent of "
          "nbins, eta, f_max, t_s ==")
    worst, where = 0.0, None
    for k_max in range(2, 9):
        for nbins in (16, 64, 1024):
            for eta in (0.1, 1.0, 2.0):
                for f0 in (1.0, 142.857, 1000.0):
                    for t_s in (0.5, 8.39, 512.0):
                        got = cheb_corner(k_max, nbins, eta, f0, t_s)
                        rel = abs(got - k_max / 2) / (k_max / 2)
                        if rel > worst:
                            worst, where = rel, (k_max, nbins, eta, f0, t_s, got)
    print(f"  worst relative deviation from k_max/2: {worst:.3e} at "
          f"(k_max, nbins, eta, f0, t_s, value) = {where}")

    print("\n== A3. per-branched-order contribution ==")
    print(f"{'order j':>8} {'chebyshev':>10} {'taylor 2^k':>11} {'2^(j-2)':>9}")
    k_max = 8
    da = psr_utils.poly_cheb_step_vec(k_max, NBINS, ETA, np.array([F0]))[0][::-1]
    dd = psr_utils.poly_taylor_step_d_vec(
        k_max, TOBS, NBINS, ETA, np.array([F0]), 0, True
    )[0][::-1]
    x = np.linspace(-1, 1, 20001)
    for j in range(1, k_max + 1):
        c = np.zeros(k_max + 1)
        c[j] = da[j - 1] / 2.0
        cheb_axis = F0 * np.max(np.abs(npcheb.chebval(x, c))) / C_VAL / TOL
        t_axis = (F0 * (dd[j - 1] / 2.0) * TOBS**j / math.factorial(j) / C_VAL / TOL)
        print(f"{j:>8} {cheb_axis:>10.6f} {t_axis:>11.6f} {2.0**(j - 2):>9.4f}")

    print("\n== A4. cell volume in a common coordinate system (Taylor d-coords) ==")
    print(f"{'k_max':>6} {'cheb / taylor(2^k)':>19} {'2^(k(k+1)/2)':>13}")
    for k_max in (2, 3, 4, 5):
        t_s = TOBS / 2
        dd_v = psr_utils.poly_taylor_step_d_vec(
            k_max, TOBS, NBINS, ETA, np.array([F0]), 0, True
        )[0]
        da_v = psr_utils.poly_cheb_step_vec(k_max, NBINS, ETA, np.array([F0]))[0]
        jac = np.zeros((k_max, k_max))
        for j in range(k_max):
            e = np.zeros(k_max + 1)
            e[j] = 1.0
            jac[:, j] = transforms.cheby_to_taylor(e, t_s)[:k_max]
        ratio = abs(np.linalg.det(jac @ np.diag(da_v))) / float(np.prod(dd_v))
        print(f"{k_max:>6} {ratio:>19.6g} {2.0 ** (k_max * (k_max + 1) / 2):>13.6g}")

    print("\n== B. covering radius (max over the cell of the distance to the "
          "NEAREST grid point) ==")
    for modconst in (False, True):
        tag = "modulo an additive constant" if modconst else "plain sup-norm"
        print(f"\n  -- {tag} --")
        print(f"  {'k_max':>6} {'taylor(2^k)':>12} "
              f"{'taylor raw':>11} {'chebyshev':>10}")
        for k_max in (2, 3, 4, 5):
            mr = range(-3, 4) if k_max <= 3 else range(-2, 3)
            nsamp = 3000 if k_max <= 4 else 1200
            dd_c = psr_utils.poly_taylor_step_d_vec(
                k_max, TOBS, NBINS, ETA, np.array([F0]), 0, True)[0][::-1]
            dd_r = psr_utils.poly_taylor_step_d_vec(
                k_max, TOBS, NBINS, ETA, np.array([F0]), 0, False)[0][::-1]
            da_c = psr_utils.poly_cheb_step_vec(
                k_max, NBINS, ETA, np.array([F0]))[0][::-1]
            m_t, m_c = _design_taylor(k_max, TOBS), _design_cheb(k_max)
            print(f"  {k_max:>6} "
                  f"{covering_radius(dd_c, m_t, nsamp, mr, modconst):>12.4f} "
                  f"{covering_radius(dd_r, m_t, nsamp, mr, modconst):>11.4f} "
                  f"{covering_radius(da_c, m_c, nsamp, mr, modconst):>10.4f}")


if __name__ == "__main__":
    main()
