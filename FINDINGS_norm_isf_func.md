# `norm_isf_func` table-edge defects — not fixed here, written up for the maintainer

Found while making `tests/test_maths.py` deterministic (commit `00e04e5`). The test
change pins the current behaviour; **no library code was touched**. These are the
reasons to consider changing it. All against `upstream/main` at `18d04b3`.

`norm_isf_func` (`src/pyloki/utils/maths.py:80-89`) linearly interpolates
`norm_isf_table`, built by `gen_norm_isf_table` at resolution `minus_logsf_res = 0.1`
over 4000 entries.

## 1. Negative input returns ~+28 sigma — wrong in the worst direction

    norm_isf_func(-0.5)  =  28.116
    norm_isf_func(-1.0)  =  28.098
    norm_isf_func(-10.0) =  27.776

`pos = minus_logsf / 0.1` goes negative, `int(pos)` truncates toward zero, and the
negative index wraps to the tail of the table. The argument is out of domain
(`sf = exp(-x) > 1`), so the honest answer is "very negative"; returning the table
maximum means a non-detection is reported as a maximal detection.

**Reachable.** `src/pyloki/detection/scoring.py:654` computes
`norm_isf_func(max(x_single, x_double))` where each `x` is
`chi_sq_minus_logsf_func(...) - log2(n_filters)`, routinely negative for noise.
Measured: 200 pure-Gaussian-noise profiles through
`MatchedFilter(widths=[1,2,4,8], nbins=64).compute_dot_double` gave **191/200 above
20 sigma, median 28.01, min -inf**.

**Mitigating:** `compute_dot_double` and `harmonic_summing_score_func` are not called
anywhere in the repo or the examples, so the live search path does not reach this.
They are public API, so a user can.

## 2. The whole first cell `[0, 0.1)` is non-finite

`norm_isf_table[0] = norm.isf(exp(-0)) = -inf`, and every interpolation in that cell
mixes it in: 10000/10000 grid points below 0.1 come back non-finite, although the true
values there are finite (about -1.66 at x = 0.05). A lower-edge defect, not an
interpolation limit.

## 3. Interpolation error peaks just above the edge

`norm.isf(exp(-x))` is steepest near zero, so linear interpolation is worst in the first
finite cell: max error **0.02687** at x = 0.14471, exceeding `decimal=2` (which is
`1.5e-2`, not `1e-2`) over `[0.1139, 0.1802]`. Decays fast: 1.0e-2 on [0.2,0.3),
5.4e-3 from 0.3, 7.5e-4 from 1.0, 2.3e-4 from 2.0.

Combined, (2) and (3) made the old `uniform(0, 10)` test fail about **1.6%** of runs,
of which ~1.0% was the non-finite cell.

## 4. Same shape in `chi_sq_minus_logsf_func` at `df = 1`

`chi_sq_res = 0.5` against a minus-logsf with infinite slope at 0: max error **0.1397**
near x = 0.12, above `1.5e-2` out to x ~ 0.47 — it would fail ~4.7% of a uniform [0,10]
draw. `df = 3` is next at 1.12e-2. Other `df` values tested are comfortable.

## Suggested fix (all in `maths.py`; (1) and (2) are the ones with teeth)

1. Guard `minus_logsf <= 0` and return a large negative sentinel rather than indexing —
   kills the negative-index wrap and the `-inf` at entry 0.
2. Start `gen_norm_isf_table` at `minus_logsf_res`, or overwrite entry 0 with
   `norm.isf(1 - eps)`, removing the non-finite cell.
3. Optional: a 10x finer grid below x = 1 costs ~100 entries and cuts the interpolation
   error about 100x (error scales as h^2), to ~3e-4.

**Precedent.** Commit "Stop trials_scheme returning -inf before the search has branched"
fixed the same `norm.isf(1) = -inf` trap in `thresholding.py`. This is that trap one
layer down.

## Also noted: library-level unseeded RNGs (out of scope, no seed argument exists)

- `simulation/pulse.py:220,228,261,338` — reached by `tests/test_prune.py` and the four
  `test_example_*` tests.
- `detection/thresholding.py:740,1046,1094` — reached by `tests/test_thresholding.py`
  and `tests/test_example_thresholds.py`.
- `sensitivity/sim_ffa.py:198`.

A `seed`/`rng` parameter on `PulseSignalConfig.generate` and `DynamicThresholdScheme`
would let the suite pin these. Library change, hence a separate ask.
