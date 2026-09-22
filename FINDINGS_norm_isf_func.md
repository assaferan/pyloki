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

**Reachable from the public API, not from the live search** — read both halves of
this together; the measurement below means nothing without the second.
`src/pyloki/detection/scoring.py:654` computes
`norm_isf_func(max(x_single, x_double))` where each `x` is
`chi_sq_minus_logsf_func(...) - log2(n_filters)`, routinely negative for noise.
Measured: 200 pure-Gaussian-noise profiles through
`MatchedFilter(widths=[1,2,4,8], nbins=64).compute_dot_double` gave **191/200 above
20 sigma, median 28.01, min -inf**.

**Read the predicate.** That is the score the scoring path returns *for pure noise* — a
false-alarm property of the function. It is **not** a measurement that a full search
emits a 28 sigma candidate end to end, since a search applies thresholds and a pruning
tree on top of this. The number must not be quoted without that sentence attached.

*Independently reproduced* by the `injection-design` session (2026-09-22) from the
docstring rather than from this script, on a different seed: 194/200 above 20 sigma,
finite median **28.01**, finite min -1.05, max 28.12, one non-finite. The count differs
by sampling; the median lands to the decimal. Note the input shape is `(nprof, nbins)`
per the `_compute_snr_double` docstring — passing `(nprof, 2, nbins)` scores the
variance row as a profile and silently doubles `n`.

**The call graph, verified twice** (independently by this branch and by
`injection-design`, 2026-09-22). `norm_isf_func` has exactly two callers,
`scoring.py:654` in `_compute_snr_double` and `scoring.py:670` in
`harmonic_summing_score_func`. Neither `compute_dot_double` nor
`harmonic_summing_score_func` is called anywhere in `src/pyloki`. The live search scores
through `scoring.snr_score_batch_func`, via each dynamic module's `score_func`, which
never reaches `norm_isf_func`.

So **running a search does not hit this**. A user reaches it only by calling the public
scoring API directly. Anyone presenting the 28-sigma measurement must carry this
sentence with it and at the same weight, or the reader will conclude their searches are
emitting false candidates. They are not.

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
