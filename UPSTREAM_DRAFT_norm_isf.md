# DRAFT — issue: `norm_isf_func` out-of-domain behaviour

**NOT POSTED.** Held for assaferan.

---

## `norm_isf_func` returns ~+28 sigma for negative input (public scoring API; not reached by a search)

`utils/maths.py:80-89`. For `minus_logsf < 0` the table index goes negative and Python
wraps it to the tail of the table, so the function returns close to its **maximum**
value where the correct answer is its minimum:

    norm_isf_func(-0.5) =  28.116
    norm_isf_func(-1.0) =  28.098
    norm_isf_func(-5.0) =  27.955
    norm_isf_func( 0.0) =  nan
    table maximum       =  28.130

### Scope first: a search does not reach this

`norm_isf_func` has exactly two callers, `scoring.py:654` in `_compute_snr_double` and
`scoring.py:670` in `harmonic_summing_score_func`. Neither `compute_dot_double` (the
public `MatchedFilter` method that reaches the first) nor `harmonic_summing_score_func`
is called anywhere in `src/pyloki`. The live search scores through
`scoring.snr_score_batch_func` and never reaches `norm_isf_func`.

**So no search output is affected.** This is reachable only by calling the public
scoring API directly. I am reporting it because it is an exported function returning a
maximally wrong answer on input that its own call site routinely produces, and the fix
is small — not because I think it has corrupted anyone's results.

Everything below should be read with that in mind.

### The mechanism

`pos = minus_logsf / minus_logsf_res` gives `int(-0.5/0.1) = -5`, and
`norm_isf_table[-5]` is the fifth entry from the end. The upper end of the range is
guarded explicitly — `if minus_logsf < max_minus_logsf` falls through to an
extrapolation branch — so the lower end looks unconsidered rather than deliberately
wrapped.

Mathematically, `minus_logsf < 0` means `sf = exp(-minus_logsf) > 1`, which is not a
survival probability. The limiting answer is `norm.isf(1) = -inf`, i.e. "as
insignificant as possible". The function returns the opposite end of the scale.

At exactly zero, `norm_isf_table[0] = norm.isf(exp(0)) = -inf`, and the interpolation
multiplies it by a zero weight, giving NaN under `fastmath`.

### Negative input is the ordinary case on that path, not an edge case

`detection/scoring.py:649-654` computes

    x_single = chi_sq_minus_logsf_func(scores_max_single, 1) - lee_penalty_single
    x_double = chi_sq_minus_logsf_func(scores_max_double, 2) - lee_penalty_double
    results[iprof] = norm_isf_func(max(x_single, x_double))

where `lee_penalty_single = np.log2(n_filters)` is a look-elsewhere penalty. Whenever
the raw significance is smaller than the trials penalty — which is what noise does — `x`
is negative. So the wrap is hit on ordinary input to that function, not a pathological
one.

Measured: 200 pure-Gaussian-noise profiles, scored by calling the public
`MatchedFilter(widths=[1,2,4,8], nbins=64).compute_dot_double` API directly rather than
by running a search, gave 191/200 above 20 sigma, median 28.01. Independently
reproduced on a different seed, from the docstring rather than from the script:
194/200 above 20 sigma, median 28.01, max 28.12.

That is a false-alarm property of the scoring function on pure noise. It is **not** a
claim that a search emits a 28 sigma candidate, which the scope section above rules out.

### Suggested fix

A domain check. I have not sent a patch because the right return value is a judgement
about the API rather than about the maths: `-inf`, a clamp of `minus_logsf` at 0, or a
raise are all defensible, and you may have a convention. Happy to prepare whichever you
prefer.

### Secondary, on the same code path — mixed log bases?

**This shares the scope above:** `lee_penalty` appears at exactly four places,
`scoring.py:619, 620, 649, 652`, all inside `_compute_snr_double` — the same function
nothing in `src/pyloki` calls. So this too cannot affect a search, and it is raised as a
consistency question rather than as a live-scoring concern.

While tracing where the negative values come from, the penalty looks like it may be in
the wrong units:

| quantity | base |
|---|---|
| `norm_isf_func(x)`, i.e. `norm.isf(exp(-x))` — verified exact at x = 1, 5, 20 | natural |
| `chi_sq_minus_logsf_func`, i.e. `-chi2.logsf(...)` — verified exact at (10, 1), (25, 2) | natural |
| `lee_penalty = np.log2(...)` | base 2 |

If that is a mismatch rather than a deliberate scaling, the penalty is larger than
intended by `1/ln(2) = 1.4427` when subtracted from a natural-log quantity, which would
also make the negative-input case above more common than it should be. The dimensional
observation is verified; only the intent is open, so I am asking rather than
asserting.
