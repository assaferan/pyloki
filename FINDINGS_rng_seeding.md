# The library never seeded anything — eight unseeded `default_rng()` sites

Companion to `FINDINGS_norm_isf_func.md`. That file records defects found *underneath*
the `test_maths.py` flake; this one records the wider defect the branch name does not
name, and the change made for it.

**Status: prepared locally, nothing posted upstream.** No issue and no PR against
`pravirkr/pyloki`, per assaferan's standing hold.

## The defect

`np.random.default_rng()` was called with no argument in eight places across three
modules, and no public entry point accepted a `seed` or an `rng`:

| site | what it is |
|---|---|
| `detection/thresholding.py:740` | `DynamicThresholdScheme.__init__` → `self.rng`, used at nine call sites in the class |
| `detection/thresholding.py:1046` | `determine_scheme()` |
| `detection/thresholding.py:1094` | `evaluate_scheme()` |
| `simulation/pulse.py:220` | `generate_simple()` — the injected noise |
| `simulation/pulse.py:228` | `generate_noise()` |
| `simulation/pulse.py:261` | `generate_old()` |
| `simulation/pulse.py:338` | `generate()` — the current path |
| `sensitivity/sim_ffa.py:198` | `TestFFASensitivity.__init__` |

Line numbers are against `main` at `18d04b3`.

Two consequences. A run could not be repeated: same inputs, same configuration,
different output. And `np.random.seed` did not help, because `default_rng` ignores the
legacy global seed — which is what made this easy to miss, since the usual reflex
appears to work and changes nothing.

Seeding closes this at seven of the eight sites. It does **not** close
`DynamicThresholdScheme.run()`, for a reason that is about parallelism rather than
seeding; that is measured and explained under **Limitation** below, and it is a real
qualification on everything else in this file.

It is a **reproducibility** defect, not a correctness one. The search itself is
deterministic downstream of the time series; the non-determinism is entirely in the
*inputs* — the noise realisation and the threshold ladder.

## Why it had to be fixed in the library, not in the tests

`tests/test_example_ep_jerk.py` contains no RNG of its own. It calls
`PulseSignalConfig.generate()` and `thresholding.determine_scheme()`, and both drew from
generators constructed inside the library. No amount of seeding in the test file could
reach them. That is the difference from the `test_maths.py` flake, where the draw was in
the test and could be removed there.

## The change

An optional `seed` on each public entry point, accepting `int`, `np.random.Generator`,
or `None`. **`None` remains the default and still draws fresh entropy**, so no existing
caller sees different numbers; this is a pure addition.

For `PulseSignalConfig` the seed is a config field and the generator is **stored on the
config**, not re-created per call. That matters: re-seeding inside each `generate*` call
would make a loop over one seeded config draw the identical realisation every iteration.
Storing it means a seed fixes the *sequence* — successive calls differ, and two runs of
the same script agree. `get_updated()` carries the seed through, and the stored
generator is excluded from the constructor round-trip and from equality (Generators
compare by identity, which would otherwise make two identically-seeded configs unequal).

`DynamicThresholdScheme` needed only its constructor changed: everything downstream
already took `rng` as an explicit parameter (`thresholding.py:50, 393, 445, 544, 610`).

The example tests now pin `SEED = 42`.

## Tests

`tests/test_rng_seeding.py`, 23 tests, ~4 s. For every site: same seed reproduces,
different seeds differ, and **unseeded still varies** — that last one pins that the fix
did not quietly make the library deterministic for callers who never asked.

`test_no_bare_default_rng_in_library` walks the AST of every module under `src/` and
fails on any `default_rng()` called with no arguments, so a newly added unseeded
generator fails the suite instead of going unnoticed. The list is enumerated from the
source, not typed into the test, so it cannot drift.

Verified separately that the two RNG-derived inputs of the `ep_jerk` pipeline — the
`2**22`-sample time series and the 64-stage ladder — are **bit-identical across three
fresh processes** at `seed=42` (sha256 `9c075133b6450c2a` and `83dbe45a9f57f833`).

## The flake, measured

The brief asked for a rate with its denominator, not the inherited "1 in 4". Predicate:
**one run of `tests/test_example_ep_jerk.py` in which at least one of its six tests
fails**, 60 independent unseeded runs, one warm numba cache throughout so the
denominator does not move.

| | failures | rate | 95% CI (Wilson) |
|---|---|---|---|
| unseeded, `ACCEL_TOL = 1.0` | **5 / 60** | 8.3% | [3.6%, 18.1%] |
| 28 seeded realisations, accel check only | 2 / 28 | 7.1% | [2.0%, 22.6%] |

The two rows do not share a predicate — the first is per file run, the second scores
only the acceleration check — but since the acceleration check was the only test that
ever failed, they are measuring the same event in practice. Both are **per-run** rates. They are not comparable to the per-test analytic rate
for `test_norm_isf_func` (0.664%, `metric-gridding` D118); if the two ever appear in one
table they need their denominators written next to them.

All five failures were `test_recovers_injected_acceleration`; no other test in the file
ever failed. The inherited "1 in 4" is not contradicted — `1/4` carries a 95% CI of
[4.6%, 69.9%], which contains 8.3% and nearly everything else. It was never a rate.

## The cause was not the RNG

Seeding made the failure *nameable*: a failing realisation could finally be re-run and
examined, which was impossible before. What it showed is that the RNG was not the defect.

The recovered acceleration error is **bimodal, not continuous**. The best candidate
either lands in the correct acceleration cell or in an adjacent one:

| | error | seeds |
|---|---|---|
| correct cell | 0.0059 – 0.0082 | 26 of 28 |
| one cell off | 8.2787, 8.2915 | 2 of 28 (seeds 7 and 11) |

Nothing falls in between, and `8.285` is exactly `daccel` — the uncertainty the search
itself reports for that parameter.

The test carried `ACCEL_TOL = 1.0`, commented *"measured error <= 0.0082, i.e. >100x
margin"*. That reasoning assumes a continuous error distribution. The real distribution
has an atom at one grid step, so a tolerance of 1.0 — **8x tighter than the uncertainty
the search quotes** — silently required the peak to land in exactly the right cell. The
five-realisation stability sample recorded in the module docstring happened to draw five
exact-cell hits and missed the mode entirely.

The check now runs against the reported `daccel`, which keeps it non-vacuous: the
existing guard already bounds `daccel` under 10% of the injected acceleration. Seeds 7
and 11, both previously failing, now pass end-to-end through pytest.

**This is the grid, not a defect in the search.** `metric-gridding` measured the same
structure from the other side (D48): the nearest template to a signal sits ~0.5 cells
away once you minimise over neighbouring centres, so landing one cell off is an ordinary
outcome for a cell-based search — there is no continuum because there are no templates
between cell centres. 26 in-cell against 2 one-cell-off is a healthy ratio for that
geometry. So the tolerance was not merely tight, it was the wrong *kind* of quantity: an
absolute figure in accel units, checked against an error whose support is
`{~0} ∪ {multiples of daccel}`. Any value below one grid step demands an exact-cell hit
however it is chosen, and any two values between one and two steps are equivalent.
Sizing against `daccel` is the only formulation that tracks the configuration when
`prune_poly_order`, `tobs`, `eta` or the segment count move the grid — the same reason
`metric-gridding` uses the metric mismatch `m <= 1.0` as its recovery predicate rather
than a distance in parameter units.

`JERK_TOL = 0.75` against `djerk = 0.556`, and `FREQ_TOL = 2e-3` against
`dfreq = 9.2e-5`, were already sized above one grid step and so survive a one-cell miss.
Neither failed in any of the 60 unseeded runs, which is the tell — though 60 runs bound
those rates only loosely (0/60 is consistent with anything under ~6%).

## What this does and does not fix

- It **does** let every generator in the library be steered from its caller, and it
  makes `determine_scheme`, `evaluate_scheme` and all four `PulseSignalConfig.generate*`
  paths reproduce exactly, on 14 threads as well as on one.
- It does **not** make `DynamicThresholdScheme.run()` reproducible under parallel
  execution — see the limitation section below. A seed there fixes the stream but not
  the order threads consume it in. This is the one site where seeding is necessary but
  not sufficient, and it is the site the original draft called the most important one.
- It **does** fix the `ep_jerk` flake, and by the tolerance rather than by the pin — see
  the unseeded re-measurement below, which holds the new tolerance and lets the noise
  vary.
- It does **not** claim any published result is wrong. Ratcheting noise into a ladder
  changes which run you get, not whether the search works.
- It does **not** make the search more accurate. A one-cell acceleration miss still
  happens on roughly 1 realisation in 12; the test now accepts it because the search
  never claimed better, not because it stopped happening.
- It does **not** claim seeding was needed to make paired comparisons possible. Those
  were already possible by persisting the time series, which is `injection-design`'s own
  §6.1 finding. Pairability and reproducibility are different goals; a seed is now an
  additional route to the first and the only route to the second.

## Re-measured after the fix

Same predicate, same method, a fresh tree holding the new tolerance with the noise still
**unseeded** — so this measures the tolerance fix alone, not the pin.

| | failures | rate | 95% CI (Wilson) |
|---|---|---|---|
| before: unseeded, `ACCEL_TOL = 1.0` | 5 / 60 | 8.3% | [3.6%, 18.1%] |
| after: unseeded, tolerance vs `daccel` | **0 / 60** | 0% | [0.0%, 6.0%] |

**Do not read that as a significant difference on its own.** Fisher exact on 5/60 against
0/60 gives p = 0.057, which does not clear 0.05. Sixty runs cannot separate 8.3% from 0%.

The evidence is the mechanism, not the p-value. Across all 88 realisations measured here
the largest acceleration error was 8.2915, and the new tolerance is
`1.5 x daccel = 12.43`. Every realisation observed passes by construction, and a failure
now requires a **two**-cell miss, which was never observed. The 0/60 run is a
consistency check on that reasoning, not independent proof of it.

Both measurements used one warm numba cache per tree throughout, so no denominator moved
mid-measurement.

## Limitation: a seed is not sufficient at `DynamicThresholdScheme.run()`

Found by the reproducer after the fix was written, and it corrects a claim this file
made in an earlier revision.

`run_stage_legacy` is `@njit(cache=True, parallel=True)` and draws from the shared
generator **inside a `prange`** (`thresholding.py:597`, loop at `:616`, `rng` consumed at
`:639`). A seed fixes the stream; it does not fix the order threads consume it in. So a
seeded `DynamicThresholdScheme.run()` still produces different ladders:

| seed | numba threads | distinct ladders | max per-stage spread |
|---|---|---|---|
| none | 14 | 5 / 10 | 0.798 |
| none, after `np.random.seed(42)` | 14 | 5 / 10 | 0.798 |
| **42** | **14** | **2 / 5** | **0.080** |
| 42 | 1 | 1 / 3 | 0.000 |

Seeding removes most of the spread (0.798 → 0.080) but not all of it, and the residue is
thread scheduling: forced single-threaded, the same seed reproduces exactly. That
localises the remaining non-determinism precisely and rules out the seed as its source.

**Every other seeded path reproduces under parallelism.** `determine_scheme`,
`evaluate_scheme`, and all four `PulseSignalConfig.generate*` methods were each checked
6x on 14 threads and gave bit-identical output. The limitation is this one site.

Closing it needs per-iteration generators inside the kernel — deriving a child generator
from `(seed, ibeam_cur)` so the result is independent of which thread runs which
iteration — which is a change to a numba kernel rather than a signature, and is **not
done here**.

### How this got past the first round of testing

The original regression test asserted that `scheme.rng` produced the same draws for the
same seed. That is a **proxy**: it tests that the seed reached the generator, not that
the ladder reproduces. It passed while the thing it was supposed to guarantee was false.

`test_run_reproduces_single_threaded` now exercises `run()` itself. Getting it
non-vacuous took a second correction: comparing `states["threshold"]` also passes
trivially, because thresholds come off a fixed `np.linspace` and are identical even
unseeded. The field that carries the noise is `success_h0`. Checked both ways — two
unseeded single-threaded runs differ in it, two seeded runs do not.

The parallel case is deliberately **not** asserted anywhere. A test asserting that two
runs *differ* would itself be flaky, which is the defect this branch exists to remove.

## Note on the site count

The library had eight bare `default_rng()` calls; the fixed tree has five seeded ones.
The four in `pulse.py` collapsed into the single construction in
`__attrs_post_init__`, since the generator is now built once per config instead of once
per `generate*` call. Nothing was dropped — `rng_reproducibility.py` enumerates the
sites from the source, so the count is checkable rather than asserted.
