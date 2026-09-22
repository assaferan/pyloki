# DRAFT — PR description: seedable RNGs

**NOT POSTED.** Held for assaferan. Branch: `assaferan:fix-flaky-norm-isf`.

---

## Make every RNG in the library seedable

### The problem

`np.random.default_rng()` is called with no argument in eight places, and no public
entry point accepts a `seed` or an `rng`:

    simulation/pulse.py        220, 228, 261, 338
    detection/thresholding.py  740, 1046, 1094
    sensitivity/sim_ffa.py     198

(line numbers against `main` at `18d04b3`)

So a run cannot be repeated: same inputs, same configuration, different output. And
`np.random.seed` does not help, because `default_rng` ignores the legacy global seed —
which is what makes this easy to miss, since the usual reflex appears to work and
changes nothing.

The site that matters most is `DynamicThresholdScheme`, because a threshold ladder is
an **input** to a search, not an output of one. Two searches a user believes are
identically calibrated are not. Measured on a small 8-stage configuration: 20 identical
constructions gave 6 distinct ladders, with a maximum per-stage spread of 0.798 in S/N
(thresholds run 0.1 to `snr_final`), and seeding the legacy global RNG before each
construction changed nothing at all.

This is a reproducibility defect, not a correctness one. Everything downstream of the
time series is deterministic; the non-determinism is entirely in the *inputs* — the
noise realisation and the threshold ladder.

### The change

Each entry point takes an optional `seed`, accepting an `int`, a `np.random.Generator`,
or `None`. **`None` remains the default and still draws fresh entropy**, so no existing
caller sees different numbers. It is a pure addition.

Two details that are not just plumbing:

**`PulseSignalConfig` stores the generator** rather than re-creating it per call. A seed
therefore fixes the *sequence*: successive `generate*()` calls differ, and two runs of
the same script agree. Re-seeding inside each call would make a loop over one seeded
config draw the identical realisation every iteration.

**`DynamicThresholdScheme.run()` needed more than a constructor seed.** Both kernels
that draw randomness — `run_stage_legacy` and `pre_simulate_stage_folds` — are
`@njit(parallel=True)` and consume the generator inside a `prange`. A single shared
generator there leaves the result dependent on thread scheduling: seeded, on 14
threads, 5 runs produced 2 distinct ladders — enough to demonstrate that it happens,
and not offered as a rate. Each parallel iteration now gets its own generator, from
`SeedSequence(entropy, spawn_key=(istage,)).spawn(n)`, indexed by the loop variable —
so iteration `i` always uses generator `i`, and the result is independent of thread
order by construction.

### It also makes `run()` faster

Removing the shared generator removed a contention point in the parallel region.
`run()` over 32 stages, unseeded, best of 3:

| mode | before | after |
|---|---|---|
| `legacy` | 8.589 s | **2.172 s** |
| `improved` | 0.668 s | **0.391 s** |

We did not set out to change performance and would not have predicted the size of it;
it is reported as measured.

### Tests

`tests/test_rng_seeding.py`. For every site: the same seed reproduces, different seeds
differ, and **unseeded still varies** — that last one pinning that the default did not
quietly become deterministic. `run()` is covered in both `legacy` and `improved` modes,
and across thread counts.

`test_no_bare_default_rng_in_library` walks the AST of every module under `src/` and
fails on any argument-less `default_rng()`, so a newly added unseeded generator fails
the suite rather than going unnoticed. The list is enumerated from the source, not typed
into the test.

`rng_reproducibility.py` at the repo root is a standalone reproducer — no fixtures, no
data files, about 25 s. It demonstrates the defect, the fix, and the thread-count
independence, and enumerates the call sites from the installed source.

### Why the example-test tolerance fix is in the same PR

The example tests now pin a seed, and this PR also fixes a tolerance in
`test_example_ep_jerk.py` that was failing **5 times in 60 runs of that test file**
(unseeded). Those are not two changes. Pinning a seed does not merely *expose* that
tolerance bug — it would **entomb** it. A 1-in-12 flake becomes deterministic the moment
the seed is fixed, so a seed that happens to pass freezes the defect permanently out of
sight and leaves the suite looking healthier than it is. We know which seeds are hard
because we measured them, so this is concrete rather than hypothetical. Fixing the
tolerance is what makes the seeding honest rather than cosmetic.

The tolerance itself: `ACCEL_TOL = 1.0` was commented as a ">100x margin" against the
observed error spread, but it sits 8x *below* one `daccel` grid step. The recovered
error is bimodal — the best candidate lands in the correct acceleration cell (error
~0.008) or an adjacent one (~8.285, which *is* `daccel`) — so an absolute tolerance
under one grid step silently required an exact-cell hit. The check now runs against the
reported `daccel`: 5 failures in 60 runs of the file before, 0 in 60 after, with the
noise left unseeded in both measurements so it is the tolerance being tested and not the
pin.
