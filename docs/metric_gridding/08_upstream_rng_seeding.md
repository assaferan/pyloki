# A pyloki run cannot be reproduced: `default_rng()` is unseeded in eight places

> ## ⚠ SUPERSEDED — owned by `fix-flaky-norm-isf`, 2026-09-22
>
> This draft and `rng_reproducibility.py` were written here before the fix existed, and
> the whole RNG story now lives on **`fix-flaky-norm-isf`**: the implementation, an
> extended five-check reproducer, and `FINDINGS_rng_seeding.md`. **Read those, not this.**
> Kept unedited below for provenance.
>
> **One claim in the text below is now known false.** It says shape (1) was implemented
> across all eight sites. Running this reproducer against the fix found that
> `DynamicThresholdScheme.run()` still produces different ladders when seeded:
> `run_stage_legacy` is `@njit(parallel=True)` (`thresholding.py:597`) and draws from the
> shared generator inside a `prange`, so a seed fixes the stream but not the order the
> threads consume it. Seeding cuts the spread from 0.798 to 0.080 and forcing one thread
> removes it entirely, which locates the residue in scheduling rather than in the seed.
> The honest summary is **seven of eight closed; the eighth needs per-iteration
> generators inside the kernel** — not "fixed".
>
> The framing that survives, and the reason this draft is kept at all:
> **`DynamicThresholdScheme` is the site that matters most, because a ladder is an input
> to a search.** That is still true, and it is now the site with a known-open gap.

**Status: draft, NOT POSTED, and now superseded.** Held pending assaferan's review, like
every other upstream item from this branch.

*All line numbers refer to `main` at `18d04b3`, and are enumerated from the source by the
reproducer rather than typed, so they cannot drift.*

## Summary

`np.random.default_rng()` is called with no argument in **eight places** across three
modules, and no public entry point accepts a `seed` or `rng`. Two consequences:

1. **A run cannot be repeated.** Same inputs, same configuration, different output.
2. **`np.random.seed` does not help**, because `default_rng` ignores the legacy global
   seed. This is what makes the defect easy to miss — the usual reflex appears to work
   and changes nothing.

The site that matters most is `DynamicThresholdScheme.__init__`
(`detection/thresholding.py:740`), because a threshold ladder is an **input** to a
search, not an output of one. Two searches a user believes are identically calibrated
are not.

This is a reproducibility defect, not a correctness bug — see **Scope** below.

## Evidence

`DynamicThresholdScheme` built with identical arguments, run 20 times, ladder
backtracked at `P_d = 0.1` each time (8 stages, `B(s) = 2`, `ntrials = 1024`,
`nthresholds = 100` — small so the check runs in seconds; nothing depends on the
configuration):

|  | distinct ladders / 20 | max per-stage spread | mean spread |
|---|---|---|---|
| as shipped | **6** | 0.798 | 0.110 |
| after `np.random.seed(42)` before each construction | **6** | 0.798 | 0.110 |

Seeding the legacy global RNG changes nothing — not the ladders, not even the summary
statistics. The figures themselves vary between invocations of the reproducer, which is
the symptom rather than a caveat.

The spread is on a deliberately tiny configuration, chosen for speed, and we have not
measured how it scales to a production one. The claim rests on the table as it stands:
**the ladders differ, and no user-accessible seed makes them not differ.** How much that
moves a detection probability in practice will depend on `ntrials` and on the branching
pattern, and we are not asserting a magnitude here.

## The eight sites

Enumerated from the installed source by the reproducer:

| site | what it is |
|---|---|
| `detection/thresholding.py:740` | `DynamicThresholdScheme.__init__` → `self.rng`, used at nine call sites within the class |
| `detection/thresholding.py:1046` | `determine_scheme()`, module-level |
| `detection/thresholding.py:1094` | `evaluate_scheme()`, module-level |
| `simulation/pulse.py:220` | `generate_simple()` — the injected noise |
| `simulation/pulse.py:228` | `generate_noise()` |
| `simulation/pulse.py:261` | `generate_old()` |
| `simulation/pulse.py:338` | `generate()` — the current path, and the amplitude is calibrated against the realised noise (`calibrate_scale_on_folds`), so the *signal* varies too |
| `sensitivity/sim_ffa.py:198` | the sensitivity simulation's `__init__` |

**It is a family, not a site.** Fixing `thresholding.py:740` alone leaves a library that
still cannot reproduce a run. We mention this because the threshold site is the one we
hit hardest and it would be the natural thing to patch in isolation.

## Why it is worth fixing rather than working around

Every workaround is available and we have used all of them — persist the time series and
re-load it, generate the ladder once and commit the array, cache `determine_scheme`'s
output. They work. But:

- **They are per-caller and undiscoverable.** We arrived at this defect independently
  three times, in three separate pieces of work, at three different sites, each time
  after a comparison behaved oddly rather than from reading the code.
- **The failure is silent and looks like data.** A comparison of two configurations run
  against independently-drawn ladders looks exactly like a comparison of two
  configurations. Nothing errors, nothing warns; the difference is absorbed into the
  result.
- **It blocks the ordinary regression test.** "Run this configuration, assert the
  output" is not expressible today. The same root cause is already in the test suite:
  `tests/test_maths.py:10` draws from an unseeded `default_rng`, so
  `test_norm_isf_func` samples a different point every run and fails spuriously when it
  lands where the tolerance does not hold.

## Scope — what this does NOT claim

- It does **not** claim any published result is wrong. Ratcheting noise into a ladder
  changes which run you get, not whether the search works.
- It does **not** claim the search itself is non-deterministic. It is not: everything
  downstream of the time series is deterministic, and we have verified that two processes
  given a byte-identical time series produce bit-identical sorted score vectors. The
  non-determinism is entirely in the *inputs* — noise realisation and threshold ladder.
- It does **not** claim a fix is needed to make paired comparisons possible. They are
  already possible by persistence, and we said so in our own notes before proposing this.
  The argument here is about reproducibility for a user who does not already know the
  workaround, which is a different goal from pairability.

## Possible shapes, none of them obviously right

We have no strong view, and the choice interacts with API taste:

1. **A `seed` or `rng` keyword on the public entry points**, defaulting to `None` and
   preserving today's behaviour. For `DynamicThresholdScheme` the plumbing already
   exists — everything downstream of `__init__` already takes `rng` as an explicit
   parameter (`thresholding.py:50, 393, 445, 544, 610`), so this is a small change at
   that site.
2. **Remove the randomness where it is not needed.** When we made
   `test_norm_isf_func` deterministic in our own tree we did it this way — replacing the
   sampling with fixed evaluation points rather than seeding it — and it was clearly the
   better fix at that site. Some of the eight may be the same.
3. **A module-level default generator** that a user can set once. Least invasive at the
   call sites, worst for concurrent use.

~~**We have implemented (1) in our own tree** across all eight sites~~ — **corrected:
seven of eight.** A `seed: int | np.random.Generator | None = None` keyword defaulting to
`None`, so current behaviour is unchanged for every existing caller; but
`DynamicThresholdScheme.run()` is not made reproducible by it, for the parallel-kernel
reason in the banner above. It is not proposed here as a patch
because the choice between the three shapes is yours; say which you prefer and we will
open it as a PR, with a regression test asserting that two identical constructions agree.

One design point that came out of doing it, in case it informs the choice:
`PulseSignalConfig` holds its `Generator` for the object's lifetime, so a seed fixes the
*sequence* rather than each call — two `generate()` calls on one config still differ,
while two configs built from the same seed agree call-for-call. That is the right
semantics for common random numbers, but it is a behaviour worth naming in the docstring
rather than leaving a caller to discover.

## Reproducing

    python docs/metric_gridding/rng_reproducibility.py

Self-contained — no data files, no fixtures, a few seconds. It prints the table above and
enumerates the eight sites from the installed source, so the list cannot go stale against
the code. Script:
[`rng_reproducibility.py`](https://github.com/assaferan/pyloki/blob/injection-design/docs/metric_gridding/rng_reproducibility.py).

**Limitations, stated plainly.** The measured table is one small configuration at one
`P_d`. It establishes that identical arguments give different ladders and that
`np.random.seed` does not prevent it; it does not quantify how much the resulting
detection probability moves on a production configuration, and nothing here should be
read as such a quantification.
