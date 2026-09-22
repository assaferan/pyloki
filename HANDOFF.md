# flaky-test — brief for an incoming agent

Read this first. `FINDINGS_norm_isf_func.md` is the record of what the branch already
did; this is what is left and why the scope is wider than the branch name.

## What this branch already did

Two commits, both about `tests/test_maths.py`:

- `00e04e5` replaced random sampling of `norm_isf_func` with fixed off-node points and
  tolerances measured on a 10^6-point grid.
- `500a99a` recorded the table-edge defects found while doing it.

Note it used **two different fix shapes**, and both are legitimate:
*remove* the randomness (the fixed points) and *seed* it (`test_maths.py:215`,
`default_rng(42)`). Which one is right depends on whether the randomness is buying
coverage. Fixed points that are chosen well beat a seeded draw; a seeded draw beats
fixed points when the input space is large and you want the generator to explore it.

## The actual scope: the library never seeds anything

Eight unseeded `np.random.default_rng()` constructions in `src/`, verified on this
branch today:

    src/pyloki/simulation/pulse.py:220, 228, 261, 338
    src/pyloki/detection/thresholding.py:740, 1046, 1094
    src/pyloki/sensitivity/sim_ffa.py:198

None of them takes a `seed` or `rng` argument from its caller. `np.random.seed` does
**not** help — `default_rng` ignores the legacy global seed, which is what makes this
bug non-obvious to anyone reporting it.

This has been diagnosed independently on three branches and fixed on none:
`injection-design` (`docs/metric_gridding/05_injection_design.md` §6.1),
`Chebyshev` (`HANDOFF.md:170` and `:490`), and `metric-gridding`
(`docs/metric_gridding/README.md:166`). Each one **worked around** it — persist the
time series, commit the ladder, cache the scheme. Nobody proposed a library change.
There is no issue and no PR upstream; `#15 #12 #9` open and `#8` closed are not this.

## The finding that sets the scope

`tests/test_example_ep_jerk.py::test_recovers_injected_acceleration` was seen to fail
once in four runs. That test file contains **no RNG of its own** — grep it. It calls
`PulseSignalConfig.generate(shape="gaussian")` at line 96, and the noise is drawn inside
the library at `pulse.py:261`/`:338`.

So that flake **cannot be fixed at the test level**. No amount of seeding in the test
file reaches the generator that produced the noise. That is the argument for fixing the
library rather than the tests, and it is checkable in two greps.

`test_maths.py` was fixable at the test level because the draw was in the test. This one
is not. Do not assume the rest of the suite resembles `test_maths.py`.

## Suggested approach, not binding

Thread an optional `seed: int | None = None` (or an `rng` argument) from each public
entry point to its `default_rng` call. The plumbing in `thresholding.py` already exists
— everything downstream of `__init__` takes `rng` as an explicit parameter
(`:50, :393, :445, :544, :610`), so only the constructor needs to change. `pulse.py` has
four sites and may need more thought; check whether they should share one generator.

Default behaviour must not change: `seed=None` keeps today's fresh entropy, so no
existing caller sees different numbers. That makes this a pure addition, which is the
easiest kind of change for a maintainer to accept.

Add a regression test that constructs the same object twice with the same seed and
asserts the outputs are identical. A fix without that test will silently rot.

## Constraints

- **Nothing goes upstream.** Not an issue, not a PR to `pravirkr/pyloki`. Prepare, commit
  locally, stop. Pushing a branch to the `assaferan` fork is ordinary and allowed;
  opening anything on the maintainer's tracker is not, until assaferan has read it.
- Do not use the issue tracker to talk to other sessions. `SendMessage` after
  `ListAgents`, or the shared memory store.
- `pyloki` in this worktree resolves to the **main** checkout unless you set
  `PYTHONPATH=/Users/assaferan/Documents/GitHub/pyloki/worktrees/flaky-test/src`.
  Verify it, do not trust it. There is no bare `python`, no `uv`, no `timeout`;
  use `/Users/assaferan/Documents/GitHub/pyloki/.venv/bin/python`.
- A stale numba cache under `src/pyloki/**/__pycache__` can make `test_example_ffa.py`
  error at setup. Clear it before believing a phantom failure.
- No Claude attribution trailers on commits — no `Co-Authored-By: Claude`, no session
  link, no "Generated with" footer.
- Memory: this branch writes only `wt-fix-flaky-norm-isf-*` facts, via
  `/Users/assaferan/Documents/GitHub/pyloki/.claude/memctl`. Never hand-edit `MEMORY.md`.

## Measure the flake before and after

The rate is the claim. `test_example_ep_jerk` was 1 failure in 4 — that is an estimate
with enormous error bars, not a rate. Run it enough times to bound it, state the
denominator, and re-measure after the fix. This project's recurring defect is quoting a
fraction whose denominator or predicate has moved; see `upstream-posts-held-for-review`
in the shared memory before writing any number down.
