# flaky-test — brief for an incoming agent

The library-wide RNG work this file used to describe as "what is left" is **done**.
Read `FINDINGS_rng_seeding.md` for the record, `FINDINGS_norm_isf_func.md` for the
earlier `test_maths.py` work, and this file for what remains.

## State of the branch

Four commits, in two pieces of work.

**`test_maths.py`** — `00e04e5` replaced random sampling of `norm_isf_func` with fixed
off-node points and tolerances measured on a 10^6-point grid; `500a99a` recorded the
table-edge defects found underneath it. Both `norm_isf_func` and
`chi_sq_minus_logsf_func` now use fixed points, so the shared class-attribute generator
is gone and no test in the suite draws from an unseeded RNG.

**The library** — `d99ae3e` gave all eight `default_rng()` sites an optional `seed`
(`int`, `Generator`, or `None`), with `None` still the default and still fresh entropy,
so it is a pure addition. `4598f9a` pinned `SEED = 42` in the six test modules that drive
those generators, and fixed the tolerance that the seeding exposed.

`tests/test_rng_seeding.py` guards it, including
`test_no_bare_default_rng_in_library`, which walks the AST of every module under `src/`
and fails on any argument-less `default_rng()`. A new unseeded generator fails the suite.

## The finding worth carrying forward

The `ep_jerk` flake was **not** an RNG defect. Seeding only made it reproducible; the
defect was a tolerance sized against the spread of a few passing runs rather than against
the grid. The recovered acceleration error is bimodal — right cell (~0.008) or one cell
off (~8.285, which *is* `daccel`) — so `ACCEL_TOL = 1.0` read as a ">100x margin" while
sitting 8x below one grid step.

Generalise this before sizing any other recovery tolerance in this repo: compare it to
the reported per-parameter uncertainty, not to the observed spread. See
`wt-fix-flaky-norm-isf-accel-tolerance-is-bimodal` in shared memory.

## What is left

- **Nothing has gone upstream**, and nothing should until assaferan has read it. There is
  a drafted note on the `injection-design` branch
  (`docs/metric_gridding/08_upstream_rng_seeding.md`) with its own reproducer; it was
  written before this fix existed and offers to prepare shape (1), which is what was
  built here. If assaferan clears an upstream post, that draft and these two FINDINGS
  files are the material, and the offer should become "here is the PR".
- `FINDINGS_norm_isf_func.md` records real library defects still unfixed —
  `norm_isf_func(-0.5)` returning ~+28 sigma, reachable from `scoring.py:654` on ordinary
  noise, and `norm_isf_func(0)` returning NaN. Those are correctness bugs, unlike
  anything in this file, and they are the stronger upstream item.
- The accel one-cell miss still happens on roughly 1 realisation in 12. The test now
  accepts it because the search never claimed better. Whether the search *should* do
  better at this configuration is unexamined.

## Constraints

- **Nothing goes upstream.** Not an issue, not a PR to `pravirkr/pyloki`. Pushing a
  branch to the `assaferan` fork is ordinary and allowed.
- Do not use the issue tracker to talk to other sessions. `SendMessage` after
  `ListAgents`, or the shared memory store.
- `pyloki` in this worktree resolves to the **main** checkout unless you set
  `PYTHONPATH=/Users/assaferan/Documents/GitHub/pyloki/worktrees/flaky-test/src`.
  Verify it, do not trust it. There is no bare `python`, no `uv`, no `timeout`;
  use `/Users/assaferan/Documents/GitHub/pyloki/.venv/bin/python`.
- A stale numba cache under `src/pyloki/**/__pycache__` can make `test_example_ffa.py`
  error at setup. Clear it before believing a phantom failure — but never clear it
  *during* a rate measurement, or the denominator moves underneath you.
- No Claude attribution trailers on commits.
- Memory: this branch writes only `wt-fix-flaky-norm-isf-*` facts, via
  `/Users/assaferan/Documents/GitHub/pyloki/.claude/memctl`. Never hand-edit `MEMORY.md`.

## Measuring a rate here

Full suite: 159 tests. One run of `test_example_ep_jerk.py` is ~15 s warm, ~48 s cold.
State the denominator and the predicate with any fraction you quote; see
`upstream-posts-held-for-review` in shared memory. The scripts used for the measurements
in `FINDINGS_rng_seeding.md` are in this session's scratchpad, not committed — they run
a tree copy so `src/` can be edited while a measurement is in flight, which is worth
reproducing rather than measuring against a tree you are also changing.
