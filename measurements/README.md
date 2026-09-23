# Measurement harness

The scripts behind the numbers in `FINDINGS_rng_seeding.md`. Committed on the working
branch only — they are **not** part of PR #16 or #17 and should not be.

Run everything with the venv python and an explicit `PYTHONPATH`; `import pyloki` in a
worktree silently resolves to the MAIN checkout otherwise, which for RNG work is its own
trap.

    PY=/Users/assaferan/Documents/GitHub/pyloki/.venv/bin/python
    WT=/Users/assaferan/Documents/GitHub/pyloki/worktrees/flaky-test

## `margin_sweep.py` — which realisations are hard

The tool that chose every seed in the `SEED_SWEEP` tuples. Reports error/tolerance per
predicate per seed, so a ratio of 1 is a failure and the top of the ranking is where the
tolerance is actually being tested.

    PYTHONPATH=$WT/src:$WT $PY measurements/margin_sweep.py accel $(seq 1 40)
    PYTHONPATH=$WT/src:$WT $PY measurements/margin_sweep.py ffa $(seq 1 40)
    PYTHONPATH=$WT/src:$WT $PY measurements/margin_sweep.py circular $(seq 1 20)

Note it needs the worktree root on `PYTHONPATH` as well as `src`, because it imports the
test modules to reach their `_run_search` / `_run_case` entry points.

**Use this before adding a seed sweep to any test.** Picking seeds without measuring
costs the same runtime and covers less; the value of a sweep is concentrated entirely in
the near-miss realisations.

## `run_flake.sh` — the flake rate of a test file

    measurements/run_flake.sh <treedir> <nruns> <outfile>

Runs a test file N times and records pass/fail per run. It takes a **tree directory**
rather than running in place, which is the point: copy the tree first
(`git archive HEAD | tar -x -C <dir>`) so `src/` can be edited while a measurement is in
flight without moving the denominator underneath it.

Two rules that cost time to learn:

- Never clear the numba cache mid-measurement. One warm cache per tree, throughout.
- State the predicate with the number. These runs measure "one run of the file in which
  at least one test failed" — a per-file-run rate, not a per-test one.

## `bench.py` — `run()` timing

    PYTHONPATH=$WT/src $PY measurements/bench.py legacy
    PYTHONPATH=$WT/src $PY measurements/bench.py improved

Backs the 8.589 s -> 2.172 s figure in #16. Compare against a tree at the pre-fix commit
by pointing `PYTHONPATH` at that tree's `src`.
