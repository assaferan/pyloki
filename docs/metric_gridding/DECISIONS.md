# DECISIONS.md — metric-based grid refinement

Conventions and choices, fixed once and referred to thereafter. Append a session log
entry at the end of every session.

## Fixed conventions (Phase 0)

These are **observed from the code**, not chosen, and any metric code must match them.

| # | convention | value | source |
|---|---|---|---|
| C1 | Axis ordering | reverse `k`: `[d_kmax, ..., d_2, d_1]` for the `poly_order` branchable axes | `core/taylor.py:52-64`, `:106-109` |
| C2 | Leaf layout | `(n_leaves, poly_order+2, 2)`; row `[-1]` = `f_0`/basis-flag, `[-2]` = `d_0` (never branched), `[-3]` = `d_1`, `[:-3]` = `[d_kmax..d_2]` | `poly_taylor_seed` |
| C3 | Column meaning | col 0 = value, col 1 = per-axis **half-width** | `core/common.py:50`, `core/taylor.py:109` |
| C4 | `Δt` sign | `delta_t = t_new − t_old`; `t_mat` lower-triangular, **unit diagonal** | `transforms.py::shift_taylor_params`, callers at `taylor.py:222,352` |
| C5 | Units | `d_k` in m·s⁻ᵏ, `f_0` in Hz, `C_VAL` in m/s | `utils/misc.py` |
| C6 | Phase units | **cycles**, and `phase = f0 * d / C_VAL` with **no 2π** | `poly_taylor_resolve_batch:228-234` |

## Choices made (Phase 0)

- **D1 — Base branch.** `metric-gridding` is cut from `6b11aba` = `upstream/main`, not
  from the `Chebyshev` working branch. Rationale: the plan's ground rule is "do not
  modify existing behaviour", so a clean upstream base keeps the eventual diff reviewable
  and upstreamable. Consequence: pruning segfaults on this base until PR #3 lands
  (see open question O1).

- **D2 — `2π` handling.** The codebase carries phase in cycles with no `2π`
  (C6). Rather than track a loose factor, `g` will be **defined so that
  `m = δᵀ g δ` is directly the dimensionless fractional S/N loss**, with any `(2π)²`
  folded into `g` itself. `poly_phase_metric` will state this in its docstring and
  `mismatch()` will return a pure number. *To be confirmed against Phase 1 test 4
  (empirical S/N loss) before being treated as settled.*

- **D3 — How to run anything in this worktree (important).** The repo's `.venv` lives in
  the **main** checkout and its editable install resolves `import pyloki` to
  `/Users/assaferan/Documents/GitHub/pyloki/src/pyloki` — i.e. **the main checkout's
  source, on whatever branch it happens to be**, not this worktree's. Running `pytest`
  here without care silently tests the wrong code. Verified. Always shadow it:

  ```sh
  PYTHONPATH=$PWD/src /Users/assaferan/Documents/GitHub/pyloki/.venv/bin/python -m pytest
  ```

  Confirmed to resolve to the worktree source and give `25 passed` at `2bc0455`.
  (A dedicated venv in the worktree would also work but duplicates ~1 GB of deps.)

## Open questions (carried from Phase 0, need a human answer)

- **O1 — Phase 3 blocker.** Every `prune_dyp_tree` call segfaults on this base
  (numba cannot box the heterogeneous stats dict; fixed on `fix-prune-segfault`,
  upstream PR #3, open). Phases 0-2 are unaffected. Merge the fix branch when Phase 2
  ends, or wait for the PR to land?

- **O2 — `m_max` bridge.** Match the ellipsoid to the `η`-box along the
  best-conditioned axis, or by matched volume? Neither is canonical, because the current
  criterion is sup-norm and the new one is mean-square and the conversion factor depends
  on `poly_order`. Proposal: `m_max_from_eta` computes **both**, Phase 1 test 5 reports
  both ratio sets, and the choice is made from the numbers.

- **O3 — Does `"metric"` inherit the `2**k` coarsening?** `poly_taylor_step_f`'s
  coarsening also feeds `B(s)` and the base grid `G0`, so it cannot simply be deleted.
  Reading of intent: `"metric"` should bypass it for leaf sizing (the metric *is* the
  sizing rule) while `G0` keeps whatever the FFA built. Needs confirmation — it changes
  what Phase 3's "Figure 7 analogue" is comparing against.

- **O4 — Harmonic weighting.** The plan asks which harmonic weighting best matches the
  boxcar scoring in `detection/scoring.py`. Deferred to Phase 1 test 4, which measures
  it rather than assuming.

## Session log

```
## 2026-09-09 — Phase 0
Done:
  - Worktree + branch `metric-gridding` created off 6b11aba (= upstream/main).
  - metric_PLAN.md committed to the branch.
  - docs/metric_gridding/00_pipeline_notes.md: all nine plan "Repository facts"
    verified against the code; two additions, one correction.
  - Full consumer trace of leaf column 1 (10 call sites, with file:line).
  - B(s) / threshold path documented, including that tiling_strategy already
    feeds B(s) so "metric" forces recalibration.
  - Current vs target invariant written down, with the sup-norm/mean-square
    non-equivalence made explicit.
Conventions fixed: C1-C6 above (all observed, not chosen); D1, D2, D3.
Corrections to the plan:
  - `validate` is a no-op in the Taylor path and consumes column 1 not at all.
  - `world_tree.py` never interprets column 1 (only column 0 and whole rows),
    so Phase 2 decision (a) vs (b) is less constrained than the plan assumes.
  - Phase 3 is blocked on PR #3 (pruning segfaults on this base).
Carried in from sibling branch (quantifies the plan's premise):
  - the `2**k` coarsening makes the current box 1.5x-15.5x looser than `eta`,
    growing with poly_order;
  - `aggressive` vs `conservative` is worth +0.584 of score in Taylor
    (5/5 replicates, p=0.005) for ~9-11x less pruning time.
Open questions: O1-O4 above.
Next session starts at: Phase 1 — src/pyloki/core/metric.py + tests/test_metric.py,
  after the human resolves O2/O3 (O1 and O4 can wait).
```
