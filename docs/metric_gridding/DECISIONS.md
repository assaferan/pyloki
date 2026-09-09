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

- **D3 — This worktree has its own `.venv`; always use it.**

  ```sh
  .venv/bin/python -m pytest        # from worktrees/metric-gridding
  ```

  Baseline at `fb3d413`: **25 passed**.

  *Why it needs saying.* The repo's top-level `.venv` (in the parent checkout) is an
  editable install that resolves `import pyloki` to the **parent checkout's** `src/` —
  i.e. whatever branch that checkout happens to be on, not this worktree's. Using it
  from here silently tests the wrong source. Verified, not hypothetical.

  An earlier revision of this entry recommended `PYTHONPATH=$PWD/src <parent>/.venv/...`
  instead. Do not use that: the worktree-isolation guard refuses inline `PYTHONPATH`
  before a `python` invocation, so the advice was unusable. Hence the local venv
  (built from `/opt/homebrew/bin/python3.13`, which has tkinter —
  `detection/schemes.py` imports it at module scope). `.venv` is already gitignored.

## Choices made (Phase 1)

- **D4 — `m_max_from_eta` returns every bridge, not one number** (resolves O2, per the
  human's answer "compute both"). It returns a `MMaxBridge` with:
  `per_axis` (the `m_max` that matches the `eta`-box exactly on each axis),
  `axis_tight` = `min(per_axis)` (ellipsoid **inscribed** in the box — never wider than
  the box on any axis), `axis_loose` = `max(per_axis)` (box inscribed in the ellipsoid),
  and `volume` (matched ellipsoid/box volume). Nothing in the library picks one; callers
  and Phase 3 state which they used. Rationale: the old criterion is sup-norm and the new
  one mean-square, so no single bridge is canonical (notes §5).

- **D5 — `"metric"` does not inherit the `2**k` coarsening** (resolves O3, human's
  answer). Consequences, to hold to in Phase 2: leaf sizing under `"metric"` comes from
  the metric alone; anywhere the metric path needs the pre-existing box for comparison it
  uses `poly_taylor_step_d_vec(..., use_cheby=False)`. The coarsened box remains the
  reference for the *existing* strategies, and `G0` keeps whatever the FFA built. Phase 1
  test 5 therefore reports ratios against **both** the un-coarsened box (the like-for-like
  comparison) and the coarsened one (what the code actually uses today), since they differ
  by `2**k` per axis.

- **D6 — Mismatch convention: `m = 1 − A/A₀` (amplitude), so `g = 2π² · Cov`.**
  For a single harmonic, `A/A₀ = |⟨e^{2πiΔΦ}⟩| ≈ 1 − ½·Var(2πΔΦ) = 1 − 2π²·Var(ΔΦ)`
  with `ΔΦ` in cycles. Defining `g = 2π² (f₀/c)² Cov(τ^{kᵢ}/kᵢ!, τ^{kⱼ}/kⱼ!)` makes
  `m = δᵀgδ` the fractional **amplitude** (S/N) loss directly, with the `2π²` folded in
  as D2 requires.

  **The plan is internally inconsistent here and this deviates from one half of it.**
  `poly_phase_metric`'s docstring says "fractional S/N loss" (amplitude), but Phase 1
  test 4 says to compare `(S/N_mismatched / S/N_true)²` against `1 − m`, which is the
  **power** convention and larger by a factor of 2 (`m_power = 4π²·Var`). I chose
  amplitude because the plan cites Owen (1996) / Allen et al. (2013), where the match is
  the normalised overlap and `m = 1 − match` is amplitude-like. `test_metric.py` checks
  the amplitude relation and *also* reports the power one, so the factor of 2 is visible
  rather than buried. **Flag for the human:** if the intended convention is power, change
  the `2 * np.pi**2` constant in `metric.py` to `4 * np.pi**2` and halve every `m_max`.

- **D7 — `nbins` is accepted but unused by `poly_phase_metric`.** `g` does not depend on
  `nbins`; it enters only through `m_max_from_eta` (via `eta/nbins`). The argument is
  kept because the plan's API specifies it. It could be dropped in Phase 4.

- **D8 — Phase 1 test 5 ratios (the exit criterion), and what they say.**
  Ellipsoid axis extent ÷ box half-width, per axis, in leaf order
  `[d_kmax ... d_1]`, at `eta=1`, `nbins=64`, `f0=142.86 Hz`, interval `[0, 67.109] s`,
  `t_ref=0`, against the **un-coarsened** box (D5):

  | poly_order | `axis_tight` (inscribed) | `volume`-matched |
  |---|---|---|
  | 2 | `[0.97, 1.00]` | `[2.22, 2.29]` |
  | 3 | `[0.66, 1.00, 0.43]` | `[4.86, 7.40, 3.18]` |
  | 4 | `[0.50, 1.00, 0.67, 0.16]` | `[10.5, 21.2, 14.1, 3.47]` |
  | 5 | `[0.40, 1.00, 0.91, 0.36, 0.058]` | `[22.5, 56.5, 51.2, 20.2, 3.26]` |

  No units or ordering bug: everything is within a factor of a few at `poly_order=2`,
  where box and ellipsoid nearly coincide. The interesting part is how fast they diverge.
  At `poly_order=5` a volume-matched ellipsoid reaches **56x** the box half-width on one
  axis, while an inscribed one is **0.058x** on another. Both readings say the same
  thing: the mismatch ellipsoid is strongly correlated (ill-conditioned) and an
  axis-aligned box is a poor proxy for it, increasingly so with `poly_order`. That is the
  quantitative case for this project — a box either wastes volume along the long axes or
  gaps along the short ones, which is exactly the `conservative`-vs-`aggressive`
  dilemma, and why the gap cost measured in the notes grows with `poly_order`.

- **D9 — Harmonic weighting is NOT optional; measured (bears on O4).** Plan test 4's
  duty-cycle diagnostic, boxcar-scored folded profile vs the single-harmonic metric, at
  `m_target = 4e-3`, `poly_order=4`:

  | ducy | `A/A₀` | actual loss ÷ `m` |
  |---|---|---|
  | 0.50 | 0.999108 | 0.22 |
  | 0.20 | 0.991343 | 2.16 |
  | 0.10 | 0.969917 | **7.52** |
  | 0.05 | 0.851051 | **37.2** |

  A fundamental-only metric **under-predicts the real loss by ~7.5x at `ducy=0.1`** (the
  plan's `ref_ducy`) and ~37x at 0.05, and over-predicts by ~4x for a broad `ducy=0.5`
  pulse. The scaling is roughly `ducy**-2.2`, consistent with an `n**2`-weighted sum over
  the ~`1/ducy` harmonics a narrow boxcar carries.

  Consequence: `m_max` cannot be a single global constant. Either `g` carries a harmonic
  weighting factor set by the search's `ducy_max`, or `m_max` is chosen per duty cycle.
  **This needs a decision before Phase 2 sizes any leaf** — see O5.

  Caveat on the measurement: folding with wrong parameters is modelled here as circular
  convolution of the profile with the histogram of `dPhi` over the observation, which
  assumes the drift is slow compared with one rotation. It is a good approximation in
  this regime but is not the full FFA path, so treat the factors as indicative to ~10%,
  not exact.

## Open questions (carried from Phase 0, need a human answer)

- **O5 — How should the harmonic weighting enter?** (new, from D9). Options: (a) fold a
  `ducy`-dependent factor into `g` inside `poly_phase_metric`, so `m` predicts boxcar
  loss directly for the configured `ducy_max`; (b) keep `g` single-harmonic and make
  `m_max` a per-duty-cycle budget chosen by the caller; (c) sum explicitly over
  harmonics with the profile's power spectrum, which is most faithful and most
  expensive. Recommendation: (a) for Phase 2, because leaf sizing needs one number and
  the search already has `ducy_max` in `PulsarSearchConfig`; revisit if Phase 3 shows
  the `ducy` spread matters.

- **O1 — Phase 3 blocker.** Every `prune_dyp_tree` call segfaults on this base
  (numba cannot box the heterogeneous stats dict; fixed on `fix-prune-segfault`,
  upstream PR #3, open). Phases 0-2 are unaffected. Merge the fix branch when Phase 2
  ends, or wait for the PR to land?

- ~~**O2 — `m_max` bridge.**~~ **RESOLVED (human, 2026-09-09): compute both.**
  See D4.

- ~~**O3 — Does `"metric"` inherit the `2**k` coarsening?**~~
  **RESOLVED (human, 2026-09-09): no, `"metric"` must not inherit it.** See D5.

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

## 2026-09-09 — Phase 1
Done:
  - O2 and O3 answered by the human; recorded as D4 and D5.
  - src/pyloki/core/metric.py: poly_phase_metric, mismatch,
    ellipsoid_axis_extents, shift_matrix, transform_metric, cholesky_factor,
    m_max_from_eta (returning MMaxBridge). Pure NumPy; nothing wired into the
    search, no existing behaviour touched.
  - tests/test_metric.py: 41 tests covering all five plan tests. Full suite
    66 passed.
Conventions fixed: D4 (m_max bridge returns all variants), D5 (no 2**k
  inheritance), D6 (amplitude convention, g = 2*pi**2 * Cov), D7 (nbins unused),
  D8 (test-5 ratios), D9 (harmonic weighting measured).
Deviations from the plan, deliberate:
  - D6: the plan is self-inconsistent on the mismatch convention (docstring says
    amplitude, test 4 says power, factor of 2). Chose amplitude per Owen/Allen;
    the test reports both so the factor is visible. Flagged for the human.
  - shift_matrix() added to the API (not in the plan) so the T construction is
    in one place and testable against transforms.shift_taylor_params directly.
    The restriction to the leading poly_order block is exact, not an
    approximation, because the full matrix is lower triangular with d_0 last;
    there is a test for that.
Findings worth acting on:
  - D8: box and mismatch ellipsoid diverge fast with poly_order (56x on one axis
    at poly_order=5). Quantifies why the tiling dilemma worsens with order.
  - D9: single-harmonic m under-predicts real boxcar loss by ~7.5x at ducy=0.1.
    m_max cannot be one global constant. Blocks Phase 2 leaf sizing.
Open questions: O1 (Phase 3 needs PR #3), O5 (new, harmonic weighting) — O5
  should be answered before Phase 2 starts.
Next session starts at: Phase 2, after O5. Phase 2 design decisions to make
  first are listed in metric_PLAN.md (metric storage, lattice, meaning of
  column 1); note Phase 0 found world_tree.py does not interpret column 1, so
  option (a) is cheaper than the plan assumed.
```
