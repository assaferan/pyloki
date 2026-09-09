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
  and upstreamable. The consequence that pruning segfaulted on this base is **gone**:
  PR #3 merged upstream as `2b4b80c` and this branch is rebased onto it (O1).

- **D2 — `2π` handling.** The codebase carries phase in cycles with no `2π`
  (C6). Rather than track a loose factor, `g` will be **defined so that
  `m = δᵀ g δ` is directly the dimensionless fractional S/N loss**, with any `(2π)²`
  folded into `g` itself. `poly_phase_metric` will state this in its docstring and
  `mismatch()` will return a pure number. **Confirmed in Phase 1** (D6, D10): `m` is
  the fractional amplitude loss and matches a direct measurement to 0.5% once the
  harmonic weight is applied.

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
  rather than buried. **CONFIRMED by the human (2026-09-09): amplitude is fine.** D6 is
  settled; no change to the `2 * np.pi**2` constant.

- **D7 — `nbins` is used only when `ducy` is given.** `g` is independent of `nbins` for
  the single-harmonic metric, and `poly_phase_metric` raises if `ducy` is passed without
  it. Superseded in part by D10: with a `ducy`, `nbins` sets the harmonic cutoff and the
  argument is no longer dead.

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

- **D9 — Harmonic weighting is NOT optional (bears on O4).** A fundamental-only metric
  badly under-predicts the real S/N loss, because folding with wrong parameters
  multiplies harmonic `n` by `kappa_n = <exp(2*pi*i*n*dPhi)>`, so the penalty grows as
  `n**2`. Measured under-prediction at `poly_order=4`, ideal matched filter, against the
  single-harmonic `m`: a factor of `~25` at `ducy=0.1`, tracking `harmonic_weight` to
  within 15%. So `m_max` cannot be a single global constant. Resolved as O5(a) / D10.

  > **CORRECTED.** The first version of this entry gave a table of `loss/m` = 0.22 /
  > 2.16 / 7.52 / 37.2 for `ducy` = 0.5 / 0.2 / 0.1 / 0.05, and those numbers were
  > **wrong**. They came from modelling mis-parameterised folding as circular convolution
  > of the profile with a *histogram* of `dPhi`, which has a hard resolution floor of one
  > phase bin: the rms `dPhi` there was ~0.9 bins, right at the limit, and at smaller
  > offsets the kernel collapses to a delta function and the model reports *exactly zero*
  > smearing. The replacement computes `kappa_n` directly from the time samples with no
  > binning (`tests/test_metric.py::matched_filter_loss`). Two further errors were found
  > and fixed while chasing this: the comparison must mean-subtract `dPhi`, because
  > `poly_phase_metric` projects out the constant-phase mode (otherwise `Re(kappa_n)`
  > charges the mean as loss and the ratio comes out ~1.45x high, and integer bin-shift
  > maximisation does not fix it because the mean is a small fraction of a bin); and the
  > boxcar filter is the wrong validation target (see D10).

- **D10 — `harmonic_weight(nbins, ducy, weighting=...)`, defaulting to `"power"`**
  (resolves O5(a), the human's answer). `g` is scaled by a single scalar `<n**2>`:

  | ducy | `"power"` (default) | `"cross"` |
  |---|---|---|
  | 0.50 | 1.59 | 0.685 |
  | 0.20 | 6.99 | 2.44 |
  | 0.10 | 25.4 | 12.1 |
  | 0.05 | 102.2 | 37.0 |

  - `"power"` weights by `|p_n|**2` (ideal matched filter). **Validated**: with the
    corrected measurement, `loss / m` is **0.995-0.998** across `ducy` 0.05-0.5.
  - `"cross"` weights by `Re(b_n^* p_n)`, the cross-spectrum with the boxcar filter the
    search actually scores with, and is 2-3x smaller — nominally closer to the pipeline.

  `"power"` is the default despite `"cross"` being nominally more faithful, for two
  reasons. It is the one that can be validated: a boxcar has sinc sidelobes, so
  `Re(b_n^* p_n)` is **negative** for some `n`, attenuating those harmonics can
  *increase* the score, and boxcar loss is therefore **not monotonic in smearing** —
  measured directly it returns a small *negative* loss at `ducy=0.1`. And for a
  *covering* criterion over-predicting the loss is the safe direction, since it yields
  smaller leaves. Phase 3's recalibration is the place to decide whether the 2-3x that
  `"cross"` would claim is worth taking.

## Open questions (carried from Phase 0, need a human answer)

- ~~**O5 — How should the harmonic weighting enter?**~~
  **RESOLVED (human, 2026-09-09): (a), fold it into `g`.** Implemented as D10. Note the
  implementation ended up equivalent to the plan's option (c) in spirit — the factor
  *is* computed from the profile's spectrum — but collapsed to one scalar `<n**2>` up
  front, so there is no per-harmonic cost at branch time.

- **O6 — `"power"` or `"cross"` weighting for production?** (new, from D10). `"power"`
  is the validated, conservative default and `"cross"` claims a 2-3x looser grid but
  cannot be validated in isolation because boxcar loss is non-monotonic in smearing.
  Phase 3 should settle this against real injection-recovery rather than a model. Not
  blocking Phase 2: it is one keyword.

- ~~**O1 — Phase 3 blocker.**~~ **RESOLVED (2026-09-09): upstream PR #3 merged**
  (`2b4b80c`), so the fix arrived by fast-forward rather than needing a branch merge.
  This branch was rebased onto it. Pruning now runs on this base and `tests/test_prune.py`
  comes in from upstream: the suite is **79 passed**. Phase 3 is unblocked.

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
Open questions: O6 (new: power vs cross weighting, not blocking).
Resolved after the fact: O1 (PR #3 merged upstream as 2b4b80c; branch rebased,
  79 passed) and D6 (human confirmed amplitude).
Next session starts at: Phase 2. Phase 2 design decisions to make
  first are listed in metric_PLAN.md (metric storage, lattice, meaning of
  column 1); note Phase 0 found world_tree.py does not interpret column 1, so
  option (a) is cheaper than the plan assumed.
```
