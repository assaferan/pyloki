# metric_PLAN.md — Metric-based grid refinement for pyloki

## Purpose

Replace the axis-aligned box representation of a template's validity region
with the parameter-space metric (mismatch tensor), so that re-centring the
Taylor expansion at each Extreme Pruning stage no longer forces a choice
between redundant coverage (conservative/quadrature tiling) and coverage gaps
(aggressive tiling).

Reference: Kumar & Zackay (2026), "Coherent Signal Detection with Pruning I",
Sections 3.4, 5.2.2–5.2.4, and the closing remark of 5.2.4 ("A complete
solution likely requires replacing fixed coordinate spacings with a local
metric-based mismatch criterion."). Background on the metric formalism:
Owen (1996); Allen et al. (2013). Verify citations before relying on them.

## Ground rules

- Work in `pyloki` (Numba reference implementation), not `loki` (C++/CUDA).
- Do not modify existing behaviour. Every change goes behind a new
  `tiling_strategy="metric"` option in `config.py`. The existing
  `aggressive`, `quadrature`, and `conservative` strategies must produce
  bit-identical results before and after this work.
- Tests before integration. Phase 1 must be green before Phase 2 starts.
- One phase (or less) per session. Update `DECISIONS.md` at the end of every
  session with any convention chosen (ordering, signs, units).
- When unsure about a convention already in the codebase, read the code and
  write down what it does. Do not guess.

## Repository facts (verified by reading the code; re-verify in Phase 0)

- A leaf is `np.ndarray` of shape `(n_leaves, poly_order + 2, 2)`.
  Column 0 = parameter value, column 1 = per-axis half-width (`dparam`).
  Rows `[:-2]` = kinematic params in reverse order `[d_kmax, ..., d_2, d_1]`,
  row `[-2]` = `d_0`, row `[-1]` = `f_0` (value) and basis flag (col 1).
- Branching is per-axis and independent:
  `core/taylor.py::poly_taylor_branch_batch` calls
  `utils/psr_utils.py::branch_param_padded` on each axis, then a padded
  Cartesian product (`np_utils.cartesian_prod_padded`).
- Step sizes come from `psr_utils.poly_taylor_step_d_vec` →
  `poly_taylor_step_f`, including the Chebyshev coarsening factor
  (`use_cheby=True`).
- The three tiling strategies live in
  `utils/transforms.py::shift_taylor_errors`. The shift matrix `T(Δt)` is
  built there from `delta_t**powers / fact(powers)` (lower triangular).
- `core/taylor.py::poly_taylor_transform_batch` calls
  `transforms.shift_taylor_full`, which shifts values and errors together.
- Dynamic-programming plumbing: `dynamic/dyn_poly_taylor.py` (StructRef
  wrappers for load/seed/branch/validate/resolve/shift_add/score/transform/
  ascend/report). `utils/world_tree.py` holds the candidate buffer.
- Threshold schemes and the branching factor `B(s)`:
  `detection/schemes.py`, `detection/thresholding.py`.
- Existing tests are minimal (`tests/cheby.py`, `mat_inv_tests.py`,
  `test_maths.py`). Assume no coverage of branching.

---

## Phase 0 — Orientation (no code changes)

Deliverables: `docs/metric_gridding/00_pipeline_notes.md`,
`docs/metric_gridding/DECISIONS.md` (initially empty except conventions).

Tasks:
1. Confirm every repository fact above by reading the code. Correct any that
   are wrong in `00_pipeline_notes.md`.
2. Trace every consumer of leaf column 1 (the half-width): `branch`,
   `transform`, `resolve`, `validate`, `ascend`, `report`, and anything in
   `world_tree.py`, `prune.py`, `io/cands.py`. List them with file:line.
3. Document how `B(s)` is computed for threshold calibration and where
   `eta`, `nbins`, and `tiling_strategy` enter it.
4. Write down the current gridding invariant in one paragraph:
   "For any true λ* inside a leaf's box, `max_t |Φ(t;λ*) − Φ(t;λ)| ≤ η/N_b`
   over the current accumulated interval." Note that this is a sup-norm
   (worst-instant) criterion.
5. Write down the target invariant: "For any true λ* inside a leaf's
   ellipsoid, `m(λ*, λ) = δλᵀ g δλ ≤ m_max`, where `m` is fractional S/N
   loss to second order." Note that this is a mean-square criterion, and
   that the two are not equivalent. Record in `DECISIONS.md` how `m_max`
   will be tied to `η` for the comparison in Phase 3.
6. Record the ordering convention (reverse `d_k`), sign convention of `Δt`
   in `T(Δt)`, and units (`d_k` in m s⁻ᵏ, `f_0` in Hz, `C_VAL`).

Exit criterion: notes reviewed by the human; no open questions about how
the current code represents a leaf.

---

## Phase 1 — The metric, standalone

Deliverable: `src/pyloki/core/metric.py` (pure NumPy; Numba later),
`tests/test_metric.py`.

### API

```python
def poly_phase_metric(
    t_ref: float, t_start: float, t_end: float,
    poly_order: int, f0: float, nbins: int,
) -> np.ndarray:
    """Mismatch tensor g (n_params × n_params) for the kinematic Taylor
    basis [d_kmax, ..., d_2, d_1] about t_ref over [t_start, t_end].

    Φ(t) in cycles: Φ = f0 * [(t - t_ref) - d(t)/c], so
    ∂Φ/∂d_k = -(f0 / c) * (t - t_ref)^k / k!.
    g_ij = <∂_iΦ ∂_jΦ> - <∂_iΦ><∂_jΦ>   (time average over the interval).
    Constant-phase mode is projected out by the subtraction.
    Ordering matches leaf arrays (reverse k).
    Scale so that m = δdᵀ g δd is the fractional S/N loss for a single
    harmonic; document the factor of (2π)^2 and the harmonic weighting
    choice for narrow pulses in DECISIONS.md.
    """

def mismatch(g: np.ndarray, delta: np.ndarray) -> float | np.ndarray

def ellipsoid_axis_extents(g: np.ndarray, m_max: float) -> np.ndarray:
    """Bounding-box half-widths of {δ : δᵀ g δ ≤ m_max} along each axis:
    sqrt(m_max * diag(inv(g))). Directly comparable to the current dparam."""

def transform_metric(g: np.ndarray, t_mat: np.ndarray) -> np.ndarray:
    """g' = T^{-T} g T^{-1} for the coefficient map d' = T d.
    Reuse the T construction from transforms.shift_taylor_errors."""

def cholesky_factor(g: np.ndarray) -> np.ndarray   # g = L Lᵀ

def m_max_from_eta(eta: float, nbins: int, poly_order: int,
                   t_ref, t_start, t_end, f0) -> float:
    """Criterion bridge: choose m_max so the ellipsoid matches the current
    η-box along the best-conditioned axis (or by matched volume — decide and
    record). Used only for before/after comparison, not as a design rule."""
```

### Tests (all must pass before Phase 2)

1. **Brute-force mismatch.** Random λ, small δλ; compute `m` via the
   metric and via direct evaluation: sample both phase polynomials on a
   fine time grid, take variance of the difference, scale identically.
   Agreement to second order (relative error → 0 as |δλ| → 0).
2. **Covariance under epoch shift.** For random Δt,
   `transform_metric(g, T(Δt))` equals `poly_phase_metric` recomputed
   directly about `t_ref + Δt` over the same absolute interval.
3. **Ellipsoid invariance.** Sample points on the boundary
   `δᵀ g δ = m_max`; map them with `T`; check they lie on the boundary of
   the transformed ellipsoid. This is the property that dissolves the
   tiling dilemma; test it explicitly.
4. **Empirical S/N loss.** Using `simulation/`, inject a pulsar with λ,
   fold with λ and with λ + δλ using the existing fold + `Score` path,
   measure `(S/N_mismatched / S/N_true)^2`, compare with `1 − m(δλ)` in the
   small-mismatch regime for several duty cycles. Record the harmonic
   weighting that best matches boxcar scoring in `DECISIONS.md`.
5. **Sanity against current spacing.** For a symmetric interval,
   `ellipsoid_axis_extents` with `m_max_from_eta` should be within a small
   factor of `poly_taylor_step_d_vec` output; report the ratio per axis.
   Large disagreement means a units or ordering bug.

Exit criterion: tests 1–4 green; test 5 ratios understood and documented.

---

## Phase 2 — Metric-aware leaves behind `tiling_strategy="metric"`

Deliverables: changes in `config.py`, `utils/transforms.py`,
`core/taylor.py`, `utils/world_tree.py` (if needed),
`dynamic/dyn_poly_taylor.py` (dispatch only), `tests/test_branch_metric.py`.

### Design decisions to make first (record in DECISIONS.md)

- **Metric storage.** All leaves at a stage share the same interval; `g`
  depends on the leaf only through `f0`. Prefer one Cholesky factor per
  stage per `f0`-bin over one per leaf. Only fall back to per-leaf storage
  `(n_leaves, n_params, n_params)` if the `f0` dependence cannot be binned
  acceptably. Memory is the binding constraint in EP; justify the choice.
- **Lattice.** Start with a hypercubic lattice in whitened coordinates
  (unit spheres of radius 1 covering; lattice spacing 2/√n for covering
  radius 1 in n dims, or simpler: spacing chosen so the sphere covering
  radius equals 1). Leave a hook for A_n* later. Do not optimise the
  lattice until Phase 3 shows the covering redundancy is the bottleneck.
  **Amended 2026-09-10:** the redundancy is already measured, at Phase 2, and
  it is large: 45x the volume lower bound at `poly_order=4` and past a 500k cap
  at `poly_order=5`, against 512 children for `aggressive` on the same step. The
  overhead factorises exactly as `4.935` (covering thickness of `Z^4`) x `9.2`
  (retention dilation). A_4* buys 2.8x of the first factor; anisotropic spacing
  in the metric's eigenbasis buys most of the second, since two of the four
  parent semi-axes are already thinner than one child's covering radius. Note
  ~1000 children at `poly_order=4` is the *volume bound* and so irreducible for
  any complete covering — `aggressive` sits essentially at it. Phase 3 no longer
  needs to establish whether redundancy is the bottleneck; it needs to decide
  whether the metric strategy's detection gain is worth ~2x branching at best.
- **What column 1 means under `metric`.** Either (a) store the ellipsoid
  axis extents in column 1 so downstream readers that expect a half-width
  keep working approximately, or (b) leave column 1 unused and guard every
  consumer found in Phase 0 step 2. Choose (a) unless a consumer needs
  exact box semantics.

### Changes

1. `config.py`: add `"metric"` to `tiling_strategy` choices; add
   `m_max` (or `metric_eta_bridge: bool`) and `metric_lattice: str`.
2. `transforms.py::shift_taylor_errors`: add the `"metric"` branch. Under
   this strategy the function must transform the Cholesky factor (or `g`)
   exactly with `transform_metric` and return the new axis extents for
   column 1. No inflation, no diagonal truncation.

   **DONE 2026-09-14**, but *not* in `shift_taylor_errors` (DECISIONS.md D18).
   As written the step is not implementable — D12 stores only the ellipsoid's bounding
   box in column 1, and a bounding box does not determine the ellipsoid — so the exact
   transform was lifted one level up, to where it is not only possible but cheap:
   `delta_t` and the stage metric are both per-stage, so
   `taylor.metric_transform_extents` computes the re-centred extents once per level at
   `f0 = 1` (same exact `1 / f0` law as the covering, D11) and `transform_extents`
   joins the shared `transform` signature. `poly_taylor_transform_batch` gained the
   `"metric"` branch: values via `shift_taylor_params` exactly as before, column 1 from
   the table, `d_0` left at zero as the projected-out mode. `transforms.py` itself is
   untouched, so the box strategies are bit-identical.

   The other half of this step, not in the original text: `B(s)`.
   `generate_bp_poly_taylor` simulates *box* branching and is `@njit`, so `"metric"`
   gets `generate_bp_poly_taylor_metric`, which is exact rather than averaged (D19).
   Without it the threshold scheme could not be built and nothing could run.
3. `core/taylor.py::poly_taylor_branch_batch` (new function
   `poly_taylor_branch_metric_batch`, dispatched by strategy). **Amended
   2026-09-10:** split into two functions instead of one, see step 6:
   `metric_branch_tables` (plain Python, once per stage) and
   `poly_taylor_branch_metric_apply` (`@njit`, once per batch):
   - Compute `g_s` for the new interval (stage-s metric) and its Cholesky
     factor `L_s`.
   - For each parent: whiten the parent ellipsoid with `L_s`
     (parent metric `g_{s−1}` → `L_s⁻¹ g_{s−1} L_s⁻ᵀ`), lay down the
     lattice, keep lattice points with whitened parent mismatch ≤ `m_max`
     (plus a margin so the parent is fully covered, not just its centres),
     un-whiten to Taylor coordinates.
   - Output the same leaf array shape as the existing function so
     `resolve`, `shift_add`, `score` are unchanged.
   - Preserve `batch_origins` semantics.
4. `core/taylor.py::poly_taylor_resolve_batch`: children no longer sit on
   the rectangular base grid `G0`. Add a diagnostic (behind a debug flag)
   that computes the metric mismatch between each child and the nearest
   base-grid point it resolves to, using the *base-segment* metric.
   Log the distribution. If the p95 exceeds the per-stage budget, stop
   and discuss (options: denser `G0`, or accept and fold into thresholds).

   **DONE 2026-09-14, and the gate TRIPPED. See DECISIONS.md D22.**
   `taylor.metric_resolve_mismatch` + `Pruning._log_resolve_mismatch` behind
   `cfg.metric_resolve_diagnostic`. Measured p95 of **454 / 598 / 5770** against
   `m_max = 0.2` over three levels. The diagnostic is sound — the same code on
   `aggressive` children returns p50 = 7e-10.

   Neither of the plan's two anticipated options was the answer, because the cause was
   not `G0` being too coarse. `metric_branch_tables` inferred the parent's region from
   the accumulated baseline alone and never consulted the region the parent actually
   occupies; on short baselines that `m_max` ellipsoid exceeds the **whole search
   space** by 9.2e7x in jerk and 4.1e6x in accel, so children landed outside
   `param_limits` and `resolve` could only clamp them.

   **FIXED the same day (D23), by carrying the region explicitly** — the human chose
   the principled option of the four offered. A region is a form `{d : d^T A d <= 1}`,
   seeded from the leaf's FFA cell, kept as run state, and advanced per level; plus a
   guard, `region_fits_in_one_child`, that emits a single child when the stage cannot
   resolve anything finer — the exact counterpart of the box strategy's
   `shift_bins < eta`. The gate now passes at the control's own floor: p95
   **7.0e-10 / 7.0e-10 / 2.8e-9**, matching `aggressive` exactly.

   The cost findings are retired, not merely reframed: the recorded 869 / 171 / 93
   children per parent measured a region the parent never had. On the smoke config the
   honest branching factor is **1 at every level** (a 1-second observation cannot
   refine an FFA cell at `m_max = 0.2`); on a 64-segment schedule refinement starts at
   level 9 and runs 9 / 33 / 33 / 29 / 27 / 27. Phase 3 must re-measure.
5. `validate`, `report`, `ascend`, `io/cands.py`: audit each consumer from
   Phase 0 step 2; make `metric` strategy either use axis extents or skip.

   **DONE 2026-09-14 (D27): axis extents everywhere, no code change needed.**
   Of the ten consumers Phase 0 traced, `validate` is a no-op, `ascend` and
   `world_tree.py` never read column 1, `branch` reads the per-stage tables instead,
   and `circular`/`chebyshev` refuse `"metric"` outright (D16). What is left is the
   publishing path, which already receives the region's per-axis full span from
   D12/D18/D24. Verified end to end: both strategies publish comparable spans for the
   same candidate. Two pre-existing caveats are recorded in D27 — `report` combines
   column 1 in quadrature as if the axes were independent, and the published span is a
   bounding box that does not tile.
6. `dynamic/dyn_poly_taylor.py`: dispatch only. No logic here.

   **AMENDED 2026-09-10 — this step was blocked as originally written.**
   `dyn_poly_taylor.py::branch_func` is `@njit(cache=True, fastmath=True)`, so
   it cannot call a covering that needs `eigh`, Cholesky solves and a recursive
   enumeration. As written, step 6 depended on the numba-isation that this plan
   schedules for Phase 4 — an ordering error, since Phase 4 is gated on Phase 3,
   which is gated on this step.

   The fix does **not** require moving Phase 4 forward, because of D11: the
   covering depends only on the *stage* (`t_obs_prev`, `t_obs_cur`, `nbins`,
   `ducy`, `poly_order`, `m_max`) and never on the leaves — the metric is exactly
   proportional to `f0**2`, so a leaf's `f0` enters only as a `1/f0` rescale. So
   split the branch along that seam:

   - `metric_branch_tables(...) -> (offsets_unit, extents_unit)` — plain Python,
     all the linear algebra and the enumeration, **once per stage**.
   - `poly_taylor_branch_metric_apply(leaves_batch, offsets_unit, extents_unit)`
     — `@njit`, a broadcast add and a `1/f0` rescale, **once per batch**.

   Verified: the njit half compiles in nopython mode, is callable from an `@njit`
   caller and from inside a `prange`, and agrees with its `py_func` to ~1 ulp
   (`fastmath` FMA contraction). See `tests/test_branch_metric.py::TestNjitDispatch`.

   **Step 6b DONE 2026-09-14, route (i)** (human's choice; DECISIONS.md D15).
   `Pruning._metric_branch_tables` builds the tables once per level and memoises them
   on `(t_obs_prev, t_obs_cur)`; `branch_offsets` / `branch_extents` were added to the
   shared `branch` signature and threaded through `pruning_iteration_batched`. The
   Chebyshev and circular bases accept and ignore them, and `_setup_pruning` now
   refuses `"metric"` on those bases outright (D16) rather than ignoring it silently.
   Two config knobs came with it (D17): `metric_branch_max` (the box `branch_max` is a
   per-axis width and cannot serve as a total) and `metric_ducy`.
   `aggressive` is unchanged, and a test asserts junk tables cannot perturb it.

### Tests

- Existing strategies unchanged: run any available end-to-end script with
  `aggressive` before and after; outputs byte-identical.
- Coverage: sample random points inside a parent ellipsoid (in the
  stage-(s−1) metric); assert every sample is within `m_max` of at least
  one child (in the stage-s metric). This is the "no corner gaps" test.
- Redundancy: report mean number of children covering a random point.
  Compare against the box-strategy equivalent for the same parent.
- Shape/`batch_origins` contract identical to the existing branch function.

Exit criterion: coverage test green; redundancy known; `aggressive` path
untouched. Steps 1, 2, 3, 4 and 6 are done and `tiling_strategy="metric"` completes a
real `prune_dyp_tree` run over three pruning levels, threshold scheme included;
`aggressive` completes the identical run unchanged. Step 5 (the consumer audit) remains.

Step 4 (D22) found that the coverage test was self-consistent rather than correct — it
asked whether children cover the same over-large region the branch assumed. D23 fixed
that by carrying the region explicitly, and the coverage test now checks the region the
branch was actually given. What remains before Phase 3 is **re-measuring the cost**: the
recorded child counts and the 45x redundancy are retired, and `B(s)` has to be taken
again on a schedule long enough for the guard to release.

---

## Phase 3 — Verification against the paper's diagnostics

Reproduce each experiment first with `aggressive` (validates setup), then
with `metric`. Config for all: 18-min observation, 128 segments, `η = 1.0`,
`N_b = 64`, circular-orbit search with `P_orb^min = T_obs`,
`m_c,max = 10 M⊙`, `m_p,min = 1.2 M⊙`, spin period 7 ms, unless stated.

1. **Figure 7 analogue** — cumulative `∏ B(s)` versus stage. Expectation:
   `metric` lies near the aggressive curve. If it is near the quadrature
   curve, the covering is too redundant; revisit lattice / margin.
2. **Recalibrate thresholds.** `B(s)` has changed. Rerun the Viterbi
   optimisation in `detection/` for the `metric` strategy at
   `P_d = 0.1` before any sensitivity comparison. Save the scheme.
3. **Figure 8 / 11 analogue** — injection-recovery versus injected
   significance for `n_run ∈ {1, …, 32}`, circular-orbit panel. The claim
   under test: the near-threshold rightward shift relative to the
   independent-trial prediction shrinks. Also compare the EP score vs
   Ascend-restored score (top panel) — corner gaps should reduce the
   pre-Ascend degradation.
4. **Figure 12(a) analogue** — `P_d[q]` versus anchor segment at
   `S/N = 15`. Expect fewer/shallower dropouts away from the `ḟ = 0`
   phases. The nodal-phase traps are a different mechanism
   (snap–acceleration ill-conditioning) and are not expected to move.
5. **Cost.** Wall-clock and per-stage candidate counts for both strategies
   at equal recalibrated `P_d`. Report the trade honestly, including any
   extra cost in `branch`.

Deliverable: `docs/metric_gridding/03_results.md` with figures and a
short verdict: proceed / iterate lattice / abandon.

---

## Phase 4 — Only if Phase 3 is positive

- Numba-ise `metric.py`; profile `poly_taylor_branch_metric_batch`.
  **Amended 2026-09-10:** no longer a prerequisite for anything, and mostly
  unnecessary. The per-batch half of the branch is already `@njit`
  (`poly_taylor_branch_metric_apply`); what remains in Python — `eigh`, Cholesky,
  Fincke-Pohst — runs once per stage, not per leaf, so it is not on the hot path.
  Numba-ise it only if profiling shows per-stage setup actually matters.
- Circular-orbit extension: `g` from the exact sinusoidal phase model once
  `T_s ≳ 0.2 P_orb`, matching where the code switches to exact circular
  propagation (`core/circular.py`, `dynamic/dyn_circular_taylor.py`).
- Open a discussion issue on `pravirkr/pyloki` with the Phase 3 results
  before porting to `loki`.

---

## Session log template (append to DECISIONS.md each session)

```
## <date> — Phase N
Done:
Conventions fixed:
Open questions:
Next session starts at:
```
