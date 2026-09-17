# DECISIONS.md — metric-based grid refinement

Conventions and choices, fixed once and referred to thereafter. Append a session log
entry at the end of every session.

## Fixed conventions (Phase 0)

These are **observed from the code**, not chosen, and any metric code must match them.

| # | convention | value | source |
|---|---|---|---|
| C1 | Axis ordering | reverse `k`: `[d_kmax, ..., d_2, d_1]` for the `poly_order` branchable axes | `core/taylor.py:52-64`, `:106-109` |
| C2 | Leaf layout | `(n_leaves, poly_order+2, 2)`; row `[-1]` = `f_0`/basis-flag, `[-2]` = `d_0` (never branched), `[-3]` = `d_1`, `[:-3]` = `[d_kmax..d_2]` | `poly_taylor_seed` |
| C3 | Column meaning | col 0 = value, col 1 = per-axis **full span** (see D24; this row said *half-width* and was wrong) | `psr_utils.branch_param_padded`, `core/taylor.py:109` |
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

  > **PARTLY SUPERSEDED by D20 (2026-09-14).** The table above is on a *one-sided*
  > window, `t_ref = 0` over `[0, T]`, with the epoch at the window's edge. The search
  > actually expands about the **centre** of a symmetric window. Re-measured there, at
  > the same `max|tau|`, the divergence is far milder:
  >
  > | poly_order | `axis_tight` (inscribed) | `volume`-matched |
  > |---|---|---|
  > | 2 | `[1.00, 0.52]` | `[1.57, 0.81]` |
  > | 3 | `[1.00, 0.51, 0.65]` | `[2.43, 1.23, 1.59]` |
  > | 4 | `[1.00, 0.50, 0.89, 0.33]` | `[3.72, 1.87, 3.32, 1.23]` |
  > | 5 | `[0.88, 0.44, 1.00, 0.39, 0.25]` | `[5.62, 2.82, 6.41, 2.53, 1.63]` |
  >
  > The **qualitative** claim survives — the ellipsoid is still ill-conditioned, still
  > worsening with `poly_order`, and a box still cannot be both tight and gap-free. But
  > the headline "**56x** on one axis at `poly_order=5`" was an artefact of the window,
  > not a property of the metric: the honest figure is **6.4x**. The quantitative case
  > for this project is real but roughly an order of magnitude weaker than recorded.
  > Phase 3 must quote the symmetric numbers.

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

## Choices made (Phase 2)

- **D11 — Children are stored as offsets from the parent centre, computed once per
  stage at `f0 = 1` and rescaled by `1 / f0` per leaf.** (Phase 2 decision "metric
  storage".) `poly_phase_metric` is *exactly* proportional to `f0**2` — the phase
  derivative is `-(f0/c)(t - t_ref)^k / k!` and the harmonic weighting depends only on
  `(nbins, ducy)` — verified to machine precision (max relative deviation 2.3e-16, and
  the offsets follow the resulting `1 / f0` law to 7e-14). So no metric needs storing in
  the leaf array at all: the Fincke–Pohst enumeration, much the most expensive part,
  runs **once per stage** rather than once per distinct `f0` in the batch.

- **D12 — Column 1 (`dparam`) holds `ellipsoid_axis_extents(g_child, m_max)`, the
  bounding box of the child's mismatch ellipsoid.** (Phase 2 decision "meaning of
  column 1"; the plan's option (a), which Phase 0 already found cheaper than assumed
  because `world_tree.py` does not interpret column 1.) It is the honest per-axis
  half-width and keeps any consumer of that column reading something meaningful, but it
  is deliberately **not** a covering box: an ellipsoid's bounding box does not tile.
  Nothing in the metric path derives spacing from it.

- **D13 — `poly_taylor_branch_metric_batch` takes `coord_prev` explicitly, and is not
  `@njit`.** The box strategy never needs the parent's interval because a leaf's own
  `dparam` column encodes its spacing; a metric does need it, since per-axis half-widths
  cannot represent the parent ellipsoid's *orientation*. `coord_prev` is already
  threaded to `branch_func` and currently unused, so it was available. The alternative
  considered — squatting on the free leaf slot `[:, -2, 1]` — was rejected because
  `core/chebyshev.py:60` already uses that slot for a velocity width in its own layout,
  so writing a parent interval there would couple the two layouts silently.

- **D14 — The metric branch is split at the stage/batch seam, which unblocks Phase 2
  step 6 without moving Phase 4 forward.** (Plan amended 2026-09-10.)

  The blocker: `dynamic/dyn_poly_taylor.py::branch_func` is
  `@njit(cache=True, fastmath=True)` and cannot call a covering that needs `eigh`,
  Cholesky solves and a recursive enumeration. As the plan was written, step 6 depended
  on the numba-isation scheduled for **Phase 4** — an ordering error, since Phase 4 is
  gated on Phase 3, which is gated on step 6.

  Moving Phase 4 forward is not necessary, because D11 already says the covering depends
  only on the *stage* and never on the leaves. So the branch splits along that seam:

  - `metric_branch_tables(...)` — plain Python, all the linear algebra and the
    enumeration, **once per stage**.
  - `poly_taylor_branch_metric_apply(leaves, offsets_unit, extents_unit)` — `@njit`, a
    broadcast add and a `1 / f0` rescale, **once per batch**.

  Verified rather than assumed: the njit half compiles in nopython mode, is callable
  from an `@njit` caller *and* from inside a `prange` (the context that produced the
  earlier SIGABRT), and agrees with its `py_func` to ~1 ulp — `fastmath` contracts the
  multiply-add, so the two are **not** bit-identical, and the test asserts `rtol=1e-14`
  rather than equality. `TestNjitDispatch` pins all of this.

  Consequence for Phase 4: numba-ising `metric.py` is no longer a prerequisite for
  anything, and is mostly unnecessary — what stays in Python is per-stage, not per-leaf,
  so it is off the hot path.

## Choices made (Phase 2, step 6b)

- **D15 — Step 6b takes route (i): the covering tables are explicit arguments to
  `branch`.** (Human's choice, 2026-09-14, from the two routes this section used to
  list.) `pruning_iteration_batched` now passes `branch_offsets, branch_extents` to
  `prune_funcs.branch(...)`, and the two arrays were added to the shared signature:
  three `branch_func`s, six proxy `.branch()` methods and six `@overload_method` stubs
  across `dyn_poly_taylor.py`, `dyn_poly_cheby.py` and `dyn_circular_taylor.py`. The
  Chebyshev and circular halves accept and ignore them (`ARG001/ARG002` are already
  per-file-ignored there); a comment at each says why.

  Explicit dataflow was preferred over the mutable-structref alternative because the
  table is *stage* state passing through a per-batch call: as a structref field it
  could silently go stale, and nothing in the type system would say so. The cost is
  that two other bases carry a Taylor-only concept in their signature.

  `dyn_poly_taylor.branch_func` is the only one that reads them, under
  `if self.tiling_strategy == "metric"`. It **raises** on an empty table rather than
  branching, because `poly_taylor_branch_metric_apply` with zero children would
  otherwise return an empty leaf array and end the run silently.

- **D16 — `"metric"` is refused on the Chebyshev and circular bases.**
  `Pruning._setup_pruning` raises when `tiling_strategy="metric"` is paired with
  `poly_basis="chebyshev"` or with `prune_poly_order=5` (which is how `config.py`
  selects the circular search). `core/metric.py` builds the mismatch tensor on the
  Taylor kinematic basis (C1/C2) and those two paths have their own leaf layouts, so
  the strategy would otherwise be accepted and then quietly ignored. Refusing is the
  honest behaviour and it is what makes the "unused here" comments in their
  `branch_func`s true rather than aspirational.

- **D17 — Two config knobs that step 6b forced into the open.**

  - `metric_branch_max = 500_000`. The metric cap **cannot** reuse `branch_max`: that
    is a per-axis padding width defaulting to **16**, while the metric cap is a total
    (already noted in `poly_taylor_branch_metric_batch`'s docstring), and 16 as a total
    would raise on even the `poly_order=3` covering of 869 children. A separate field
    is the only correct reading.
  - `metric_ducy = 0.0` meaning "follow `ducy_max`", resolved in `__attrs_post_init__`
    the same way `bseg_brute`/`bseg_ffa` already resolve their `0` defaults.

  The `metric_ducy` default is the **least conservative** choice available and is
  flagged as O7 below, not settled here.

- **D18 — Step 2 takes the fourth option: the re-centred extents are a per-stage
  table, threaded down the seam step 6b already built.** (Human's choice, 2026-09-14,
  with (c) as the stated fallback; the fallback was not needed.)

  The plan asks the transform to "transform `g` exactly with `transform_metric` and
  return the new axis extents for column 1", and as written that is impossible: D12
  stores only `ellipsoid_axis_extents` in column 1, and a bounding box does not
  determine the ellipsoid it bounds. But it *is* possible one level up. `delta_t` and
  the stage metric are both properties of the **stage**, so
  `ellipsoid_axis_extents(transform_metric(g, T(delta_t)), m_max)` is a per-stage
  constant, and — since `g` is proportional to `f0**2`, hence `inv(g)` to `f0**-2` —
  it obeys the same exact `1 / f0` law as the covering offsets (verified: 2e-12).

  So `taylor.metric_transform_extents` computes it once per level at `f0 = 1` next to
  `metric_branch_tables`, `Pruning._metric_stage_tables` returns all three arrays, and
  `transform_extents` joins the shared `transform` signature the way `branch_offsets`
  joined `branch`. The result is the exact transform the plan wanted, with no `eigh`
  inside `@njit`, no per-batch cost, and no new structref state. Option (a) was
  rejected for putting per-stage linear algebra back on the per-batch path; option (b)
  for being inexact when exactness turned out to be affordable.

  Under `"metric"`, `poly_taylor_transform_batch` shifts the *values* exactly as the
  box strategies do (`transforms.shift_taylor_params`) and takes column 1 entirely
  from the table — the old column 1 is not an input. Row `[-2]` (`d_0`) keeps a zero
  half-width: it is the constant-phase mode `poly_phase_metric` projects out, so it
  has none.

- **D19 — `B(s)` under `"metric"` is computed, not simulated.**
  `generate_bp_poly_taylor` models branching as a product of independent per-axis
  counts — exactly the axis-aligned assumption the covering replaces — and is `@njit`,
  so it could not call the enumeration anyway. Under `"metric"` the answer is both
  simpler and *exact rather than averaged*: the covering depends only on the stage
  (D11), so every parent at level `s` emits the same number of children and there is
  no `f0` spread to average over. `core/taylor.py::generate_bp_poly_taylor_metric`
  walks the `MiddleOutScheme` and returns `len(metric_branch_tables(...)[0])` per
  level; `config.generate_branching_pattern` dispatches to it.
  `generate_branching_pattern_approx` **raises** under `"metric"` rather than return a
  worst-case per-axis factor that has no meaning for a covering.

- **D20 — CORRECTION: `coord[1]` is a half-width, not an interval endpoint.**
  (Found 2026-09-14 while scoping step 4, which needs "the base-segment metric".)

  `metric_branch_tables` called `poly_phase_metric(0.0, 0.0, coord[1], ...)`, and
  `metric_transform_extents` inherited it. That is wrong twice over.
  `MiddleOutScheme.get_coord` returns `(ref, scale)` with `ref` the **centre** of the
  accumulated window and `scale = ref - min` its **half-width**; the box code's name
  for it, `t_obs_minus_t_ref`, means *max |tau|*, which is exactly what a sup-norm step
  size needs. The metric caller reused that quantity as an *averaging window*, so it
  put the epoch at the window's **edge** and averaged over **half** its true length.

  The correct windows, relative to the leaf's epoch (the previous centre, since
  `coord_cur = (prev_ref, cur_scale)` under a moving grid):

  | quantity | window |
  |---|---|
  | parent | `[-t_half_prev, +t_half_prev]` (the epoch *is* its centre) |
  | child, while scoring | `[delta_t - t_half_cur, delta_t + t_half_cur]` |
  | child, after the transform | `[-t_half_cur, +t_half_cur]` about the new epoch |

  with `delta_t = coord_next[0] - coord_cur[0]` — already in `_metric_stage_tables`'s
  cache key, so nothing structural changed. `metric_transform_extents` got *simpler*:
  rebuilding `g` about the new epoch over its symmetric window is the transported
  metric, so the `shift_matrix` / `transform_metric` round trip is gone. Phase 1 test 2
  is the statement that the two agree, and `TestTransformExtents` now pins it in the
  exact configuration the loop uses.

  **Why the tests could not catch it.** Phase 1 is innocent: `poly_phase_metric` takes
  `t_ref` and `[t_start, t_end]` independently and `test_metric.py` exercises it with
  the epoch centred. Phase 1 test 5 could not catch it either, because
  `m_max_from_eta` *calibrates* `m_max` against the box, so any global scale error in
  `g` is absorbed by construction. Only the callers were wrong.

  **What it changes, and what it does not.**
  - Extents move by **4x to 64x**, growing with `poly_order`. So the absolute meaning
    of `m_max`, the reported column-1 widths, and D8's ratios all shift (see D8).
  - Child counts barely move — 869 -> 871, 171 -> 171, 85 -> 93 on the smoke config —
    and the **Phase 2 findings table is essentially unchanged** (39 / 871 / 46691
    against volume bounds 8.0 / 63.9 / 1021, i.e. overheads 4.9 / 13.6 / 45.7 versus
    the recorded 4.9 / 13.6 / 45.0). The covering depends on the parent-to-child
    *ratio*, which the window error largely cancels out of. The redundancy analysis,
    its factorisation into lattice thickness x retention dilation, and the conclusion
    that ~1000 children at `poly_order=4` is the intrinsic volume bound all stand.
  - `TestTransformExtents::test_uses_the_symmetric_window_not_the_half_width_as_an_endpoint`
    pins the convention so it cannot quietly revert.

- **D21 — `m_max_from_eta` sized its box from `t_end - t_start`.**
  That is `max|tau|` only when the epoch sits at an endpoint; for the centred window
  the search uses it is 2x too long. Now `max(t_end - t_ref, t_ref - t_start)`, which
  is what `poly_taylor_step_f`'s `(tobs - t_ref)` means. No recorded number changes —
  every existing call has the epoch at an endpoint — but the bridge is now correct for
  the configuration Phase 3 will actually quote.

## Verified end-to-end (2026-09-14)

`tiling_strategy="metric"` now completes a real `prune_dyp_tree` run. On the
`tests/test_prune.py` fixture at `poly_order=3`, `bseg_ffa = nsamps // 4` (4 segments,
3 pruning levels), the branching pattern builds and the run finishes through ascend and
report, with three distinct per-stage coverings:

| level | `t_obs` prev -> cur (s) | children per parent |
|---|---|---|
| 1 | 0.131 -> 0.262 | 871 |
| 2 | 0.262 -> 0.393 | 171 |
| 3 | 0.393 -> 0.524 | 93 |

The counts fall with level because the covering depends on the *ratio* of the
accumulated intervals, which drops from 2.0 to 1.33 as the baseline grows.
`aggressive` completes the identical run, unchanged.

- **D22 — STEP 4 RESULT: the covering has no idea where its parent actually is, and
  the plan's gate has tripped. STOP AND DISCUSS.** (2026-09-14.)

  Step 4 asked for the mismatch between each child and the base-grid point it resolves
  to, with the instruction "if the p95 exceeds the per-stage budget, stop and discuss".
  Measured on the 3-level smoke run at `m_max = 0.2`:

  | level | children | p50 | p95 | max | verdict |
  |---|---|---|---|---|---|
  | 1 | 7839 | 76.2 | 454 | 688 | **2270x over** |
  | 2 | 10944 | 99.6 | 598 | 1180 | **2990x over** |
  | 3 | 5952 | 315 | 5770 | 6664 | **28900x over** |

  **The diagnostic is not at fault**; the control settles that. Run on `aggressive`
  children over the identical step, the same code returns **p50 = 7.0e-10**, as it
  should: box children sit on a rectangular refinement of `G0` by construction. Two
  unit tests pin the forward map independently -- a child placed exactly on a cell
  centre costs `< 1e-18`, and a half-cell offset matches a hand computation of
  `g[1,1] * dv**2` to 1e-5.

  **Root cause, and it is not rounding.** `metric_branch_tables` defines the parent's
  region as `{d : d^T g_parent d <= m_max}`, built from the accumulated baseline alone.
  Nothing in that expression refers to the region the parent *actually* occupies. On
  the short baselines of the early levels the two are wildly different. At
  `t_half = 0.131 s`, against the search space:

  | axis | `m_max` ellipsoid / search half-span |
  |---|---|
  | `d_3` (jerk) | 9.2e7 |
  | `d_2` (accel) | 4.1e6 |
  | `d_1` (freq) | 0.66 |

  A 0.13 s baseline cannot constrain jerk or acceleration at all, so the `m_max`
  ellipsoid on those axes runs millions of times past the entire search range. The
  covering dutifully tiles it. Children then land at 82-204 Hz against a 141.9-143.9 Hz
  band and at accelerations of +/-9.6e8 against a range of +/-4; `resolve` clamps them
  to the grid edge, and the clamped distance is what the table above is measuring.

  **What is NOT wrong.** The metric itself is sound: on the same baseline the ellipsoid
  tracks the eta-box to within a factor of a few (`< 100x`, and about 4-11x in
  practice), consistent with D8-as-corrected. The covering, the retention bound, the
  lattice and the `1/f0` law are all fine. The defect is one of *scope*: the metric
  branch sizes leaves **absolutely**, from `m_max` and the baseline, where the box
  strategy only ever **subdivides** the region it was given -- and keeps a single child
  when `shift_bins < eta` says refinement is not yet warranted. The metric branch has
  no equivalent of either the starting region or that guard.

  This is why the Phase 2 coverage tests pass and missed it: they ask whether the
  children cover *the stage-(s-1) `m_max` ellipsoid*, which is the same wrong region
  the branch assumes. Self-consistent, and not the question.

  **It also reframes the Phase 2 cost findings.** The recorded 869 / 171 / 93 children
  per parent, and the "45x the volume bound" redundancy, are the cost of covering a
  region the parent never occupied. The true branching factor cannot be known until the
  scope is fixed; it should fall, possibly a great deal, since most of that volume is
  outside the search space.

- **D23 — The parent region is carried explicitly, as run state.** (Human's choice,
  2026-09-14: option (3), the principled one, from the four D22 offered.)

  A region is represented as a form `A` with `{d : d^T A d <= 1}`. The recursion:

  | | region |
  |---|---|
  | seed | the FFA cell it came from, from leaf column 1 (`taylor.metric_seed_region`) |
  | level that refines | the stage's mismatch ellipsoid, `g_s / m_max` |
  | level that cannot refine | **unchanged** — the parent's region is inherited |
  | after the transform | `T^-T A T^-1`, the same map a metric takes |

  `Pruning._region_form` holds it, seeded in `initialize` and advanced once per level in
  `_metric_stage_tables`. No per-leaf storage: it is one region per stage, which is what
  D11 already established for everything else on this path.

  **The guard is the other half, and it is the exact counterpart of `shift_bins < eta`.**
  `metric.region_fits_in_one_child` asks whether the parent already sits inside a single
  child's ellipsoid; in the child's whitened coordinates that is exactly
  `lambda_min(a_mat) >= 1`. When it holds the branch emits **one** child and the region
  is untouched. Its absence is what let the covering tile regions no leaf occupied.

  **Result: the step 4 gate now passes, at the control's own value.** p95 goes from
  454 / 598 / 5770 to **7.0e-10 / 7.0e-10 / 2.8e-9** against `m_max = 0.2` — identical
  to what `aggressive` scores on the same step, which is the floor the diagnostic can
  return. Children no longer leave `param_limits`.

  **What the fix reveals about cost.** On the 1-second smoke config the metric strategy
  now emits **one child at every level**: over a 0.26-1.05 s baseline, `m_max = 0.2` is
  looser than a single FFA cell, so no refinement is possible and the honest branching
  factor is 1. That is correct, not inert — the guard releases as soon as the baseline
  can resolve something. On a 64-segment schedule refinement starts at **level 9** and
  then runs 9 / 33 / 33 / 29 / 27 / 27 children per level; with `m_max = 1e-6` on the
  short config it starts at level 3. Both are pinned.

  **This retires the Phase 2 redundancy findings.** The recorded 869 / 171 / 93 children
  per parent, and the "45x the volume bound", were the cost of covering the over-large
  region. Real branching factors must be re-measured on a schedule long enough to
  refine, and Phase 3 cannot quote the old ones.

- **D24 — CORRECTION: leaf column 1 is a full span, not a half-width.** C3 said
  half-width; `psr_utils.branch_param_padded` reads it as `param_cur -/+ dparam / 2`,
  and a seed's jerk entry is 16.0 for a search range of `[-8, 8]`. So it is the full
  cell span, and C3 above is corrected.

  D12 and D18 had been writing `ellipsoid_axis_extents` — a half-width — into it, a
  factor of 2 too small. The metric path now writes `2 * region_axis_extents`, and
  `metric_seed_region` halves column 1 when reading a cell back. Inside the search loop
  this changed nothing (column 1 is write-only there under `"metric"`), but it was
  wrong in everything `report` and `io/cands.py` print, and it now feeds the region
  recursion, where a factor of 2 is not cosmetic.

- **D25 — A child's region is parent ∩ ellipsoid, not the ellipsoid.** Carrying
  `g_child / m_max` alone let the region *grow* on every axis the stage did not
  actually refine; the ellipsoid then shrinks only slowly, so the containment guard
  never fired again and the branch re-tiled at every level.
  `metric.region_intersection` returns the minimum-volume member of `t*A + (1-t)*B`,
  which is a sound outer bound for every `t` and can only beat whichever input is
  smaller. Worth **4.3e33 -> 1.1e26** on `prod B(s)`. Per-axis the bound can still
  exceed the narrower input — unavoidable for an ellipsoidal approximation of a lens,
  and harmless, since the covering consumes the form and not its bounding box.

## Branching-factor comparison (2026-09-14) — the Figure 7 analogue

> **READ THE SECOND HALF FIRST.** The comparison below was initially run on one small
> config and reported as "the metric strategy is not competitive". **That conclusion was
> wrong** — it generalised from the single regime where the box strategy is at its best
> and where there is no problem to solve. Corrected under "Regime dependence" below.

Config: `T_obs = 16.8 s`, 32 segments, `poly_order = 3`, `eta = 1`, `nbins = 64`,
`ducy_max = 0.2`, `m_max = 0.2`. Config-only, no FFA run needed.

| strategy | `prod B(s)` | first refining level |
|---|---|---|
| aggressive / quadrature / conservative | **27** | 2 |
| aggressive, un-coarsened (D5 like-for-like) | 27 | 2 |
| metric, as implemented | **1.1e26** | 8 |

**The metric strategy is not competitive as designed, by 24 orders of magnitude, and
the cause is not the criterion.** Three measurements pin it:

1. **The box is already inside the metric budget.** Its leaf carries mismatch
   **0.002-0.024** against `m_max = 0.2` at every level — 1% to 12% of budget. So the
   metric criterion demands *less* total refinement than the box already performs. The
   entire excess is covering waste.
2. **The overhead is paid per level and compounds.** From the first refining level the
   volume ratio between the parent region and the child ellipsoid is only **1.5-1.9**,
   yet the covering emits **9-33 children**. Covering an ellipsoid with ellipsoids of
   nearly its own size costs the lattice thickness (4.9 at `n = 3`) whatever the volume
   gain is. Over ~24 levels that is `~10^26`. The box never pays it: it splits per axis
   in integer factors, and `shift_bins < eta` defers until a whole factor accrues.
3. **The floor for any covering scheme** is `V_total * thickness^(branch events)`. To
   be competitive it must branch **rarely** — which is in direct tension with pruning,
   since pruning needs branching early to have anything to prune against.

**Deferral rescues it.** Branch only when the parent region overhangs a child ellipsoid
by at least a factor `R` in the worst direction (`1 / sqrt(lambda_min(a_mat))`, the
natural counterpart of `shift_bins < eta`):

| defer `R` | `prod B(s)` | branch events |
|---|---|---|
| 1.0 (containment only — current) | 1.1e26 | 24 |
| 1.5 | 1.0e21 | 19 |
| 2.0 | 4.4e4 | 4 |
| **3.0** | **21** | **1** |
| 4.0+ | 1 | 0 (never refines) |

At `R = 3` the metric strategy costs **21 against the box's 27** — competitive, and
slightly cheaper. Note how sharp the cliff is: `R = 2` is 2000x worse and `R = 4` stops
refining altogether, so `R` is not a knob that can be set carelessly.

**The catch Phase 3 has to settle.** `prod B(s) = 21` at `R = 3` comes from a *single*
branch event: the tree stays one leaf wide until that level, then fans out 21x. The box
instead branches 3x at level 2 and 3x at level 6. Cost parity therefore says nothing
about detection parity — pruning works by having diversity to threshold against, and a
single late fan-out may well detect worse. That is an injection-recovery question, and
it is the real Phase 3 experiment.

Deferral is **not implemented**; the table above is a simulation over the same schedule.

- **D26 — `metric_defer_factor`: branch only once the region overhangs by `R`.**
  Generalises the D23 guard from exact containment (`lambda_min >= 1`) to
  `region_overhang <= R`, where the overhang `1 / sqrt(lambda_min(a_mat))` is how far
  the parent reaches past one child in the worst direction — the direct counterpart of
  the box's `shift_bins < eta`.

  **Default 1.0, i.e. exact containment**, so the shipped behaviour keeps the `m_max`
  guarantee for every leaf at every level. `R > 1` trades it for cost: a leaf may then
  carry up to `R**2 * m_max` of mismatch between branch events. That price is steep and
  easy to overlook — at `R = 3, m_max = 0.2` a leaf reaches `m = 1.8`, outside the
  regime where the quadratic mismatch model is even valid — so it is not defaulted on.

  Deferral is nonetheless **necessary**, not optional: re-covering costs the lattice
  thickness whatever the volume gain, so at `R = 1` the overhead compounds over levels.

## Regime dependence — the corrected comparison (superseding two earlier attempts)

> **This section has been wrong twice.** First it reported "the metric strategy is not
> competitive" from a single small config. Then the correction over-swung: it compared
> the metric against a baseline computed with **`conservative`** error propagation while
> labelling it "the box", which flattered the metric by up to 6 orders of magnitude.
> `conservative` is a *different, more expensive strategy*, not an honest accounting of
> `aggressive`. The table below is the third and, as far as I can check, correct
> version: the hand-rolled box loop now reproduces `generate_branching_pattern` exactly
> for `aggressive` (5.314e5 both ways at 67 s / `po=4`), which is the validation the
> earlier attempts lacked.

`prod B(s)` at **equal worst-case leaf mismatch**, against **both** box strategies. The
metric column is optimised over the deferral factor `R` to reach `conservative`'s
mismatch.

| config | `aggressive` | `conservative` | metric (at `conservative`'s `m`) |
|---|---|---|---|
| 17 s, 32 seg, `po=3` | 27 @ m=0.024 | 27 @ m=0.024 | 147 (`R=16`) |
| 67 s, 64 seg, `po=3` | 243 @ m=0.038 | 1.77e5 @ m=0.047 | 5.05e4 (`R=12`) |
| 67 s, 64 seg, `po=4` | 5.31e5 @ m=0.291 | 7.41e15 @ m=0.368 | 1.17e9 (`R=8`) |
| 268 s, 64 seg, `po=4` | 1.51e12 @ m=0.521 | 1.48e32 @ m=0.533 | 1.35e19 (`R=16`) |

**This is the tiling dilemma, quantified — and the metric lands between the horns.**
`aggressive` is cheap and gappy; `conservative` is gap-free and, by `po=4`, absurd
(`1e32`). They diverge by ten to twenty orders of magnitude, and that gap *is* the
problem the project exists to attack. The metric sits consistently about **a third of
the way** from `aggressive` to `conservative` on a log scale — 33% at 67 s/`po=4`, 34%
at 268 s/`po=4`.

So the honest headline is neither "the metric wins" nor "the metric loses":

- against `conservative`, the only gap-free box option, the metric is **6.3e6 x cheaper**
  at 67 s/`po=4` and **1.1e13 x cheaper** at 268 s/`po=4`;
- against `aggressive`, it is **2200x** and **8900x** *dearer* at the same two points.

Whether that trade is worth taking depends entirely on **how much detection
`aggressive` actually loses to its gaps** — which no branching-factor calculation can
answer. That is Phase 3 step 3, and it is now unambiguously the decisive experiment.

Two things the table does not say. These are branching factors, not detections: EP caps
candidates and thresholds hard, so `1e32` is never realised and `conservative`'s real
failure is that it is unrunnable, while `aggressive`'s is lost sensitivity. And equal
*worst-case* mismatch at m ~ 0.5 is a poor operating point for anybody; the comparison
is a like-for-like, not a recommendation.

## Phase 3 reduced-grid validation (2026-09-14)

Run before committing to a full injection grid, and it changed the config twice.

- **D30 — Phase 3 runs at 268 s / 64 segments / `poly_order=4`, chosen by measurement.**

  | config | FFA fold array | base grid `[snap, jerk, accel, freq]` |
  |---|---|---|
  | 67 s, 64 seg | 0.01 GB | `[1, 1, 1, 12]` — **degenerate** |
  | 268 s, 64 seg | 0.1 GB | `[1, 1, 5, 467]` |
  | 537 s, 128 seg | 5.0 GB | `[1, 1, 22, 3383]` |
  | **1074 s, 128 seg (the plan's)** | **3778 GB** | `[1, 2, 526, 53151]` |

  Two hard facts. The plan's 18-minute config needs a **3.8 TB** fold array and is not
  runnable on this machine; 537 s / 128 segments needs 5 GB and is a spot check at
  best. And below ~268 s the base grid is degenerate — at 67 s only frequency is
  gridded, so all three strategies recover the *same* candidate with an identical
  mismatch of 8.716e-4 and nothing is discriminated. My earlier choice of 67 s was
  therefore wrong, and any 67 s number in this document should be read with that in
  mind. 268 s is the smallest config that both fits and discriminates.

- **D31 — the metric's apparent recovery failure is thresholding, not coverage.**
  With a *shared* placeholder threshold ramp (1.5 → 8.0) and `max_sugg = 2**11`:

  | strategy | recovered | best mismatch | best score | candidates |
  |---|---|---|---|---|
  | `aggressive` | yes | 0.0015 | 19.5 | 1090 |
  | `conservative` | yes | 0.167 | 13.5 | 1730 |
  | `metric` | **no** | 2.8e6 | 4.9 | 217 |

  Re-run with zero thresholds and `max_sugg = 2**14`, the metric recovers cleanly —
  best mismatch **0.0037** at score **19.0**, against `aggressive`'s 0.024 at 17.3 on
  the same setting. So the covering contains the injected signal; the true track was
  simply thresholded away.

  That is exactly the failure the plan's ordering guards against: "Recalibrate
  thresholds ... **before any sensitivity comparison**", because `B(s)` differs between
  strategies and a scheme tuned for one is meaningless for another. It is now
  demonstrated rather than assumed, and no sensitivity number should be quoted until
  step 2 is done per strategy.

  `conservative`'s degraded 0.167 / 13.5 is the same effect: its `B(s)` is enormous, so
  it overruns `max_sugg` and thresholds hard.

Harness: `docs/metric_gridding/injection_recovery.py`. Recovery is judged by mismatch to
the truth in the full-baseline metric, not by parameter distance — a tolerance in Hz
would favour whichever grid happens to be finer in frequency.

## Phase 3 step 2 — threshold recalibration (2026-09-14)

- **D32 — ~~BLOCKER~~ RESOLVED 2026-09-14: the Viterbi threshold optimiser was broken
  on this base.** Two independent upstream bugs, both now fixed, tested and filed:
  **issue #8** (zero-filled state records) and **issue #9** (`trials_scheme` returning
  `-inf`). Branch `fix-viterbi-zero-states`, cut from PR #7 since the crash it fixes
  masks both. Details of the original diagnosis below.
  `DynamicThresholdScheme.run()` completes without error but leaves **every state empty
  from stage 2 onward**, so `backtrack_best` raises "No non-empty final states to plot"
  and no scheme can be extracted. Stage 0 has 36 non-empty states, stage 1 has ~155,
  stage 2 has 0.

  Not our doing, and not the config:
  - it reproduces on the **example notebook's own branching pattern**
    (`examples/optimal_thresholds.ipynb`, 127 stages, max `B = 8`), with the notebook's
    own settings, so it is nothing to do with the metric strategy;
  - it happens in **both** `legacy` and `improved` modes;
  - it is not the PR #7 cherry-pick (D28): `boxcar_snr_2d_serial` is **bit-identical**
    to `boxcar_snr_2d` on the same input, verified directly;
  - the beam is not starving it — stage 2 offers 42 candidate thresholds.

  So the failure is in the state propagation between stages 1 and 2. Diagnosing it is
  work in `detection/`, unrelated to metric gridding, and probably belongs upstream as a
  sibling to PR #7.

- **D33 — Fallback: per-strategy schemes from `determine_scheme`, not optimised.**
  `thresholding.determine_scheme(probs, bp, ...)` with `probs = 1 / B(s)` ("expect one
  survivor per branch") works, and is per-strategy, which is the property D31 says is
  indispensable. It is **not** P_d-optimal, so every number below is provisional.

  | strategy | `prod B(s)` | `P_d` | log2 complexity | **log2 cost = complexity / P_d** |
  |---|---|---|---|---|
  | `aggressive` | 1.51e12 | 0.040 | 7.28 | **11.93** |
  | `conservative` | 5.29e32 | 1.3e-6 | 8.81 | **28.36** |
  | `metric` | 1.35e19 | 0.141 | 14.10 | **16.92** |

  First evidence on the axis that actually matters. The metric reaches **3.5x the
  detection probability of `aggressive`** (0.141 vs 0.040) because it branches at only
  5 of 63 stages and can therefore afford much lower thresholds (1.59-4.02, against
  `aggressive`'s 2.11-6.83). It pays for that in complexity, and on the combined
  figure of merit it lands at **2^5 = 32x worse than `aggressive`** and **2^11 =
  87 000x better than `conservative`** — the same "between the horns, nearer the cheap
  one" position the branching comparison found.

  **Caveats.** `P_d` is noisy: the scheme's RNG is unseeded and `ntrials = 1024`, and
  repeat runs gave 0.040-0.25 for the same strategy. These are single draws, the
  schemes are not optimised, and the ranking of `aggressive` vs `metric` on cost is
  well inside neither. Treat the table as an order-of-magnitude sighting shot, not a
  result. The real comparison needs D32 fixed, several seeds, and step 3.

- **D34 — The second bug, found only after fixing the first.**
  `schemes.trials_scheme` is `norm.isf(1 / cumprod(B))`. A branching pattern that starts
  with `B(s) = 1` has a cumulative trial count of 1, so the leading entries are
  `norm.isf(1) = -inf`. `DynamicThresholdScheme` centres its threshold beam on that
  path, so the beam selects nothing and every state is empty from stage 0 — the same
  end-user symptom as #8, a completely different cause, and it survives the #8 fix.

  This is exactly the case the metric strategy produces and the box strategies do not:
  `aggressive` starts at `B = 4`, while the metric starts with four unbranched stages
  (D26's guard doing its job). Fixed by flooring the path at zero — no branching means
  no trials pressure, so no threshold is required — which leaves already-branching
  patterns numerically unchanged.

## Phase 3 step 2 — COMPLETE (2026-09-14)

Viterbi-optimised scheme per strategy, `P_d` target 0.1, on the 268 s config:

| strategy | `prod B(s)` | `P_d` achieved | log2 complexity | log2 cost | thresholds |
|---|---|---|---|---|---|
| `aggressive` | 1.51e12 | 0.1031 | **8.16** | 11.44 | 2.30 – 7.70 |
| `metric` | 1.35e19 | 0.1031 | **18.22** | 21.50 | 1.40 – 8.40 |
| `conservative` | 5.29e32 | 0.1031 | **49.99** | 53.27 | 1.70 – 9.10 |

**All three now hit the same detection probability**, so for the first time the
comparison is at equal `P_d` rather than equal worst-case mismatch — which is what the
plan asked for and what makes the complexity column meaningful.

At `P_d = 0.103` the metric costs `2^10.1` = **1070x** `aggressive` and `2^31.8` =
**3.7e9 x less** than `conservative`. On a log scale it again sits about a quarter to a
third of the way from the cheap-and-gappy strategy to the safe-and-unusable one — the
same position three independent measurements have now found.

Note this is the metric at its *cost-optimal* operating point, `(m_max, R) =
(0.00208, 16)`, chosen to match `conservative`'s worst-case mismatch. Other operating
points move that row and nothing else; settling `(m_max, R)` jointly with O6 and O7 is
still open.

These supersede the provisional D33 numbers, which came from the non-optimised
`determine_scheme` fallback and had `P_d` varying by 6x between runs.

## Reading the paper (2026-09-14) — a correction and two findings

The plan's ground rules say "Verify citations before relying on them." I did not, until
now, and asserted that D36 "contradicts Kumar & Zackay (2026) §5.2.4" on the strength of
the plan's one-line paraphrase. The source is in `paper/pruning1.tex`. Having read it:

- **D37 — RETRACTION: D36 does not contradict §5.2.4. It sharpens one sentence of it
  and confirms another.**

  §5.2.2 (the tiling trilemma) describes the gaps as *per-tile* — `aggressive`'s inner
  AABB misses its own sheared corners, which is exactly what Figure `grid_tiling`
  panel (d) shows, and is true. It then says: *"These gaps are partially mitigated in
  practice by **natural overlap between neighbouring templates**."* D36 measures that
  collective cover and finds the mitigation is not partial but **exact**, for a
  lower-triangular unit-diagonal `T`. That is a strengthening of the paper's own
  caveat, not a contradiction.

  §5.2.4 then anticipates this project's cost result outright:

  > "An alternative approach ... is transporting a non-orthogonal lattice that **exactly
  > tracks the coordinate shear** ... To preserve the `eta` bound, the sheared cell must
  > be subdivided along its elongated axes. The resulting template count therefore
  > approaches the same scaling as a conservative AABB cover. Consequently, **tracking
  > exact geometric shear provides no practical computational advantage** over bounding
  > box methods when a rigid `eta` constraint is enforced."

  The measured 1070x is that prediction coming true. The paper also adopts *aggressive
  Taylor as the operational default* and says its gaps "must then be controlled
  empirically, for example by tightening `eta`" — the cheaper alternative.

  The closing remark the plan is built on — *"A complete solution likely requires
  replacing fixed coordinate spacings with a local metric-based mismatch criterion.
  **Further work is required to quantify the sensitivity loss**"* — is a hedged
  speculation and a call to measure, not a claim that a metric would be cheaper. The
  plan's Purpose reads it as "the metric dissolves the dilemma"; §5.2.4 argues the
  opposite.

- **D38 — §3.4: the box grid is already metric-derived, which reframes the whole
  comparison.** The paper states the problem this project set out to solve, in §3.4:

  > "the valid parameter search volume is a highly elongated hyper-ellipsoid (a
  > 'needle') rather than a hyper-rectangle. A simple rectangular grid ... is therefore
  > highly redundant, as the grid axes do not align with the principal axes of the
  > **parameter metric**."

  And it already fixes it analytically: *"we retain the physically intuitive Taylor
  coefficients for the search coordinates but define the grid density based on an
  orthogonal basis analysis ... utilizes Chebyshev polynomials to **diagonalize the
  parameter metric**"*. The output is exactly the `2**(k-1)` coarsening factor
  (Appendix `app:optimal_gridding`).

  So the shipped `aggressive` grid is **not naive** — it already carries a
  metric-diagonalisation correction. D5 decided `"metric"` would not inherit that
  coarsening, so our covering has been competing against a box that is metric-corrected
  on the diagonal while ours pays full price. What a metric covering can still add over
  that is only the **off-diagonal** correlation structure, which is a much smaller prize
  than the plan assumed — and D8-corrected already measured that residual anisotropy at
  1.6x-6.4x, not the 56x originally recorded.

  §5.2.3 adds that Chebyshev alone buys ~10x fewer templates from isotropy, and that
  "applying the coarsening factor brings the Taylor branching pattern down to a similar
  scaling as the relevant Chebyshev profiles" — i.e. the coarsening already captures
  most of the orthogonal-basis gain.

**Bottom line.** The paper anticipated the cost result, already applies a diagonal
metric correction, recommends the strategy our measurements favour, and asks for the one
thing nobody has produced: a *quantified sensitivity loss*. That last item is the only
part of the plan's premise still standing, and it is not what this branch has been
building.

## Still to do in Phase 2

- **Step 5** — the consumer audit. Phase 0 traced ten consumers of column 1; the ones
  that still matter under `"metric"` are `poly_taylor_report_batch`,
  `periodogram.add_run` and `io/cands.py`, all of which turn column 1 into *reported*
  uncertainties. D24 means those were a factor of 2 out until now. Worth confirming that
  an ellipsoid bounding box is what they should print.

  Note that column 1 is no longer write-only under `"metric"`: D23 reads a seed's cell
  back out of it, so it is load-bearing at the start of every run.

- **Choose the operating point.** `metric_defer_factor` exists (D26) but defaults to
  1.0, which is correct but expensive. The winning `R` is regime-dependent (8-12 in the
  target regime) and trades the `m_max` guarantee for cost at `R**2`. Phase 3 has to
  settle `(m_max, R)` jointly, together with O6 and O7 — all four are scale factors on
  the same criterion.

- **An upstream oddity found while wiring step 2, deliberately not fixed.**
  `report_func` (both `dyn_poly_taylor.py` and `dyn_circular_taylor.py`) calls
  `poly_taylor_transform_batch` under `not use_moving_grid` and **discards the result**,
  reporting `leaves_batch` unshifted — the call is already a no-op on `main`. `"metric"`
  skips it rather than raise for a per-stage table that would be thrown away. Fixing the
  discard is an upstream change, out of scope for this branch.

## Phase 2 findings

- **The retention bound was 167x loose; fixed, and the fix is the interesting part.**
  The first implementation (`_retention_radius`) dilated the parent ellipsoid
  *uniformly* by `1 + sqrt(lambda_max(A))`. That equals `1 + r / a_min`, i.e. it is
  tight for the **shortest** semi-axis and scales every longer one by the same factor —
  the bounding-box mistake in spherical form, and with D8's axis ratios (up to 335 at
  `poly_order=4`) it cost 167x the volume bound.

  The obvious tightening — dilate each semi-axis `a_i -> a_i + r` — is **unsound**.
  Containment of `E + B(r)` needs `|v| sqrt(sum a_i^2 v_i^2) <= sum a_i v_i^2`, which is
  Cauchy–Schwarz the wrong way round; sampling `E + B(r)` puts points at **1.40** in
  that form. `tests/test_branch_metric.py::test_per_axis_dilation_would_be_unsound`
  pins that so the argument is not quietly re-derived wrongly.

  What is used instead is the standard S-procedure external ellipsoid of a Minkowski
  sum: `c_i**2 = a_i**2 / t + r**2 / (1 - t)`, sound for every `t in (0,1)`, exactly
  `c = a + r` when the semi-axes are equal, with `t` chosen to minimise the volume.
  Child counts fell **3.7x** (170235 -> 45909 at `poly_order=4`) with coverage and worst
  mismatch **bit-identical** — the discarded points were provably never nearest.

- **Child counts, and the whole overhead accounted for.** `m_max=0.2`, one segment
  doubling 16.78 s -> 33.55 s:

  | `poly_order` | volume bound | metric children | overhead | box strategy (`aggressive`) |
  |---|---|---|---|---|
  | 2 | 8.0 | 39 | 4.9 | 4 = `[4,1]` |
  | 3 | 63.9 | 869 | 13.6 | 32 = `[8,4,1]` |
  | 4 | 1021 | 45909 | 45.0 | 512 = `[16,8,4,1]` |
  | 5 | — | >500000 (cap) | — | — |

  The overhead factorises exactly: `1021 x 4.935 x 9.2 = 46300` vs measured **45909**,
  where 4.935 is the **covering thickness of `Z^4`** (`V_n (sqrt n / 2)^n`, an intrinsic
  property of the cubic lattice, independent of the region — and it shows up directly as
  the measured mean coverings-per-point, 4.94) and 9.2 is the residual retention
  dilation. Nothing is unexplained.

- **CORRECTED: the box strategy branches ~512 per parent at `poly_order=4`, not
  27–81.** An earlier session reported 27–81 and I carried that forward when first
  reading the 170235 as a catastrophe. Measured directly on the identical parent→child
  step, `poly_taylor_branch_batch` emits 4 / 32 / 512 at `poly_order` 2 / 3 / 4 — i.e.
  essentially **at** the volume bound (8 / 64 / 1021), short of it only by a factor ~2
  on the last axis, which `shift_bins < eta` suppresses. Pinned in
  `TestAggressiveUntouched::test_box_branch_counts_unchanged`.

  This reframes the redundancy question. **~1000 children per parent at
  `poly_order=4` is intrinsic**, not a defect of the metric approach: no covering of any
  kind can beat the volume ratio, so the metric strategy is *at best* ~2x the box
  strategy's branching, and the box strategy is not leaving free headroom on the table.
  The real gap is the 45x, and it is addressable rather than fundamental:

  1. **Use a thin lattice.** `A_4*` has covering thickness 1.766 against the cubic
     lattice's 4.935 — a 2.8x saving, and the plan already anticipated this for Phase 3.
  2. **Anisotropic spacing in the metric's eigenbasis.** At `poly_order=4` the parent's
     semi-axes in units of the covering radius are `[237.8, 6.43, 0.94, 0.71]`: **two of
     four axes are thinner than a single child's covering radius**, so no refinement is
     needed along them, yet an isotropic cubic lattice still spends ~5 points on each.
     Per-axis spacing — what the box strategy does in Taylor axes — would recover most
     of the 9.2 dilation factor. This is likely the larger win of the two.

  The measurement is therefore available *now*, at Phase 2, rather than at the Phase 3
  point where the plan intended to gate the lattice decision.

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

- **O7 — Which duty cycle should set the harmonic weighting? (new, from D17.)**
  `cfg.metric_ducy` defaults to `ducy_max`, the **widest** pulse the boxcar bank scores,
  and that is the least conservative reading available. The harmonic weight rises
  steeply as the pulse narrows (D10: 1.59 at `ducy=0.5`, 102 at `ducy=0.05`), so a
  signal narrower than `ducy_max` has its S/N loss under-predicted and is given a
  coarser grid than it deserves — precisely the failure this project exists to remove.
  The safe reading is the *narrowest* duty cycle the search is sensitive to, which is
  set by the smallest boxcar width (~`1/nbins`) and would cost a large factor in
  children per parent.

  Defaulted to `ducy_max` rather than decided, because the honest answer needs the same
  injection-recovery measurement as O6 and the two interact (both are just scale factors
  on `g`). Phase 3 should settle them together. Cheap to change: one config field.

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

## 2026-09-10 — Phase 2
Done:
  - config.py: tiling_strategy gains "metric"; m_max=0.2, metric_lattice="cubic".
    Defaults unchanged, so the shipped path is untouched (test asserts this).
  - core/metric.py: cubic_lattice_spacing, _retention_form,
    _enumerate_lattice_in_ellipsoid (Fincke-Pohst), lattice_children.
  - core/taylor.py: poly_taylor_branch_metric_batch, same contract as
    poly_taylor_branch_batch. Not @njit (see the Phase 2 blocker).
  - tests/test_branch_metric.py: 34 tests. Full suite 113 passed.
Decisions fixed: D11 (offsets, one enumeration per stage via the exact f0^2
  law), D12 (column 1 = child ellipsoid bounding box), D13 (coord_prev passed
  explicitly, not the free leaf slot).
Exit criterion: coverage green (0 uncovered at poly_order 2-4, worst mismatch
  0.90-0.98 of m_max, guaranteed by construction not by sampling); redundancy
  known and factorised exactly into lattice thickness x retention dilation;
  aggressive path untouched.
Phase 2 step 6 unblocked (D14, plan amended): branch_func being @njit does
  NOT require moving Phase 4 forward. Split the branch at the stage/batch seam
  -- metric_branch_tables (Python, per stage) + poly_taylor_branch_metric_apply
  (@njit, per batch) -- since D11 already showed the covering depends only on
  the stage. Verified callable from @njit and from inside a prange. Step 6b
  (the per-stage precompute hook in prune.py) is still to do, so
  tiling_strategy="metric" remains inert end-to-end.
Also done: installed ruff (was in the dev extra but absent) and made src/ and
  the metric tests clean; added tests/**/*.py per-file-ignores for S101/T201,
  without which pytest asserts fail lint -- upstream's own test_prune.py does.
  ruff format left alone: nothing in the repo is format-clean and ruff warns
  COM812 conflicts with the formatter.
One thing the human needs to decide:
  - Redundancy is a Phase 3 design question that is already answerable: the
    metric covering costs 45x the volume bound and ~90x the box strategy at
    poly_order=4, and poly_order=5 exceeds a 500k cap. A_4* buys 2.8x;
    anisotropic spacing in the metric eigenbasis buys most of the rest (two of
    four axes are thinner than one child's covering radius). Neither was done
    here -- the plan puts the lattice choice in Phase 3.
Corrected: the box strategy emits 4/32/512 children per parent at poly_order
  2/3/4, not the 27-81 an earlier session reported. ~1000 at poly_order=4 is
  intrinsic (it is the volume bound), so the metric approach cannot be cheaper
  than ~2x the box strategy, and the box strategy is not leaving headroom.
Open questions: O6 (unchanged, non-blocking).

## 2026-09-14 — Phase 2, step 6b
Done:
  - Route (i) chosen by the human (D15): `branch_offsets` / `branch_extents` added to
    the shared `branch` signature -- 3 `branch_func`s, 6 proxy methods, 6 overload
    stubs across the three `dynamic/` modules, plus the call in
    `prune.py::pruning_iteration_batched`. Chebyshev and circular accept and ignore
    them, with a comment saying why.
  - `dyn_poly_taylor.branch_func` dispatches on `self.tiling_strategy == "metric"` to
    `taylor.poly_taylor_branch_metric_apply`, and raises on an empty table rather than
    branching every parent into nothing.
  - `prune.py::Pruning._metric_branch_tables`: the per-stage precompute hook. Builds
    the tables once per level in Python and memoises them on
    `(t_obs_prev, t_obs_cur)`; returns empty arrays for every other strategy.
  - `_setup_pruning` now refuses `"metric"` on the Chebyshev and circular bases (D16).
  - config.py: `metric_branch_max` (500k) and `metric_ducy` (0 = follow `ducy_max`),
    both required by the hook -- see D17 for why `branch_max` could not be reused.
  - tests/test_branch_metric.py: +9 tests (`TestStrategyDispatch`, `TestPruneWiring`,
    `TestMetricConfigDefaults`). Full suite **126 passed**, `src` and the metric tests
    ruff-clean (`SLF001` added to the tests per-file-ignores: the wiring hangs off
    private members of `Pruning` on purpose).
Decisions fixed: D15 (route (i)), D16 (refuse metric off the Taylor basis),
  D17 (the two new config knobs).
Measured:
  - Per-stage table cost is negligible: ~1 ms at `poly_order=3` (869 children) and
    ~18 ms at `poly_order=4` (46051), so the memoisation is insurance, not a
    necessity. Child counts depend only on the *ratio* `t_cur / t_prev`, not on the
    absolute interval -- 869 at both 0.262->0.524 s and 16.78->33.55 s -- as the
    `f0**2` proportionality of D11 predicts.
Status of "metric" end-to-end: **still inert, but for a different reason.** With
  step 6b in place a real `prune_dyp_tree` run at `poly_order=3` now branches,
  validates, resolves, shift_adds and scores metric children successfully, and dies at
  `transform` -- because **Phase 2 step 2 (`transforms.py`) was never implemented**.
  `shift_taylor_errors` and `shift_taylor_full` still raise on `"metric"`, which also
  breaks `generate_branching_pattern` and so the threshold scheme. Step 2 is not
  mechanical: D12 stores only the ellipsoid's bounding box in column 1, and a bounding
  box does not determine the ellipsoid, so the plan's "transform `g` exactly" cannot be
  done from leaf state alone. Options written up under "Still to do in Phase 2".
Open questions: O6 (unchanged), O7 (new: which ducy sets the harmonic weighting; the
  default is the least conservative choice and should be settled with O6 in Phase 3).
Next session starts at: Phase 2 step 2 (`transforms.py`), then steps 4 and 5.

## 2026-09-14 (b) — Phase 2, step 2
Done:
  - Fourth option chosen by the human (D18), with (c) as an unused fallback:
    `taylor.metric_transform_extents` computes the re-centred child extents once per
    level at `f0 = 1`; `transform_extents` joins the shared `transform` signature
    (3 `transform_func`s, 6 proxies, 6 overload stubs) and `_metric_branch_tables`
    became `_metric_stage_tables`, returning all three arrays keyed on
    `(t_obs_prev, t_obs_cur, delta_t)`.
  - `poly_taylor_transform_batch` gained the `"metric"` branch: values via
    `shift_taylor_params` as always, column 1 from the table, `d_0` left at zero.
  - `core/taylor.py::generate_bp_poly_taylor_metric` + dispatch in
    `config.generate_branching_pattern` (D19), so the threshold scheme can be built.
    `generate_branching_pattern_approx` refuses `"metric"`.
  - tests/test_branch_metric.py: +12 (`TestTransformExtents`, `TestTransformDispatch`,
    `TestMetricBranchingPattern`). Full suite **138 passed**, ruff clean.
Decisions fixed: D18 (per-stage transform extents), D19 (exact B(s), not simulated).
Verified: `tiling_strategy="metric"` completes a real `prune_dyp_tree` run over three
  pruning levels, branching pattern included — see "Verified end-to-end" above.
  `aggressive` completes the identical run unchanged.
Found, not fixed: `report_func` discards the result of its own transform call upstream.
Open questions: O6, O7 (both unchanged, both for Phase 3).
Next session starts at: Phase 2 steps 4 and 5, then Phase 3.

## 2026-09-14 (c) — Phase 2, window-convention correction (D20/D21)
Done:
  - Found while scoping step 4: `metric_branch_tables` (and `metric_transform_extents`,
    which copied it) passed `coord[1]` to `poly_phase_metric` as an interval endpoint.
    `coord[1]` is a **half-width** and the epoch is the window's **centre**, so the
    metric was averaged over half the right length with the epoch at the edge.
  - Corrected all three windows (parent symmetric, child offset by `delta_t`, post-
    transform symmetric about the new epoch); `metric_transform_extents` no longer
    needs `delta_t` or the shift matrix at all.
  - `m_max_from_eta` box span -> `max(t_end - t_ref, t_ref - t_start)` (D21).
  - Full suite **141 passed**, ruff clean; `"metric"` and `aggressive` both still
    complete the 3-level end-to-end run.
Decisions fixed: D20 (the correction), D21 (the bridge's box span).
Re-measured:
  - Phase 2 redundancy table essentially unchanged (overheads 4.9 / 13.6 / 45.7 vs the
    recorded 4.9 / 13.6 / 45.0) — the covering depends on the parent/child ratio, which
    the error largely cancels out of.
  - D8 is the casualty: the recorded "56x on one axis at poly_order=5" was an artefact
    of the one-sided window. On the symmetric window it is **6.4x**. The qualitative
    argument for this project stands; the quantitative case is ~10x weaker than
    recorded. Phase 3 must quote the corrected numbers.
Not done: step 4, which is what I was scoping when this turned up.
Open questions: O6, O7 (unchanged).
Next session starts at: Phase 2 step 4 (the resolve diagnostic), on a metric that now
  means what it says.

## 2026-09-14 (d) — Phase 2, step 4
Done:
  - `taylor.metric_resolve_mismatch`: replays resolve's forward map, looks up the cell
    centre actually loaded, and prices the residual in the base-segment metric
    (`poly_order=2`, since `G0` grids only accel and freq).
  - `Pruning._log_resolve_mismatch` behind `cfg.metric_resolve_diagnostic`, sampling 64
    parents per level. Deliberately outside the `@njit` loop: it re-branches the sample
    in Python rather than instrumenting the hot path.
  - 3 tests: on-grid children cost < 1e-18; a half-cell offset matches a hand
    computation to 1e-5; and the step 4 finding itself is pinned.
  - Full suite **144 passed**, ruff clean.
Result: **the gate tripped, hard** — p95 of 454 / 598 / 5770 against `m_max = 0.2`.
  See D22. The diagnostic is sound: the same code on `aggressive` children returns
  p50 = 7.0e-10.
Root cause: the metric branch sizes the parent region from the baseline alone
  (`m_max` ellipsoid of the stage metric) and never consults the region the parent
  actually occupies. On short baselines that ellipsoid exceeds the whole search space
  by 9.2e7x in jerk and 4.1e6x in accel, so children land far outside `param_limits`
  and `resolve` clamps them. The metric is fine; the *scope* is wrong.
Consequence: the Phase 2 child-count and redundancy findings measure the cost of
  covering a region the parent never had, and should fall once the scope is fixed.
Open questions: O6, O7, and now the step 4 design choice (four options under
  "Still to do in Phase 2"), which needs a human decision before Phase 3.
Next session starts at: that decision, then step 5.

## 2026-09-14 (e) — Phase 2, the D22 fix
Done:
  - Option (3) chosen by the human, the principled one: regions are carried explicitly
    as forms `{d : d^T A d <= 1}`, per stage rather than per leaf (D23).
  - `core/metric.py`: `region_form_from_box`, `region_axis_extents`,
    `region_from_metric`, `region_fits_in_one_child`, `lattice_children_for_region`.
    `lattice_children` is now a thin wrapper, so the old two-metric call is the special
    case it always looked like.
  - `core/taylor.py`: `metric_branch_tables` takes the region and returns the new one;
    `metric_transform_region` replaces `metric_transform_extents`;
    `metric_seed_region` reads a seed's FFA cell out of column 1.
  - `prune.py`: `_region_form` seeded in `initialize`, advanced per level. The stage
    cache is gone — the tables are history-dependent now, and were only ever computed
    once per level anyway.
  - `generate_bp_poly_taylor_metric` replays the same recursion, so `B(s)` reports the
    single child the guard really emits.
  - D24: column 1 is a full span, not a half-width; C3 corrected, and the metric path
    now writes `2 * region_axis_extents`.
  - Full suite **146 passed**, ruff clean.
Result: the step 4 gate passes at the control's own floor — p95 7.0e-10 / 7.0e-10 /
  2.8e-9 vs `m_max = 0.2`, matching `aggressive` exactly. Children stay inside
  `param_limits`. Both strategies complete the 3-level end-to-end run.
Consequence worth flagging: on the smoke config the metric strategy now emits ONE child
  per level, because `m_max = 0.2` over a 1-second observation is looser than an FFA
  cell. Correct, not inert: on 64 segments refinement starts at level 9 and runs
  9/33/33/29/27/27. The old 869/171/93 counts and the 45x redundancy are retired.
Open questions: O6, O7. The step 4 design question is closed by D23.
Next session starts at: step 5 (the consumer audit, now with D24 to check), then
  re-measuring cost on a schedule long enough to refine, then Phase 3.

## 2026-09-14 (f) — branching-factor comparison
Done:
  - Pushed the branch (8 commits) to origin/metric-gridding.
  - Figure 7 analogue measured, config-only: metric `prod B(s) = 1.1e26` against the
    box strategies' **27** on a 32-segment, 16.8 s, poly_order=3 schedule.
  - D25 found and fixed on the way: the child region must be parent ∩ ellipsoid, not
    the ellipsoid. Worth 4.3e33 -> 1.1e26. `metric.region_intersection` + 4 tests.
  - Full suite **153 passed**, ruff clean.
Diagnosis (all measured, see the section above):
  - The box leaf already carries only 1-12% of `m_max`, so the criterion is not the
    problem — the metric demands *less* refinement than the box performs.
  - The covering pays the lattice thickness at every level for a volume gain of only
    1.5-1.9x, and that compounds over ~24 levels.
  - Simulated deferral at R=3 gives `prod B(s) = 21` vs the box's 27 — competitive.
    The cliff is sharp: R=2 is 2000x worse, R=4 never refines.
Verdict: viable, but only with deferral, and the R=3 operating point buys cost parity
  with a *single* late branch event. Whether that detects as well as the box's early
  branching is an injection-recovery question, and it is now the central Phase 3
  experiment rather than a side one.
Open questions: O6, O7, and the deferral factor R (new).
Next session starts at: implement `metric_defer_factor`, re-measure on the Phase 3
  schedule, then step 5.

## 2026-09-14 (g) — metric_defer_factor, and a corrected verdict
Done:
  - `metric_defer_factor` implemented end to end (D26): `metric.region_overhang`,
    the generalised `region_fits_in_one_child`, threaded through `metric_branch_tables`,
    `generate_bp_poly_taylor_metric`, `config.py` and `prune.py`. Default **1.0**, so
    the shipped behaviour keeps the `m_max` guarantee. 10 tests.
  - Full suite **163 passed**, ruff clean.
CORRECTION to the previous session entry: "the metric strategy is not competitive, by
  24 orders of magnitude" was measured on ONE small config and does not generalise. At
  equal worst-case mismatch the sign flips with the regime — 103x worse at 16.8 s /
  po=3, but 3.5x / 6.3e6 x / 1.8e11 x BETTER at 67 s po=3, 67 s po=4 and 268 s po=4.
  See "Regime dependence" above for the table.
Why the small config misleads: there `aggressive` under-reports its region by 1.003x
  and its leaf carries m=0.024 at prod B(s)=27 — it is both cheap and honest, so there
  is no gap to fix and the covering overhead is pure loss. The premise of the project
  only bites once the baseline and poly_order grow, and there it bites hard.
Open questions: O6, O7, and (m_max, R) jointly — all four are scale factors on the same
  criterion and Phase 3 should settle them together.
Next session starts at: step 5, then Phase 3 with the target-regime config.

## 2026-09-14 (h) — Phase 2, step 5
Done:
  - Audited all ten Phase 0 consumers of column 1 against what the metric path writes
    (D27). Resolution: "use axis extents" everywhere; nothing needs to skip, and no
    code change was required — D12/D18/D24 had already put the right thing in column 1.
  - Verified end to end rather than by reading: a real `prune_dyp_tree` run under each
    strategy publishes comparable spans for the same candidate (jerk 15.99 vs 16.15,
    accel 8.39 vs 8.76, freq 0.0741 vs 0.2223 Hz).
  - 4 tests pinning the publishing path, including that the published number is a full
    span (D24) at the point a user reads it, and that `validate` really is a no-op.
  - Full suite **167 passed**, ruff clean.
Flagged, not fixed (pre-existing, applies to the box strategies too):
  - `poly_taylor_report_batch` combines column 1 in quadrature, treating a full span as
    a 1-sigma uncertainty and assuming axis independence. Neither holds.
  - The published span is a bounding box, not a covering box; it does not tile.
Phase 2 is now complete: steps 1-6 all done.
Open questions: O6, O7, and (m_max, R) jointly.
Next session starts at: Phase 3, on a target-regime config (67 s+, poly_order 4), with
  the four scale factors settled together.

## 2026-09-14 (i) — Phase 3 setup, and a second correction to the cost comparison
Done:
  - Phase 3 scoped (D29): the plan's circular-orbit config cannot run, since `"metric"`
    is refused on the circular basis (D16). Human chose the **Taylor analogue** at
    `poly_order=4`, 67 s, 64 segments. Figures 8/11/12 will be reproduced in spirit,
    not literally; the circular extension stays in Phase 4 where the plan put it.
  - D28: cherry-picked PR #7 (open upstream, unmerged) to unblock threshold
    recalibration. Verified `DynamicThresholdScheme.run()` now completes.
  - `docs/metric_gridding/phase3_config.py` — one config for every Phase 3 experiment.
  - `cost_at_equal_sensitivity.py` rewritten around a validated box reference.
CORRECTION, the second on this comparison: the previous entry compared the metric
  against a baseline computed with `conservative` propagation while calling it "the
  box". That is a different, dearer strategy, and it flattered the metric by up to 6
  orders. Validated the loop this time — it reproduces `generate_branching_pattern`
  exactly for `aggressive` — and rebuilt the table against both baselines.
  Corrected reading: the metric sits ~1/3 of the way from `aggressive` to
  `conservative` on a log scale, i.e. far cheaper than the only gap-free box option and
  far dearer than the gappy one. Whether that is a good trade is a detection question.
Open questions: O6, O7, (m_max, R) jointly.
Next session starts at: Phase 3 step 2 — recalibrate thresholds for both strategies at
  P_d = 0.1 on the Taylor-analogue config, then step 3 (injection-recovery), which is
  now the decisive experiment.

## 2026-09-14 (j) — Phase 3 reduced-grid validation
Done:
  - `injection_recovery.py`: inject, prune, judge recovery by mismatch in the
    full-baseline metric. Runs for all three strategies.
  - Validated on a reduced grid before committing to the full one, which was the right
    call — it moved the Phase 3 config twice and caught a false negative.
Findings:
  - D30: the plan's 18-min config needs a 3.8 TB FFA fold array and cannot run here.
    And 67 s -- which I had picked -- has a degenerate base grid `[1,1,1,12]`, where all
    three strategies return an identical mismatch and nothing is discriminated. Phase 3
    now runs at 268 s / 64 seg / po=4: 0.1 GB, FFA 2.9 s, grid `[1,1,5,467]`.
  - D31: the metric appeared to fail recovery (best m 2.8e6, score 4.9). It does not --
    with zero thresholds it recovers at m=0.0037 / score 19.0, better than aggressive on
    the same setting. The failure was a shared, uncalibrated threshold ramp. The plan's
    "recalibrate before comparing sensitivity" is now demonstrated, not assumed.
Open questions: O6, O7, (m_max, R) jointly.
Next session starts at: Phase 3 step 2 proper -- `DynamicThresholdScheme` per strategy
  at P_d = 0.1 on the 268 s config -- then step 3 on a real injection grid.

## 2026-09-14 (k) — fixing D32 upstream, and finishing step 2
Done:
  - Worked the fix on its own branch in its own worktree, off PR #7 (whose SIGABRT masks
    both bugs), with an isolated venv: `worktrees/fix-viterbi`,
    branch `fix-viterbi-zero-states`.
  - **Upstream issue #8** — `run_stage_legacy`/`run_stage_improved` wrote zero-filled
    records. Root cause is a numba miscompilation: given a jitted callee returning a
    one-element structured array, `dst[i] = make()[0]` writes zeros, `r = a[0];
    dst[i] = r` writes zeros, field-by-field loses bool fields, and only
    `dst[i:i+1] = make()` round-trips. Minimal repro in the issue.
  - **Upstream issue #9** — `trials_scheme` returns `-inf` when the pattern starts
    unbranched, emptying the beam. Found only after #8 was fixed (D34).
  - `tests/test_thresholding.py`: 11 tests. Verified they fail on the parent commit and
    pass after, for each fix separately, by reverting just the src file.
  - Both cherry-picked into `metric-gridding`; suite **178 passed**.
  - Phase 3 step 2 complete: all three strategies now reach P_d = 0.1031. See the table
    above.
Result worth stating: at equal detection probability the metric costs 1070x
  `aggressive` and 3.7e9 x less than `conservative` — the same "quarter to a third of
  the way across on a log scale" position the branching comparison found, now measured
  on the axis the plan actually cares about.
Open questions: O6, O7, (m_max, R) jointly.
Next session starts at: Phase 3 step 3 — injection-recovery with these schemes, which
  is the decisive experiment.

## 2026-09-14 (l) — Phase 3 step 3
Done:
  - Upstream PR #10 opened (fixes #8 and #9, based on #7).
  - Injection-recovery grid: 3 strategies x 6 S/N x 3 realisations, with the
    recalibrated schemes. Results and verdict in `03_results.md`.
  - **D35 — `max_sugg` is not a neutral knob for the metric strategy.** Its branching
    is spiky (17 549 children from one parent), so a shared `max_sugg = 2**14` cannot
    hold even one parent's offspring and the true track is trimmed mid-burst. At S/N 15
    the metric scored 0/3; with `2**17` and nothing else changed it scores **3/3**, at a
    median mismatch of 2.5e-4 against `aggressive`'s 1.5e-3 and an equal score.
    The second time a shared resource setting has made the metric look broken -- D31 was
    the first. Any future comparison must give each strategy the resources it needs and
    count that as part of its cost, not silently starve one of them.
Verdict (Taylor analogue): the metric works, matches `aggressive`'s detection and
  localises better, and costs ~3 orders of magnitude more -- 8x buffer, ~450x
  wall-clock, ~1070x complexity at equal P_d. Three independent measures agreeing to a
  factor of a few.
  But the config cannot answer the central question: `aggressive` under-reports its
  region by only 1.02-1.33x here, so there is no coverage gap to fix and only the
  metric's cost is visible. NOT a reason to abandon -- a reason to say the decisive
  experiment has not been run. See `03_results.md` for what would change that.
Open questions: O6, O7, (m_max, R) jointly, and now: which affordable regime, if any,
  makes `aggressive` visibly fail.
Next session starts at: scanning for a Taylor config where `aggressive` under-reports
  badly; failing that, Phase 4's circular extension or a bigger machine.

## 2026-09-14 (m) — the benefit side, measured at last
Done:
  - Costed `quadrature` at the Phase 3 config, which had never been done: **1.10e17**
    against `aggressive`'s 1.51e12 and the metric's 1.35e19.
  - **D36 — measured the coverage gap, and it is zero.** `T` is lower-triangular with
    unit diagonal (C4), so the sheared lattice of leaf centres keeps each axis's period
    and boxes of unchanged width tile exactly, for any shear. `aggressive` records
    exactly `spacing * |diag(T)| == spacing`, which is precisely sufficient.
    Uncovered fraction is 0.0000 for all three box strategies at `delta_t` from 0.5 to
    100, verified by forward substitution and independently by brute force over a
    +/-400 lattice in 2D.
    Took three attempts; the first two searched a fixed neighbour block instead of
    solving for the covering index and reported spurious holes, once claiming
    `aggressive` leaves 100% uncovered. `docs/metric_gridding/coverage_loss.py` carries
    the working version and its self-checks.
Consequence: the project's premise does not hold as stated. `quadrature` and
  `conservative` buy no coverage over `aggressive` at 1e5 x and 1e20 x the cost, and the
  metric pays ~1070x at equal P_d for a guarantee `aggressive` already provides.
  Consistent with step 3, where `aggressive` matched or beat the metric.
Caveats that could overturn it, in order: pruning is not modelled (a signal migrating
  into a thresholded-away leaf is still lost, and a better tiling does not fix that);
  this is the Taylor basis only, and the Chebyshev and circular transforms are not
  unit-diagonal triangular; and it contradicts the paper as this plan reads it.
Open questions: O6, O7, (m_max, R) -- all now moot unless D36 is overturned.
Next session starts at: checking D36 independently, and looking for the gap in the
  Chebyshev/circular transforms where the argument does not apply.

## 2026-09-14 (n) — read the paper at last
Done:
  - Read `paper/pruning1.tex` §3.4 and §5.2.2-5.2.4 (numbering verified: 5.2.1 Reference
    Frame Dilemma, 5.2.2 Moving Grid and Axis Misalignment, 5.2.3 Orthogonal Basis,
    5.2.4 Branching Strategy) and Appendix `app:optimal_gridding`.
  - D37: retracted the claim that D36 contradicts §5.2.4. It sharpens §5.2.2's
    "partially mitigated by natural overlap" to "exactly mitigated", and §5.2.4 already
    predicts that exact shear tracking gives no computational advantage.
  - D38: §3.4 shows the shipped box grid is already metric-derived -- the `2**(k-1)`
    coarsening IS the paper's analytic diagonalisation of the parameter metric. D5 had
    us drop it, so our covering competes against a metric-corrected box while paying
    full price. The residual prize is only the off-diagonal structure.
Process failure worth naming: I asserted a contradiction with a section I had not
  opened, from the plan's paraphrase, while the plan's own ground rules say to verify
  citations. The paper had been in `paper/` the whole time. Fourth unforced error of the
  session, same root cause as the others -- trusting a baseline instead of checking it.
Open questions: whether anything in the plan's premise survives D37 and D38. The one
  piece that does is the paper's own request: quantify the sensitivity loss of
  aggressive tiling. That is a different (and much cheaper) measurement than building a
  metric covering, and it is not what this branch has built.
Next session starts at: deciding whether to pivot to the sensitivity-loss measurement,
  or to stop.

## 2026-09-14 (o) — the sensitivity loss, measured; and the verdict
(Written 2026-09-16: the session that produced commit 04183d4 was cut off by a machine
reset before it wrote this entry. The code and `03_results.md` in that commit are the
record; this reconstructs the log from them.)
Done:
  - Pivoted as (n) proposed, to §5.2.4's actual request rather than more metric work.
  - `docs/metric_gridding/sensitivity_loss.py`: samples offsets uniformly inside each
    leaf and reports the incurred phase error against the promised `eta / N_b`.
    268 s / 64 seg / `po=4`, in units of `eta / N_b`:
        strategy      median stage   worst corner   median mismatch
        aggressive    3.35           16.80          0.0240
        quadrature    3.58           14.38          0.0258
        conservative  3.42           14.20          0.0251
  - **D39 — the three tilings are indistinguishable in sensitivity**, across twenty
    orders of magnitude of cost (`prod B(s)` 1.51e12 / 1.10e17 / 1.48e32).
    Mechanism: `branch_param_padded` sets `num_points = ceil(dparam_cur / dparam_new)`,
    so the child cell is `dparam_cur / num_points`, landing just under the criterion
    step `dparam_new` whatever was transported. Confirmed by the final cell
    half-widths: all within ~2x of the criterion `[9.7e-3, 0.163, 3.64, 122]` and of
    each other, the scatter being the integer `ceil`. **The transported width never
    sets the cell; the criterion does.** It sets only how much redundant subdivision
    happens on the way, so it can buy cost and never sensitivity.
  - D39 generalises to the metric strategy: a covering is just another transport rule,
    so it cannot improve sensitivity either, and it measurably worsens cost (~1070x at
    equal `P_d`, D35). Together with D36 (zero coverage gap) and D38 (the shipped box
    grid is already the paper's metric diagonalisation), this closes the question the
    branch was opened on.
  - D40 — two self-checks pin what `eta` means: the naive per-axis grid costs exactly
    one tolerance at a full step (promise kept per axis), while the shipped grid is
    coarser by exactly `2**(k-1)` on the order-`k` axis. So `eta` does *not* mean
    "phase error below `eta / N_b`".
Conventions fixed: none new.
Verdict (Phase 3 deliverable): **abandon** the metric covering. Not because the metric
  is wrong, but because the sensitivity it was meant to recover is not lost in the first
  place (D36), the box grid is already metric-derived (D38), and the cell size is set by
  the criterion rather than by transport (D39).
Worth reporting upstream, independent of all metric work: under the recommended
  `aggressive` Taylor default a typical signal costs ~3.4x the nominal phase budget and
  ~2.5% in amplitude, a corner signal ~15x. Choosing `eta` from the grid criterion alone
  over-estimates sensitivity by roughly a factor of three, most of it the `2**(k-1)`
  coarsening compounded across axes and the rest multi-axis corner addition.
Open questions: O6, O7 and `(m_max, R)` are moot under D39. Two live ones remain, both
  places D36/D39 provably do not reach:
  - the pruning interaction, unmodelled here (deterministic geometry only): a signal
    near a cell edge may sit in a leaf that was thresholded away, and a strategy
    claiming more territory keeps more such signals alive. Measurable only by injection.
  - the Chebyshev and circular transforms, which are not unit-diagonal triangular, so
    neither D36 nor D39 carries over. If the paper's gap is real, it is there.
Next session starts at: writing up the `~3.4x` / `~15x` sensitivity-loss result for
  upstream (it stands on its own), then deciding whether to test the pruning
  interaction by injection or to close the branch.

## 2026-09-16 (p) — upstream report drafted; numbers re-verified
Done:
  - Re-ran `sensitivity_loss.py` from a clean session: `self_check` passes and all six
    headline numbers reproduce exactly (3.35/16.80/0.0240, 3.58/14.38/0.0258,
    3.42/14.20/0.0251). Also re-derived the final half-widths and `prod B(s)`
    (1.51e12 / 1.10e17 / 1.48e32) independently. D39 stands.
  - D41 — corrected the D39 mechanism as stated. The claim "the cell lands just under
    the criterion" is only the `branch_param_padded` half. The shipped code also has a
    guard (`core/taylor.py:148-154`): when an axis's accumulated shift is below `eta`
    it is not branched at all and keeps `dparam_cur`, so cells can sit *above* the
    criterion, up to ~2x. Verified the reimplementation in `sensitivity_loss.py`
    matches the shipped guard. The conclusion is unchanged and slightly stronger: the
    criterion fixes the scale on both sides, and transport survives only as the `ceil`
    remainder plus the guard slack.
  - D42 — verified against the paper that the `2**(k-1)` coarsening in
    `psr_utils.poly_taylor_step_f:93-94` is faithful. §3.4 eq. `dk_optimal` carries
    `2**(2k-1)` in terms of `Tobs` while appendix eq. at line 1585 carries `2**(k-1)`
    in terms of the half-span `t_s`; these agree since `Tobs = 2 t_s`, and the code
    uses the half-span form. No discrepancy, contra a first reading of §3.4.
  - `docs/metric_gridding/04_upstream_report.md`: draft report for upstream, leading
    with the strategy-independence result rather than the metric work. Not posted.
Conventions fixed: none new.
Note: upstream has GitHub Discussions disabled, so Phase 4's "discussion issue" would
  be a plain issue on `pravirkr/pyloki`.
Open questions: unchanged — the pruning interaction (needs injections) and the
  Chebyshev/circular transforms.
Next session starts at: posting `04_upstream_report.md` once reviewed, then either the
  injection test of the pruning interaction or closing the branch.

## 2026-09-16 (q) — the pruning interaction, and D39's conclusion refuted
Done:
  - `docs/metric_gridding/pruning_multiplicity.py`: follows the true signal down the
    tree with the shipped `poly_taylor_branch_batch` / `poly_taylor_transform_batch`,
    keeping every leaf within one claimed width of the truth and transporting the truth
    with the same `T`. Measures the excursion to the **nearest** leaf centre, not to the
    centre of the cell the signal happens to sit in.
  - **D43 — D39's conclusion is wrong, and the mechanism it rests on is still right.**
    `sensitivity_loss.py` sampled offsets inside a cell and measured from *that* cell's
    centre, which is only the right quantity at covering multiplicity 1. The box
    strategies over-claim under transport, so sibling cells overlap and the signal is
    covered many times over. 12 random signal positions, `n_seed=1` (81 base cells):

    | strategy | nearest excursion (median) | worst signal | multiplicity | mismatch at nearest |
    |---|---|---|---|---|
    | `aggressive` | 1.89 | 3.09 | 1 | 0.00878 |
    | `quadrature` | 0.617 | 1.07 | 3696 | 0.00210 |
    | `conservative` | 0.615 | 1.35 | 13042 | 0.00205 |

    So redundancy *does* buy sensitivity: 3.1x in phase error (per-signal 1.8-4.9x), and
    `aggressive` alone fails to deliver `eta / N_b` (1.89 > 1) while the redundant
    strategies keep it. D39's cell-size mechanism (D41) is unaffected and still verified
    -- the cell is criterion-set -- but "it can buy cost and never sensitivity" does not
    follow from it, because effective template *density* is not the cell size.
  - **D44 — the gain is 0.34% in S/N, which settles the pruning question without
    injections.** Converting the nearest-leaf mismatch to amplitude: 0.44% loss for
    `aggressive` against 0.105% for `quadrature`, i.e. an S/N ratio of 1.0034. Against
    that, step 2's recalibrated ladders differ by far more (top threshold 7.70 for
    `aggressive` vs 9.10 for `conservative`), and at equal `P_d` `aggressive` is 2^41.8
    cheaper than `conservative`. A 0.34% amplitude edge cannot survive a threshold
    ladder that much higher, and no injection run at these statistics could resolve it.
    The pruning interaction is therefore real in mechanism and negligible in magnitude.
  - **D45 — `conservative` is strictly dominated by `quadrature`**: identical nearest
    excursion (0.615 vs 0.617) and mismatch (0.00205 vs 0.00210) for `2^15` times the
    cost (1.48e32 vs 1.10e17). Worth reporting upstream on its own.
  - **D46 — the corner of the shipped optimal grid costs exactly
    `(2**(poly_order-1) - 1/2) * eta/N_b`.** Axis `k` of eq. `dk_optimal` contributes
    `(Delta d_k/2) * (f_max/c) * t_s^k/k! = 2**(k-2) * eta/N_b`, and the sum over
    `k = 1..k_max` telescopes. Verified against `poly_taylor_step_d_vec` for orders 2-8,
    exact to machine precision and independent of `t_s` and `f_max` (both cancel):
    1.5 / 3.5 / 7.5 / 15.5 / 31.5 / 63.5 / 127.5. So 7.5x at the default `poly_order=4`
    and 15.5x at 5, which is what a circular-orbit search runs. This is the largest term
    in the whole sensitivity question and it is in the criterion, not in any tiling.
    It also **reconciles the 7.5x single-cell figure with the 16.80x** reported in (o):
    the two measure different things, and the difference is the moving grid's displaced
    validity window, not the grid. Evaluating the same nominal corner over a window
    displaced by `delta_t` reproduces 7.500 at `delta_t=0` and grows from there; the
    realized-vs-nominal cell width (D41) does *not* explain it, since nominal cells give
    a larger worst corner (19.83) than realized ones (13.94).

Process note: two bookkeeping artefacts were caught and fixed before they became
  findings. Seeding a single base cell charged the signal's drift out of that family to
  the strategy and produced a spurious 86.5x worst case for `aggressive` (3.09 once a
  81-cell block is seeded, as the real search does). A +/-3-cell tracking window was
  also unaffordable and unnecessary; results are stable between window 1 and 2, and
  against an 10x larger tracked-leaf cap.
Conventions fixed: nearest-leaf excursion, not own-cell excursion, is the sensitivity
  quantity whenever multiplicity > 1. `sensitivity_loss.py` numbers stand as own-cell
  measurements and must not be quoted as sensitivity for the redundant strategies.
Open questions: the Chebyshev/circular transforms remain untested (not unit-diagonal
  triangular, so neither D36 nor D41 carries over). The pruning interaction is closed.
Next session starts at: `04_upstream_report.md` has been rewritten around D46 (headline)
  plus D43-D45; it is reviewed-but-unposted, pending the human. The old headline
  ("buys no sensitivity") was refuted by D43 and is gone.

## 2026-09-16 (r) — RETRACTIONS: D43/D44/D45 withdrawn, D46 reinterpreted
Done:
  - **D47 — `pruning_multiplicity.py` is unsound for multiplicity > 1, and the proof is
    internal.** Widening the tracking window can only *add* candidate leaves, so the
    reported minimum must be non-increasing. Measured, 6 signal positions:

    | window | `aggressive` | `quadrature` | `conservative` |
    |---|---|---|---|
    | 1.0 | 1.6889 | 0.6173 | 0.6681 |
    | 2.0 | 1.1751 | 0.5662 | 0.7931 |
    | 3.0 | 1.0897 | 0.9957 | 4.6284 |
    | 4.0 | 1.0857 | 1.0554 | 10.4606 |
    | 6.0 | 1.0857 | 1.0314 | 6.4246 |

    `aggressive` is monotone and converges; the other two rise and wander, which a true
    minimum cannot do. Cause: `MAX_TRACKED` keeps the best 20000 leaves *by current
    excursion*, discarding leaves that would become nearest after later branching. It
    binds in **37-60 of 63 stages for `quadrature`, 58-61 for `conservative`, and 0 of
    63 for `aggressive` at window <= 3**. The earlier cap check (20000 vs 200000 at
    window 1) was not a sensitive enough test and gave false reassurance.
  - **D43 quantitative claim — WITHDRAWN.** The *principle* stands and is now sharpened:
    nearest-template rather than own-cell is the sensitivity quantity whenever the
    signal has more than one template near it -- which, per D48, is true even on the
    base grid, so it was never specific to redundant tiling. The 1.89 / 0.617 / 0.615
    table and the "3.1x gain" are gone: 1.89 was a too-narrow window (converges to
    1.09) and the other two are cap artefacts.
  - **D44 — WITHDRAWN.** It was a conversion of the withdrawn 3.1x. **This reopens the
    pruning interaction**, which (q) wrongly recorded as closed; the injection campaign
    is not cancelled.
  - **D45 — WITHDRAWN as to sensitivity** (rested on 0.615 vs 0.617, both cap
    artefacts). The cost side, `prod B(s)` 1.48e32 vs 1.10e17, is unaffected.
  - **D48 — D46's arithmetic stands; its interpretation is retracted.** The corner
    figure is the distance to the *containing cell's* centre, and a search enumerates
    every grid point, so the operative quantity is the covering radius. Verified at
    `k_max=4`: own-cell 7.5, best adjacent centre (±half-step per axis) **0.5829**,
    best over integer offsets in -2..2 **0.5000**. The coarsened coefficients are
    Chebyshev-like, so a sign choice makes the polynomial nearly cancel over the
    interval instead of adding. Also verified: the naive uncoarsened grid's own-cell
    corner is exactly `k_max/2`, so the geometric growth is *entirely* the `2^(j-1)`
    coarsening, and the mechanism is not `|T_k| <= 1`. So "`eta` is optimistic by 7.5x"
    is wrong and must not be published; on this evidence the grid covers at O(1) and the
    coarsening is economization working, not a defect.
  - **D49 — what survives as a defensible number.** `aggressive`'s nearest-template
    error converges to **~1.09 tolerances** (window-converged at 4 and 6, cap never
    binds, 6 signal positions). Under the shipped default the nearest template sits
    roughly one tolerance from the signal. Cap-free, and the one Table 2 entry to keep.
  - **D39's conclusion is back to undetermined**, *not* vindicated. There is no sound
    measurement either way for the redundant strategies.
Process failure worth naming: I published a table into a report and told a peer session
  the pruning question was "closed analytically" on the strength of a measurement whose
  convergence I had not tested. The window sweep that broke it costs seconds and should
  have been the first thing run, not the last. The peer's internal-inconsistency
  argument -- that my own 1.89 contradicted my own 7.5 headline three sections apart --
  is what exposed it; I had both numbers in one document and did not reconcile them.
Conventions fixed: any "nearest template" measurement must report a window/cap
  convergence sweep, and must state at how many stages the cap binds. A minimum that is
  not monotone in the search width is an artefact, full stop.
Open questions: a cap-free method for multiplicity > 1. The leaf set is a Cartesian
  product per axis, so branch-and-bound over the product with a per-axis bound on the
  sup-norm contribution is the likely route; materialising leaves and truncating cannot
  work. Until then the strategy comparison is open, and so is the pruning interaction.
Next session starts at: either building the branch-and-bound nearest-template search, or
  accepting D48/D49 and rewriting the report around only those.

## 2026-09-17 (s) — the sound method, and D43's direction restored
Done:
  - **D50 — `docs/metric_gridding/nearest_template.py`, exact and validated.** The
    structural fact that makes it possible: `branch_param_padded` derives a child's
    offset from `(dparam_cur, dparam_new)` alone, and every leaf at a stage shares
    those, so **child offsets are identical for every parent**. Hence the leaf set at
    stage `S` is *exactly* a Minkowski sum
    `M_0.Seed (+) M_1.O_1 (+) ... (+) O_S` of small per-stage offset sets -- 1e34 leaves
    in a representation of a few hundred points per stage, with no approximation.
    Minimising the sup-norm over that sum is done by branch and bound with the
    admissible bound `|<p,b(t)>| - sum_{i>=m} max_a |<a,b(t)>|`, so nothing is discarded
    except against a bound. Where the budget runs out it returns the incumbent, which is
    a real leaf and therefore a **valid upper bound**, flagged as not exact. It never
    passes a truncated search off as a minimum, which was D47's defect.
    `tests/test_nearest_template.py` pins parent-independence of the offsets, exact
    agreement with brute-force enumeration (9 cases, all three strategies, rel 1e-12),
    monotonicity in the seed width, and the seeded-incumbent verdicts.
  - **D51 — `aggressive`'s nearest-template error, exactly: median 1.069**, range
    [0.240, 2.115], EXACT in all 90 cells (6 signals x 15 stages). This *confirms* D49's
    ~1.09 from the old tracker, which is expected: the cap never binds for `aggressive`,
    so that one column was always sound. Under the shipped default the nearest template
    sits just above one tolerance from the signal.
  - **D52 — D43's direction is restored on a sound footing; its magnitude is not.**
    Seeding the incumbent with `aggressive`'s exact value turns the comparison into a
    decision problem, which is far cheaper than optimising. Over the same 90 cells:

    | vs `aggressive` | strictly closer (proved) | not closer (proved) | unresolved |
    |---|---|---|---|
    | `quadrature` | **89** | 1 | 0 |
    | `conservative` | 28 | 0 | 62 |

    Median proven gain where closer: **2.30x** for `quadrature`, 2.05x for
    `conservative` -- and because these come from upper bounds, each is a *lower* bound
    on the true gain. So redundancy does put a closer template near the signal, and
    `tiling_strategy` does affect sensitivity in the Taylor basis. The withdrawn 3.1x is
    replaced by ">= 2.30x"; `conservative` is expensive to resolve and there is no
    evidence against it, only absence of proof in 62 cells.
Conventions fixed: report `(value, exact)` from any nearest-template search and never a
  bare number; a budget-limited result is an upper bound and may only be used one-sided.
Still withdrawn: **D44 (0.34% in S/N) stays withdrawn** and has NOT been recomputed. The
  amplitude consequence of D52's >= 2.30x is unknown, so nothing here says whether the
  difference matters in practice, and the pruning interaction remains open -- D52 is
  geometry, and still does not model whether the closer leaf survives thresholding.
Open questions: the amplitude conversion; whether `conservative`'s 62 unresolved cells
  can be closed with a better bound; the Chebyshev/circular transforms.
Next session starts at: converting D52 into an amplitude number, carefully, or the
  injection test that D44's withdrawal reopened.

## 2026-09-17 (t) — the amplitude conversion, done without the metric
Done:
  - **D53 — `docs/metric_gridding/amplitude_loss.py`, and why it avoids the metric.**
    A phase residual does not attenuate the pulse, it *smears* it: the folded profile is
    `<p(phi - dPhi(t))>_t`, which in Fourier is exactly `P(k) * S(k)` with
    `S(k) = <exp(-2 pi i k dPhi(t))>_t`, the characteristic function of the residual
    distribution. Exact, with no second-order truncation and no harmonic-weighting
    choice -- the two things that made D44's metric route unsafe (O4 was still open, and
    the metric's mean-square mismatch was being applied to a template selected by a
    sup-norm criterion). Scored with `detection.scoring`'s own boxcar bank, cross-checked
    against an ideal matched filter.
  - **D54 — phase averaging is not cosmetic.** Scored at one pulse position the answer
    carried a discretisation artefact the size of the effect: a Gaussian smeared toward a
    flat top scores *better* against a boxcar bank than a sharp one, which produced
    **negative losses**. Averaging over 16 sub-bin positions (the signal's absolute phase
    is arbitrary and uniform) fixes it, and the boxcar and matched-filter answers then
    agree to ~30%. Pinned by `test_loss_is_monotone_in_phase_error`.
  - **D55 — the numbers.** 6 signals x 8 stages, actual offset vectors from the exact
    search, best of the 30 nearest templates, boxcar (matched in brackets):

    | pulse duty | `aggressive` loss | `quadrature` loss | `quadrature`'s S/N advantage |
    |---|---|---|---|
    | 0.05 | 5.34% (8.00%) | 2.26% (4.03%) | +1.56% (+1.38%) |
    | 0.10 | 1.73% (2.21%) | 0.75% (1.05%) | +0.47% (+0.40%) |
    | 0.20 | 0.40% (0.57%) | 0.17% (0.27%) | +0.11% (+0.11%) |

    Duty cycle is the dominant uncertainty -- an order of magnitude across the range --
    and `ducy=0.05` is marginal at `N_b=64` (sigma 0.75 bins), so the 0.10 row is the one
    to quote.
  - **D56 — sup-norm badly overestimates the loss.** At `ducy=0.10` the real offsets cost
    1.73% at a sup-norm of 1.04, where a *linear ramp* of the same sup-norm costs 4.78%.
    The residual attains its sup only briefly at the window edge, so what smears the
    profile is the residual's distribution, not its peak. Anyone converting a grid
    tolerance to an S/N loss via the peak phase error will be pessimistic by ~3x.
  - **D57 — the practical conclusion of the withdrawn D44 is restored, on a sound basis
    and with a different number.** The tiling choice is worth **0.1% to 1.6% in S/N**
    (0.47% at `ducy=0.10`), against `quadrature` costing `2^8.8` ~ 1e5 times more
    branching and, at equal `P_d`, a recalibrated top threshold of 9.10 against
    `aggressive`'s 7.70. A half-percent of amplitude cannot pay for either. So
    `aggressive` remains the right default -- not because the strategies are equivalent
    (D52 proves they are not) but because the difference is small and the price is not.
    Note the coincidence: withdrawn D44 said 0.34% and the sound answer at `ducy=0.10` is
    0.47%. It was roughly right by luck, from wrong inputs and an unsound conversion.
Caveat that must travel with D55: both losses are **upper bounds**. The template comes
  from a sup-norm search, not a loss search, so the best-scoring leaf may not be among
  the 30 supplied. The *difference* is therefore indicative rather than proven; the
  proven part remains D52's sup-norm gap of >= 2.30x.
Also worth reporting independently of any tiling question: at `eta=1`, `N_b=64`,
  `poly_order=4` the shipped default loses ~1.7% of S/N for a 10%-duty pulse and ~5-8%
  for a 5%-duty one, to grid coarseness alone. Of the 1.7%, about 0.5 points is the
  tiling choice and ~0.75% is irreducible at this `eta`.
Open questions: the pruning interaction is still open (this is geometry plus a profile
  model, and still does not model whether the closer leaf survives thresholding), and
  the Chebyshev/circular transforms.
Next session starts at: the human's read of the report, or the injection test.

## 2026-09-17 (u) — the Chebyshev port: the basis question, measured
Done:
  - **D58 — the port is sound for the same structural reason.**
    `poly_chebyshev_branch_batch` also goes through `branch_param_padded` with
    `(param_cur, dparam_cur, dparam_new)`, and `shift_cheby_errors` propagates the error
    vector with no reference to the values, so every leaf at a stage still shares
    `dparam_cur` and child offsets stay parent-independent. The leaf set is again exactly
    a Minkowski sum. `nearest_template_cheby.py`; brute-force agreement at rel 1e-12 on
    all 8 enumerable cases across the three strategies. Two convention differences:
    `poly_cheb_step_vec` is **uniform across axes** and the branch guard uses
    `|dparam_old - dparam_new|` rather than the width. The Chebyshev branch also
    transforms *inside* itself (coord_prev -> coord_cur) where the Taylor path branches
    then transforms.
  - **D59 — the Chebyshev nominal corner is `k_max/2`, exactly.** Uniform step, so the
    corner is `k_max * (1/2)`: 1.0 / 1.5 / 2.0 / 2.5 ... against the Taylor grid's
    `2^(k_max-1) - 1/2` = 1.5 / 3.5 / 7.5 / 15.5 (D46). **Linear against geometric**, and
    it confirms the conjecture relayed from the Chebyshev branch, whose 2.0 at order 4
    and 2.5 at order 5 are the `k_max=4,5` cases. The mechanism is *not* `|T_k| <= 1` as
    such: the Taylor grid with coarsening removed is also 1/2 per order, so the entire
    difference is the `2^(k-1)` coarsening factor.
  - **D60 — measured, the Chebyshev basis helps, and the tiling matters LESS there, not
    more.** Same 90 cells, `aggressive` EXACT in all of them:

    | basis | `aggressive` nearest | `quadrature` closer | proven gain | loss at ducy 0.10 |
    |---|---|---|---|---|
    | Taylor | 1.069 | 89/90 | >= 2.30x | 1.73% |
    | Chebyshev | **0.904** | 87/90 | >= 1.56x | **2.03%** |

    So Chebyshev puts a ~15% closer template near the signal *and* shrinks the tiling
    gap. This is the opposite of what the relayed EP-score result (+2.028 Chebyshev vs
    +0.584 Taylor for `aggressive`'s cost) would suggest, and the two are not
    contradictory: that measurement is a full pruning run with thresholds, mine is
    geometry. If template proximity is fine in Chebyshev but the score deficit is larger,
    the deficit is in the **pruning interaction** -- survival of the covering leaf -- not
    in template distance. Offered as a hypothesis; I have not run their pipeline.
  - **D61 — a lower sup-norm does not mean a lower loss, again.** Chebyshev wins on
    sup-norm (0.904 vs 1.069) and *loses* on amplitude (2.03% vs 1.73%). The Chebyshev
    residual is a combination of `T_k` over the whole domain and oscillates across it,
    while the Taylor residual peaks only briefly at the window edge (D56). What smears
    the profile is the residual's distribution. Third time this distinction has changed
    a conclusion; the nominal corner (D59: 2.0 vs 7.5, a 3.75x advantage to Chebyshev)
    predicts the amplitude ordering **wrongly**.
  - **D62 — `aggressive` is the only strategy that fits the shipped `branch_max`.**
    **AMENDED — the original wording was wrong in the one cell that matters, and it
    contradicted my own table.** The guard at `psr_utils.branch_param_padded:363` is
    strict (`num_points > branch_max`), and `branch_max` defaults to 16, so 16 builds and
    17 raises (verified directly). Max per-axis child count:

    | basis | `aggressive` | `quadrature` | `conservative` |
    |---|---|---|---|
    | Taylor | 7 ok | **28 raises** | **27 raises** |
    | Chebyshev | 9 ok | **16 ok, exactly at the limit** | **20 raises** |

    So three of the four non-`aggressive` configurations are unreachable at the default,
    and the D52 Taylor gain does require raising `branch_max`. But **Chebyshev +
    `quadrature` builds at the shipped default**, which makes it the only configuration
    with a proven gain (D60's >= 1.56x in 87/90 cells) that needs no config change. That
    is the actionable cell, and my first wording ("the others describe trees the code
    would refuse to build") denied it.
    The caveat attaches to the configuration, not the claim: 16 of 16 is *at* the limit
    with no headroom, and `num_points = ceil(dparam_cur/dparam_new)` moves with
    `poly_order`, `tobs`, `eta` and the segment count, so a nearby configuration turns it
    into a raised `ValueError`. It also reconciles the relayed report that
    Chebyshev+quadrature failed the guard: that run was `branch_max=128`, `ref_seg=3`, a
    different configuration, so a count above 16 there is consistent with exactly 16
    here. The boundary is configuration-dependent, which is the real point.
Conventions fixed: `min_excursion` now takes `basis_fn`, so one search serves both bases;
  `amplitude_loss` likewise. No monkeypatching.
Open questions: the pruning interaction, now the prime suspect for the Chebyshev score
  deficit (D60) and the last unmeasured mechanism; the circular basis.
Next session starts at: the human's read of the report.

## 2026-09-17 (v) — the report's figures are generated, not typed
Done:
  - **D63 — `report_numbers.py` + `report_numbers.json` + `tests/test_report_numbers.py`.**
    Every figure in `04_upstream_report.md` is now computed, committed as JSON, and bound
    to the document by tests. Aimed at this branch's actual failure mode rather than at
    typos: twice a sentence contradicted a table in the same document that was itself
    correct. So the suite (a) recomputes the cheap groups live against the JSON, (b)
    asserts each headline figure the report quotes equals the JSON value at the printed
    precision, and (c) checks the `branch_max` verdicts as **booleans** against the
    report's wording. `test_only_one_nonaggressive_configuration_builds` pins the D62
    amendment so it cannot be over-generalised again.
  - **D64 — it immediately caught a real inconsistency in the report.** The report quoted
    losses of 1.73% and 0.747% beside an advantage of 0.47%, and those do not combine:
    the ratio of the two medians is **+1.00%**, while the median of the per-cell ratios
    is **+0.47%**. The two losses are not co-monotone across cells, so the paired
    statistic is the smaller one, and it is the right one (cell-to-cell scatter cancels).
    A reader dividing the quoted medians would have got double the stated figure and
    concluded one of them was wrong. Both statistics are now in the JSON, the report
    leads with the paired one and says explicitly why division does not reproduce it,
    and a test asserts they still differ so the caveat does not become noise.
  - **D65 — and a typed figure that disagreed with the computation.** I had written
    "+0.53% Chebyshev" where the computed paired advantage formats to **0.52%**. Caught
    by the binding test on the first run. Corrected.
  - Two test bugs found in the writing, both mine, neither a defect in the report:
    matching `"7 — raises"` as a bare substring also matches inside `"27 — raises"` (now
    guarded against a preceding digit), and `_dig` split paths on `.` while several keys
    contain one (`ducy_0.10`), silently walking into missing keys. Now `/`-separated.
  - Full regeneration reproduces every live figure exactly: Taylor `aggressive` 1.06888
    (range 0.2399-2.114554, exact 90/90), `quadrature` 89/1/0 at gain 2.302674,
    `conservative` 28/0/62; Chebyshev `aggressive` 0.904452 (exact 90/90), `quadrature`
    87/3/0 at gain 1.563983.
Conventions fixed: no figure in the report is typed. `python report_numbers.py --full`
  regenerates the expensive groups (~12 min), `--only amplitude|nearest` regenerates one.
  Any paired comparison records both the per-cell median and the ratio of medians,
  because quoting one beside the inputs of the other is what D64 was.
Open questions: unchanged -- the pruning interaction and the circular basis.
Next session starts at: the human's decision on the injection campaign (relayed as
  approved in another session; not acted on, since it was not approved to me) and, if it
  goes ahead, the design doc and the power calculation before any production run.

## 2026-09-17 (w) — injection campaign: power calculation and pre-registered design
Done:
  - Campaign approved by assaferan directly (2026-09-17). Started with the power
    calculation rather than a run, since it could have cancelled the campaign.
  - **D66 — the first power calculation was wrong, and its own absurdity caught it.** It
    said `quadrature` recovers *nothing* (`P_d` 0.0000 against `aggressive`'s 0.0497 at
    S/N 10) despite having lower loss. Cause: the per-stage losses feeding it came from a
    budget-limited search, which returns an **upper bound**, and at 16 of 64 stages the
    bound was loose -- spiking to 44-60%. Survival takes a min over 64 stages, so one
    loose stage annihilated the arm. Raising the budget collapsed the spikes
    (200k -> 4M nodes: 44.75 -> 0.74, 60.74 -> 0.45, 47.68 -> 0.66), and a 16M check held
    two of three stages exactly with the third still drifting *down*.
    **The lesson generalises: a one-sided bound is harmless when comparing two numbers
    and fatal when fed into a model that takes an extremum over many of them.** Inputs to
    a model need the same convergence discipline as published measurements.
  - **D67 — inverted the calculation to report a detectable effect size**, instead of
    chasing per-stage convergence at all 63 stages (19 were still unresolved at 4M nodes
    after ~16 min). Power is computed against the loss ratio `r = loss_quad/loss_agg`,
    with `r` measured at **0.34 (IQR 0.21-0.45)** on the 41 converged stages. This is
    both cheaper and more robust, and it shows the conclusion holds across the whole
    plausible `r` range.
  - **D68 — the campaign is affordable, and pairing is what makes it so.** At the
    measured `r`, `dP_d = +1.3` points with a 1.3% discordance rate, peaking at
    **S/N 14**, needing **~290 paired injections for 80% power** (McNemar, alpha=0.05),
    and under 600 even at `r = 0.70`. My earlier guess of tens of thousands was for an
    unpaired design. Sanity check: the model reproduces `P_d` 0.05/0.18/0.45 at S/N
    10/12/15 against the Phase 3 injections' 0/3, 1/3, 3/3.
  - **D69 — Arm A, the only shipped-defaults pair, is the harder test and needs ~600.**
    Chebyshev `aggressive` vs `quadrature` is the sole pair that builds at
    `branch_max=16` (D62), but its loss ratio is nearer 0.7, so it needs roughly twice
    Arm B's injections. Worth knowing before choosing the arm rather than after.
  - `docs/metric_gridding/05_injection_design.md`: the pre-registered design --
    hypothesis, arms, pairing, stratification with a **null stratum**, positive and
    negative controls, threshold recalibration, the power table, and a pre-committed
    analysis with a fixed-n stopping rule. Written before any outcome was looked at.
Conventions fixed: any model input derived from a bounded search must carry a
  convergence check; report the detectable effect size, not just a single n.
Blocked on: `quadrature` and both Chebyshev threshold ladders are not calibrated (only
  `aggressive`, `conservative`, `metric` are cached), and the pipeline must be able to
  reproduce a noise realisation across arms for the pairing to work. Both are
  prerequisites for the first production run.
Open questions: unchanged, plus whether the noise realisation is reproducible across arms.
Next session starts at: the two prerequisites above -- threshold calibration for the
  chosen arm, and verifying paired noise -- then the run.

## 2026-09-17 (x) — stage weighting: the gain sits where detection is not decided
Done:
  - **D70 — the headline gain is not the detection-relevant gain.** Raised by the
    Overview session from the Injection Campaign session's per-stage `P_d` curve, and
    reproduced inside my own survival model: at S/N 14, **99% of first threshold failures
    occur in stages 1-10** (1% in 11-20, 0% beyond). **NOT independent confirmation** --
    corrected per that session: my model shares their Viterbi ladder and the same
    `mu(s) = snr sqrt((s+1)/nseg)` accumulation, so the agreement is a consistency check
    on a shared assumption set, not a second line of evidence. I had described it as
    independent, and it is not. My D52/D60 cells
    were `range(4, 63, 4)`, so only **2 of 15** sit in that window. Splitting the exact
    comparison:

    | basis | window | `aggressive` nearest | closer | proven gain |
    |---|---|---|---|---|
    | Taylor | 1-10 (99% of losses) | 0.731 | 47/60 | **1.71x** |
    | Taylor | 32-60 (~0%) | 1.030 | 48/48 | 2.61x |
    | Chebyshev | 1-10 | 0.582 | 47/60 | **1.15x** |
    | Chebyshev | 32-60 | 0.910 | 48/48 | 1.72x |

    So the advantage is roughly **half** the headline where it can matter, and less
    consistent: `quadrature` fails to beat `aggressive` in 13 of 60 early cells against
    0 of 48 late. `aggressive` is also already inside one tolerance early (0.731 / 0.582),
    so there is less to win. **Nothing in D52/D60 is retracted** -- those medians are
    correct over the cells they average -- but the framing was wrong, and the report now
    carries the split with the early figure marked as the detection-relevant one.
    Recorded in `report_numbers.json` as `stage_split_*` so it stays bound to the report.
  - **D71 — my own power numbers (D68/D69) were optimistic for the same reason, and are
    revised.** `r = 0.34` is an all-stage median; on stages 1-10 the measured ratio is
    **0.675** and `aggressive`'s loss there is 1.03% rather than 2.05%. Applying the
    ratio only where it matters: `dP_d = +0.50` points at S/N 14 with 0.50% discordance,
    needing **~770 pairs** rather than ~290. Arm A (Chebyshev) is weaker again (early
    gain 1.15x) and needs more still.
  - **D72 — and ~770 is a floor.** A peer session measures a stochastic term from the two
    arms scoring *different templates*, sd 0.133 against a deterministic amplitude
    difference of 0.044 -- noise ~3x signal -- which pairing on the noise realisation
    does **not** remove. My model omits it entirely, treating the arms as differing only
    by a deterministic `loss_s`. Carried into McNemar it shrinks the asymmetry and could
    raise `n` by an order of magnitude. It must be in the model before any `n` is
    committed.
  - Also worth having, from the same session and not mine to duplicate: the cached
    ladders here are Taylor and have no `quadrature` entry; regenerated on
    `poly_chebyshev_moving` both arms reach `P_d = 0.1031` with top thresholds **7.10 and
    7.00**, far closer than the Taylor 7.70/9.10 this log has been citing as the gap that
    would swamp a coverage effect. In Chebyshev that objection is much weaker. And
    `prune_on_overload_func` (`world_tree.py:527-551`) ratchets the threshold via
    `max(current, topk, median)` where `topk` depends on buffer size, silently -- a
    `max_sugg` confound that session found and owns.
Pattern worth naming: this is the fourth time a quantity has been right and its *framing*
  wrong -- peak vs distribution (D56, D61), prose vs table (D64, D65), and now aggregate
  vs stage-weighted. The measurements keep surviving review; what they *mean* keeps not.
  Every future headline needs "and where does this quantity matter?" answered in the same
  breath as the number.
Open questions: the template-noise term (D72) before any campaign `n`; Chebyshev ladder
  power; the `max_sugg` confound (owned elsewhere).
Next session starts at: folding D72 into `injection_power.py`, which decides whether the
  campaign is affordable at all.

## 2026-09-17 (y) — the early window is exact, and the floor gets worse
Done:
  - **D73 — the early-stage losses are EXACT, which closes the largest open caveat where
    it matters.** The Injection Campaign session flagged that the deterministic term
    `mu (L_A - L_B)` rests on both losses being *upper* bounds from a sup-norm search,
    that the two bounds need not be equally loose, and that if the asymmetry is
    stage-dependent it would bite precisely in the early window my stage split relies on.
    Checked directly: over stages 1-10, both arms' per-stage losses are **identical at
    4M and 16M nodes** (max change 0.0000 points), with the searches completing in 0-2 s.
    They are exhaustive, not bound-limited, so no asymmetry can enter there. The caveat
    stands for the late stages and is now retired for the early ones.
    Per-stage, stages 1-10 (%): `aggressive` 0.674 1.027 2.163 0.891 1.042 0.639 1.302
    1.657 0.896 1.467; `quadrature` 0.674 1.667 1.694 0.167 0.214 0.463 0.442 2.067
    0.362 0.918.
  - **D74 — and it shows `quadrature` is genuinely WORSE at some early stages**, not just
    unproven: 1.667 against 1.027 at stage 2, and 2.067 against 1.657 at stage 8. Since
    these losses are exact, that is a real sign reversal rather than a loose bound. It
    matches the sup-norm picture (`quadrature` fails to beat `aggressive` in 13 of 60
    early cells against 0 of 48 late) and it means the deterministic term is **not
    sign-definite in the window that decides detection**. A paired design must therefore
    expect discordance in both directions early, which is exactly what McNemar tests and
    exactly what an aggregate `P_d` comparison would hide.
  - **D75 — D72's stochastic floor is worse than recorded: sd 0.230, not 0.133.** The
    figure is `sqrt(2(1 - rho_AB))`, and it was computed from the profile-overlap `rho`
    that the same session has since corrected; with the score correlation it is
    `sqrt(2 x 0.0265) = 0.230`, 1.73x larger. Against the deterministic median of 0.044
    the noise is now ~5x the signal, not 3x. Through their model `pi` goes 0.697 -> 0.628
    and pairs 161 -> 391 at their `p_disc`. My ~770 omits the term altogether, so it
    remains a floor and the true `n` is materially higher.
Conventions fixed: do not call a check "independent" when it shares a ladder and a score
  model with the thing it agrees with (see the D70 amendment above).
Open questions: the campaign's true `n` once D75's term is in the model; whether the late
  stages' bound asymmetry matters at all now that the early window is exact (probably
  not, since the late stages decide ~nothing).
Next session starts at: assaferan's decision on who owns the power calculation. The
  Injection Campaign session declines to take it on a peer's say-so, correctly, and its
  `injection_power.py:power()` + `rho_check.py:propagate()` already implement the
  aggregation my model lacks. One owner beats two half-models.

## 2026-09-17 (z) — the model checked against outside data, and it runs ~1 S/N optimistic
Done:
  - **D76 — my sign worry (D74) was unfounded as a criticism of their model.** Checked by
    the Injection Campaign session against the code rather than memory:
    `injection_power.power()` forms `delta = mu (l_a - l_b)` with whatever sign it has and
    puts it through signed folded-normal moments
    (`E[D+] = sd phi(d) + mean Phi(d)`, `E[|D|] = 2 sd phi(d) + mean(2 Phi(d) - 1)`), so a
    stage where `quadrature` is worse contributes a small `E[D+]` against a full `E[|D|]`
    and correctly drags `pi` down. Reversals are absorbed, not violated. D74's *finding*
    stands -- the deterministic term is not sign-definite early -- but it breaks nothing,
    and it is one of the things the paired/McNemar framing was chosen for.
  - **D77 — the first non-circular check on my survival model: not contradicted, and
    biased.** (Wording corrected: an earlier draft said the model "passes" and
    "survives". At n = 24 that test has very little power, so non-rejection is an
    *absence of contradiction*, not support, and the stronger word reads as evidence it
    is not. Flagged by the Injection Campaign session.) Their pilot batch A is external to both models: **3 of 24 recovered at
    S/N 12** on real data, against the ladder's nominal `P_d = 0.1031`. My model gives
    0.1771 at S/N 12, and `P(X <= 3 | p = 0.1771) = 0.364` -- **not contradicted**, which
    at n = 24 is nearly all such a test can say.
    But the model reaches the nominal 0.1031 at **S/N 11.00** while the pilot is
    consistent with nominal at S/N 12, so it runs about **1.0 S/N optimistic**. Weak
    (n = 24) but it is the only check on the ladder from outside either model, and it is
    the one thing that is not circular in the way D70's amendment describes.
  - **D78 — WITHDRAWN by D86. Was: inject at S/N 15-17, not 14.** The model's
    discordance peak is flat over S/N 14-16 (0.50-0.53%, 723-772 pairs); shifted by the
    ~1 S/N bias that becomes **S/N 15-17**. Injecting at the model's nominal peak would
    sit below the real one, where the ladder is cleared less often than the model
    believes and the effect is smaller.
  - Also inherited, and marked as theirs: the `rho` correction propagates into
    `injection_power.strata()`, which builds per-position `pi` from the same substitution,
    so **the stratification boundaries move too** -- corrected `pi` ~0.05 lower per
    position, making the null stratum easier to find and the effect stratum harder. That
    is a change to my design doc's section 5, not just to a number, and I am waiting for
    their table rather than guessing at it.
Conventions fixed: the ladder is the shared dependency between my survival model and
  their `pi` profile, so any check that also assumes it cannot test it. Pilot batch A is
  currently the only external constraint; prefer it over internal agreement when they
  disagree.
Open questions: the campaign's `n` with D75's term in the model (theirs); the corrected
  stratum boundaries (theirs); whether the ~1 S/N bias is real or small-sample.
Next session starts at: assaferan's ruling on power-calculation ownership, and their
  strata table.

## 2026-09-17 (aa) — the effect stratum does not exist
Done:
  - **D79 — inherited from the Injection Campaign session, unverified by me: the `effect`
    stratum is empty.** Corrected decision-weighted `pi` over the same 24 positions and
    stages [2,6,10,14,20,28]:

    | `rho_AB` from | min | q25 | median | q75 | max | pi<0.55 | pi>0.70 |
    |---|---|---|---|---|---|---|---|
    | profile (was in use) | 0.440 | 0.595 | 0.645 | 0.704 | 0.741 | 12% | 29% |
    | corrected | 0.470 | 0.556 | **0.590** | 0.619 | **0.643** | 21% | **0%** |

    Their control: the profile row reproduces the cached per-position values to 5e-5 on
    all 24 positions, so the shift is the correction and nothing else. Consequences run in
    opposite directions -- the `pi > 0.70` **effect stratum is unreachable** (corrected
    max 0.643), collapsing their three-stratum design to two with the primary test on the
    mid stratum alone; while the **null stratum roughly doubles** in incidence, from ~1
    position in 8 to ~1 in 5, which makes the design's most important control cheaper to
    find. The distribution also narrows, which is what a larger common noise term does:
    it pulls every position toward 1/2.
  - **D80 — every correction since the campaign was approved has pushed `n` up, and this
    one continues it.** The primary test now runs where `pi` is ~0.590 rather than the
    0.628 that gave their 391 pairs, and pairs scale as `(pi - 1/2)^-2`, so this is
    another factor of ~2 on its own. Chain so far: stage weighting halved the gain where
    it matters (D70/D71, ~290 -> ~770 on my model), the `rho` correction raised the
    stochastic term 0.133 -> 0.230 (D75, ~5x the deterministic difference), the ~1 S/N
    optimism shifts the operating point (D77/D78), and now the effect stratum vanishes.
    **CORRECTED (D81): that chain double-counts.** The `rho` correction and the empty
    effect stratum are not independent -- they are the *same* shift in the *same*
    quantity, seen once per cell (the stochastic term) and once per position (`pi` going
    0.645 -> 0.593, and 0.697 -> 0.628 in aggregate). Listing them as two multiplicative
    factors inflates the chain. There are **two** distinct corrections here, not four:
    the stage weighting, and the `rho` correction with its two manifestations.
    I am deliberately **not** multiplying these into a single `n`: mine omits the
    stochastic term and theirs includes it, so the factors are not composable, and
    manufacturing a combined figure is exactly the kind of move that produced D64.
    The direction is unambiguous and the magnitude is theirs to compute.
  - D78's S/N 15-17 is compatible with their section 9, which already required the
    operating point to be fixed by a 3-point pilot **at the final `max_sugg`** rather than
    carried over (S/N 14 gave 13/20 and 15/20 at 2^18 against 7/20 and 4/20 at 2^14). My
    1.0 S/N offset is a reason to centre that pilot at **16** rather than 14.
  - Language corrected throughout per their caution: "not contradicted", never "passes"
    or "survives", for a non-rejection at n = 24. See the D77 amendment.
Conventions fixed: a non-rejection from a low-powered test is an absence of contradiction
  and must be worded as one. Do not compose `n` factors across two models with different
  terms in them.
Open questions: the campaign's `n` under the two-stratum design (theirs); whether a
  second pilot S/N point gets funded, which would test whether my model is shifted or
  differently shaped.
Next session starts at: assaferan on ownership, and on whether the campaign is still
  affordable given D80.

## 2026-09-17 (ab) — affordability computed: the campaign is viable
Done:
  - **D81 — my "chain of four corrections" double-counted, and I had put it to
    assaferan in that form.** Caught by the Injection Campaign session: the `rho`
    correction (D75) and the empty effect stratum (D79) are the same shift in the same
    quantity, observed once per cell as the stochastic term and once per position as
    `pi` (0.645 -> 0.593; 0.697 -> 0.628 aggregate). Counting both as multiplicative
    pushes inflates the picture. The distinct corrections are **two**: the stage
    weighting (D70/D71) and the `rho` correction with its two manifestations. In their
    numbers it appears once, correctly. My chain table overstated the erosion, and I
    escalated "is this still affordable?" partly on the strength of it.
  - **D82 — inherited, unverified by me: the campaign is affordable at the operating
    point the design actually picks.** 954 pairs, 13.2 core-hours at 50 s/run, against
    the 28 core-hours the original verdict quoted:

    | `p_disc` | mid-stratum pairs | total pairs | core-hours |
    |---|---|---|---|
    | **0.30** (measured at the operating point) | **755** | **954** | **13.2** |
    | 0.20 | 1132 | 1430 | 19.9 |
    | 0.15 | 1509 | 1906 | 26.5 |
    | 0.10 | 2263 | 2859 | 39.7 |
    | 0.08 (batch A, S/N 12) | 2828 | 3572 | 49.6 |

  - **D83 — and my `pi` was right by luck, with a term missing.** I used 0.590, the
    median over all 24 positions; the primary test runs on the mid stratum, which
    excludes the 5 positions below 0.55, so the correct value is **0.5927** -- dropping 5
    of 24 low values barely moves the median, so my factor of ~2 survived. What I missed
    is that the null stratum's 21% of screened positions never enter the primary test, so
    the pair count inflates by **1/0.79** on top. Right answer, wrong route, and one term
    short.
  - **D84 — powering against `p_disc = 0.08` would have been a category error**, and it
    is the one I was drifting toward by quoting the chain's worst end. That figure comes
    from batch A at S/N 12, where both arms recovered 3/24 -- far below the operating
    point, in a regime the design deliberately avoids. Discordance is maximised at the
    steepest part of the recovery curve, which is why the design fixes the operating
    point at on-grid `P_d ~ 0.5` in both arms, and there the pilot measured 0.30.
    Requiring power at 0.08 is requiring power in a configuration the design excludes by
    construction.
  - What *has* changed: `n = 2000` is no longer comfortable across the whole measured
    range (it covers `p_disc >= 0.145`), so their section 9 now makes the 3-point
    operating-point pilot at the final `max_sugg` a **precondition** rather than a
    formality -- if it returns `p_disc < 0.145`, `n` goes to ~3600 (50 core-hours) or the
    design is re-scoped. Theirs, as of their `db2a543`.
Conventions fixed: before presenting corrections as a chain, check they are
  *independent*; two manifestations of one shift are one correction. And power against
  the operating point the design selects, not against the worst value observed anywhere.
Open questions: the 3-point pilot result, which now gates `n`; power-calculation
  ownership (still assaferan's).
Next session starts at: their paired `max_sugg` result and the operating-point pilot.

## 2026-09-17 (ac) — D77/D78 withdrawn: the model is differently shaped, not shifted
Done:
  - **D85 — my own falsifiable prediction failed, and their new data is what tested it.**
    I predicted the pilot offset would be "roughly constant in S/N over 12-16", so that a
    second point would sit ~1 S/N left of my curve just as the first did. Checked against
    their paired run (`aggressive`, 50 realisations, S/N 14, 2^14 buffer):

    | | model | measured | |
    |---|---|---|---|
    | S/N 12 | 0.177 | 0.125 (batch A, 3/24) | model **high** by +0.052 |
    | S/N 14 | 0.368 | 0.600 (30/50) | model **low** by -0.232 |

    The discrepancy **changes sign**, and the measured slope over S/N 12-14 is 0.237 per
    unit against my model's 0.095 -- **2.0x steeper**. So no constant offset describes it:
    the real recovery curve is differently *shaped*, which is exactly the alternative I
    said one more point would distinguish. It did, and against me.
  - **D86 — D77's "~1 S/N optimistic" and D78's "inject at S/N 15-17" are withdrawn.**
    Both were an offset fitted to a single point (batch A, n=24), which cannot separate
    offset from shape. Their §9 criterion stands instead: fix the operating point where
    the *measured* on-grid `P_d` is ~0.5, because discordance peaks at the steepest part
    of the recovery curve. At S/N 14 and 2^18 they measure `P_d` = 0.76 and 0.80, so the
    3-point pilot should bracket **below** S/N 14, not above it. A direct measurement at
    the final buffer beats my survival model plus an inferred offset, and I am not going
    to defend the latter. One caveat I have asked them to confirm: I assumed batch A ran
    at 2^14; if it did not, the S/N 12 row above is not like-for-like. It does not change
    the conclusion, because the S/N 14 row alone (model 0.368 against 0.600 measured at
    2^14, 0.76 at 2^18) already shows the model is too pessimistic at the final buffer
    for its discordance peak to locate an operating point.
  - **D87 — inherited: the `max_sugg` pressure is arm-dependent and measured, and it
    blocks the run at 2^18.** Saturation > 0.9 in 1/50 (`aggressive`) against 10/50
    (`quadrature`) at 2^18, McNemar p = 0.012; at 2^14 both arms are pinned and the
    asymmetry is invisible (8/5, p = 0.58). Median saturation 0.165 against 0.724, median
    candidates 43k against 190k. Raising the buffer only ever helped: 20 flips, all one
    direction. The **outcome** bias is unresolved rather than null -- sign test
    `p = 0.34`, mean `Delta = +0.080`, 95% CI [-0.043, +0.203] against an effect of
    interest of 0.06, i.e. the point estimate **exceeds** the effect and points the same
    way. Bounding it below 0.06 needs n ~ 211 per cell.
  - **D88 — so my stage-split geometry is not currently measurable in a real search.** The
    arm the geometry favours (`quadrature`) is also the arm still sitting inside the
    candidate buffer at 2^18, and at n = 50 the two cannot be separated. That is not
    evidence the gain is unreal -- it is that this apparatus cannot see it yet. The fix is
    theirs and already stated: raise `max_sugg` until both arms clear 0.9 saturation.
    Worth recording plainly: the geometric result (D52/D60/D70) is unaffected, and its
    *relevance* now waits on an instrument change.
Conventions fixed: an offset inferred from one point is not an offset, it is a guess with
  one degree of freedom. Prefer a direct measurement at the final configuration over a
  model plus a fitted correction, and when the two disagree say which is measuring the
  real thing.
Open questions: batch A's buffer (asked); whether the outcome bias can be bounded below
  0.06; the operating point from their 3-point pilot, which now replaces D78.
Next session starts at: their operating-point pilot at the raised buffer.

## 2026-09-17 (ad) — the criteria never disagreed; and the rebase trap is benign
Done:
  - **D89 — the "unresolved disagreement" over the operating point dissolves: both
    criteria point to the same place.** Relayed by the Overview session as needing an
    owner, and worth settling from what each criterion is *for* rather than from the
    numbers. In my own model, discordance against S/N:

    | S/N | 13 | 14 | 15 | **16** | 17 | 18 | 20 |
    |---|---|---|---|---|---|---|---|
    | model `P_d` | 0.267 | 0.368 | 0.471 | **0.567** | 0.655 | 0.732 | 0.848 |
    | discordance | .00385 | .00495 | .00508 | **.00532** | .00476 | .00417 | .00273 |

    Peak discordance sits at **model `P_d` = 0.567**, i.e. exactly the `P_d ~ 0.5` their
    §9 targets. **The two criteria are the same criterion in different variables.** The
    disagreement was never conceptual -- it was entirely my `S/N <-> P_d` mapping being
    2x too shallow (D85), which misplaces *which S/N* delivers `P_d ~ 0.5` while getting
    the `P_d` at the optimum right.
    So their §9 wins for a sharper reason than "direct measurement beats model": it is
    parameterised in the **observable** rather than in the input, which makes it immune
    to precisely the error my model has. General form, worth keeping: fix an operating
    point by the quantity you can measure, not by the quantity you have to model.
  - **D90 — the rebase trap is benign, checked rather than assumed.** `upstream/main` is
    now `18d04b3` with #4, #6, #7 and #10 merged; this branch is 15 behind and carries
    three cherry-picks taken while those PRs were blocking (D28 among them).
    `git range-diff upstream/main...HEAD` reports all three as **identical**:
    `575b1f5 = 553f47d` (nested parallel), `8d42ab6 = f4ee619` (zero-filled states),
    `6d62d79 = 7864770` (trials_scheme -inf). So git will drop them silently on rebase
    and there is no `detection/` conflict to negotiate. The injection-design branch
    descends from here and inherits the clean rebase.
    **Not rebasing on my own initiative**: another session's branch is downstream of this
    one, so the timing is assaferan's to pick, not mine.
Conventions fixed: when two sessions appear to disagree about a design criterion, check
  whether they are the same criterion in different variables before either concedes.
Open questions: unchanged. The operating point now comes from their pilot at the raised
  buffer, under a criterion my own model agrees with.

## 2026-09-17 (ae) — the apparatus cannot answer the idealised question; and my external checks were void
Done:
  - **D91 — inherited, and it is decisive: `quadrature` never converges in the buffer.**
    From the Injection Campaign session, one fixed set of 10 realisations swept over four
    buffers, growth exponent `d log2(median ncand) / d log2(max_sugg)`:

    | arm | 2^14 | 2^18 | 2^19 | 2^20 | exponent |
    |---|---|---|---|---|---|
    | `aggressive` ncand | 1 446 | 16 778 | 16 778 | 16 778 | 0.00 (converged) |
    | `aggressive` sat | 0.088 | 0.064 | 0.032 | 0.016 | falling |
    | `quadrature` ncand | 10 906 | 152 007 | 390 575 | 808 116 | ~1.0 throughout |
    | `quadrature` sat | 0.666 | 0.580 | 0.745 | 0.771 | **rising** |

    `aggressive` converges at 2^18 and is identical to the unit at 2^19 and 2^20.
    `quadrature` never does: doubling the buffer doubles the count, because the ratchet
    relaxes the cut to keep the buffer full at any size. So **in the `quadrature` arm the
    scheme's thresholds have never been the operative cut**, at any buffer tested.
    Provisional on a 2^21 point they have running.
  - **D92 — CONSEQUENCE FOR MY OWN RECORD, which they did not claim and I am recording
    against myself: the external checks in D77 and D85 were void.** Both used recovery
    measured at **2^14**, where both arms saturate, so what those numbers measured was
    the ratchet, not the Viterbi ladder that `injection_power.survival` models. Therefore:
    - **D77's "not contradicted" is withdrawn as a validation.** Batch A's 3/24 at S/N 12
      was not a check on my survival model; it was a measurement of a different cut. I
      had called it "the first non-circular check" and "the only external constraint
      either of us has" -- it is neither, because it does not constrain the thing being
      modelled. The circularity I worried about (shared ladder) was the lesser problem;
      the datum was not measuring the ladder at all.
    - **D85's conclusion is restated.** The measured curve being 2.0x steeper than mine
      over S/N 12-14 is a real disagreement with the **buffer-limited** search, but
      "my model is differently shaped" was an unwarranted attribution -- the discrepancy
      may be entirely the ratchet. What stands is that my model does not describe the
      search as actually run; what does not stand is any claim about *why*.
    - **D86 and D89 are unaffected and D89 is strengthened.** Withdrawing D78 was right
      whatever the cause, and "parameterise the operating point by the observable, not by
      the quantity you model" is *more* robust in this light: measured `P_d ~ 0.5` does
      not care whether the operative cut is the ladder or the ratchet.
  - **D93 — where this leaves my result, stated once and plainly.** Nothing geometric is
    retracted: the exact search (D50), its brute-force validation, `aggressive` at 1.069
    exact (D51), the proven `>= 2.30x` / `>= 1.56x` (D52/D60), the stage split (D70) and
    the amplitude conversion (D55) all stand. What cannot be established at this
    configuration is whether that geometry **converts into detections**, because the
    thresholding that would convert it is not the thresholding the search applies in the
    arm the geometry favours. "Geometry proven, relevance not measurable with this
    apparatus" is the end state, and the apparatus -- not the geometry -- is what failed.
  - Their §7.1 split survives and only half dies: the *idealised thresholded search*
    question is unanswerable here, while "does tiling matter for what a user actually
    runs?" remains answerable at the shipped default, treating the ratchet as part of the
    system -- but that is a claim about pyloki-as-shipped, not about tiling, and must be
    reported as such.
Open question for assaferan, not for either session: a third path neither of us has
  costed -- shrink the configuration (fewer segments, smaller `prune_poly_order` or
  `branch_max`) until `quadrature`'s true candidate count is affordable, and ask the
  idealised question there. Different external validity, and whether a smaller search is
  still the search anyone cares about is a judgement call rather than a measurement.
Next session starts at: assaferan on the shrink-the-configuration option, and on the
  report, which is unaffected by all of this but should now carry D93's framing.

## 2026-09-17 (af) — a valid check exists after all, and it falsifies my survival model
Done:
  - **D94 — their caution was right in principle, wrong about D85, and led to a better
    check that kills the model.** They warned I had over-withdrawn: `aggressive`
    converges at 2^18 (D91), so anything resting on `aggressive`-at-2^18 still constrains
    the ladder and the ratchet explanation is unavailable for it. Checked how D85 was
    actually computed: **both** its measured points were at 2^14 (batch A's 3/24 at S/N
    12, and 30/50 at S/N 14), which saturates for both arms, so D85 as computed really is
    fully voided and their caution does not rescue it.
    But applying their principle properly surfaces a point I never used -- and it is the
    one that matters:

    | | value |
    |---|---|
    | `aggressive`, S/N 14, **2^18** (converged, sat 0.064) | **38/50 = 0.760** |
    | my survival model at S/N 14 | 0.368 |
    | `P(X >= 38 | p = 0.368)` | **1.8e-08** |

    So the model is **not merely undescriptive, it is decisively falsified** -- and in the
    *pessimistic* direction, by +0.392 in `P_d`. It would need **S/N 18.5** to reach the
    recovery the real search achieves at 14, i.e. it is about **4.5 S/N pessimistic** at
    the converged buffer. The measured 95% interval [0.640, 0.880] excludes it by a wide
    margin.
  - **D95 — so `injection_power.survival` is retired as a quantitative instrument.** Its
    absolute `P_d` values, its operating point and its pair counts all depend on the
    `S/N -> P_d` mapping that D94 falsifies, so D68's ~290, D71's ~770 and the `P_d`
    column of the power table are all unreliable. Their 954 pairs, computed on their own
    model against a *measured* `p_disc`, is the number to use. What survives of mine is
    the part brute force validated: the exact nearest-template search and the
    smearing-based amplitude conversion. The survival wrapper was always the weakest link
    -- a single-leaf Gaussian proxy with no beam and no competing candidates -- and it is
    now measured to be so.
  - Worth naming, because the sequence is instructive: I described this model as
    "optimistic by ~1 S/N" (D77), then as "differently shaped" (D85), and it is in fact
    **pessimistic by ~4.5 S/N**. Every one of those framings was derived from comparisons
    against buffer-saturated data. The lesson is not that I mis-estimated a bias three
    times -- it is that **three successive characterisations of an error were all computed
    against data that could not measure it**, and one valid datum settled it immediately.
    Check what the comparison datum is measuring before interpreting the residual.
  - **D96 — and over-withdrawing has a cost, which they were right to flag.** Had I left
    D85 as "no claim about why", the `aggressive`-at-2^18 constraint would have stayed
    unused and this falsification would not have been found. Retracting to the safe
    minimum is not free: it throws away real constraints along with the bad ones.
Conventions fixed: when voiding a result, void it at the granularity of the *datum*, not
  the result -- ask which measurements are still valid rather than discarding the whole
  comparison. And prefer one datum known to measure the modelled quantity over three that
  do not.
Open questions: unchanged. Nothing here touches the geometry (D93 stands).

## 2026-09-17 (ag) — why the model was pessimistic, and 954 is a range not a number
Done:
  - **D97 — their M6 diagnosis is confirmed, and it is the dominant term.** Their reading
    of D94: a single-leaf model being *pessimistic* is what you expect if real survival is
    a **max over many near-covering leaves** rather than one draw. Tested directly. The
    inconsistency in my model is now nameable: it applies thresholds **calibrated for the
    whole branching pattern** while giving the signal **exactly one leaf**, so false
    alarms get `N` leaves and the signal gets 1. Giving the signal `K` partially
    decorrelated leaves (shared data component plus a per-leaf part, since templates at
    different offsets fold differently) at S/N 14:

    | K | rho=0.99 | rho=0.95 | rho=0.90 |
    |---|---|---|---|
    | 1 | 0.363 | 0.366 | 0.366 |
    | 5 | 0.412 | 0.479 | 0.521 |
    | 20 | 0.444 | 0.546 | 0.617 |
    | 100 | 0.472 | 0.608 | **0.703** |

    against the measured **0.760** and the single-leaf **0.368**. So plausible `K` and
    `rho` recover most of the 4.5 S/N gap: M6 is not mildly conservative, it is the
    leading error. Note this is *not* the tiling multiplicity of D43 -- it counts
    near-covering leaves whose loss is slightly worse, which exist in both arms, and
    `aggressive` had 16 778 candidates alive at 2^18.
  - **D98 — and their 954 must not be promoted to authoritative; I did promote it and am
    correcting that.** They re-ran their aggregation under different stage weightings,
    changing nothing else:

    | stage weighting | pi | total pairs |
    |---|---|---|
    | ladder `succ_h1` (the quoted one) | 0.628 | 954 |
    | flat | 0.654 | 339 |
    | `sqrt(succ_h1)` | 0.642 | 404 |
    | late stages only | 0.680 | 246 |
    | early stages only | 0.554 | **2 842** |

    **A factor of 8**, larger than the `rho` correction, and the weighting rests on a
    modelled per-stage survival curve from the same single-leaf Gaussian family as the
    model D94 falsified. So `n` is 340-2 840 and 954 is one point in it. I told assaferan
    954 was "the number to use" on the strength of it being measured-`p_disc`-based; the
    `p_disc` is measured, the *weighting* is not, and I conflated the two.
  - **D99 — and my own result cannot break the tie, though it leans one way.** If the real
    search survives far better than modelled (D94), the true stage weights decay *less*
    steeply than `succ_h1`, which moves toward the flat/late rows and `n` down toward 340.
    But the claim that the decision is made early (D70, and their §4.3) comes off the same
    modelled profile and is the 2 842 corner. **Both readings derive from a curve now
    known to be wrong**, so neither can be asserted -- including my own stage-split
    framing, whose *weighting* by "where detection is decided" inherits this. D70's
    measured early/late gain split stands as geometry; what is now unsupported is the
    claim that the early window is where it matters, since that came from the falsified
    survival profile.
Conventions fixed: "computed from a measured quantity" does not make a result measured --
  check every input, not the headline one. That is what D98 got wrong.
Open questions: the real per-stage survival profile, which nothing on either branch
  measures and which would settle both `n` and D99. In principle instrumentable (record
  how many realisations still have a covering leaf alive at each prune level, not just
  end-to-end recovery), but that instrumentation does not exist and is not mine to start.

## 2026-09-17 (ah) — the per-stage survival profile, MEASURED: the decision is late
Done:
  - **D100 — `docs/metric_gridding/survival_profile.py`, and it is validated against an
    external rate.** `Pruning.execute()` is a loop over `execute_iter()` with survivors
    in `world_tree`, so the loop is driven directly and the survivor set inspected between
    levels. Nothing in `src/` is touched. Per level, leaves and truth are both expressed
    in that level's frame (`ref_cur`, `t_half`) and the residual is evaluated over the
    accumulated data; the recorded quantity is the **minimum phase excursion to the truth
    over all survivors**, not a binary.
    Validation: final-level "alive" fraction **29/40 = 0.725** against their independently
    measured recovery of **38/50 = 0.760** (`aggressive`, S/N 14, 2^18), two-sided
    binomial `p = 0.72`. The alive threshold is not assumed -- the final excursions are
    cleanly **bimodal** (12.01 against 66.89, a 5.6x gap) and `T = 28.3` is the geometric
    mean of the gap. Checked for robustness: at every `T` from 5 to 100 the earliest loss
    is at level >= 11, and at the calibrated `T` it is >= 13.
  - **D101 — the modelled survival curve is not merely miscalibrated, it is backwards.**
    Measured, 40 real pruning runs:

    | level | 1 | 10 | 20 | 30 | 40 | 50 | 63 |
    |---|---|---|---|---|---|---|---|
    | fraction alive | 1.000 | **1.000** | 0.900 | 0.850 | 0.825 | 0.725 | 0.725 |
    | median min excursion | 1.67 | 1.27 | 2.23 | 3.19 | 4.61 | 5.61 | 6.49 |

    First-loss levels: **13, 13, 17, 19, 24, 27, 36, 41, 41, 44, 47** — median 27, none
    before 13. So **0% of losses occur by level 10**, where the falsified model put
    **99%**, and 36% by level 20, 55% by 30, 100% by 50. The decision is made in the
    middle and late stages, not early.
  - **D102 — so D70's framing is reversed, and the original headline was right.** The gain
    measured in the window where covering leaves are *actually* lost (levels 13-47):

    | basis | early 1-10 | **decision 13-47** | late 32-60 | headline 4-60 |
    |---|---|---|---|---|
    | Taylor | 1.71x (47/60) | **2.47x (54/54)** | 2.61x (48/48) | >= 2.30x |
    | Chebyshev | 1.15x (47/60) | **1.57x (54/54)** | 1.72x (48/48) | >= 1.56x |

    In the decision window `quadrature` is closer in **54 of 54 cells** -- every one --
    and the gain is essentially the headline figure. So "the gain sits where detection is
    not decided" (D70, and what I told assaferan and both peer sessions) is **false**: it
    sits almost exactly where the decision is made. The headline `>= 2.30x` / `>= 1.56x`
    was approximately right all along, because stages 4-60 happens to cover the decision
    window better than my "early" window did.
  - **D103 — consequence for the campaign's `n`: the expensive corner is excluded.** The
    340-2 840 range (D98) was driven by the stage weighting; the early-only weighting
    (`pi = 0.554`, **2 842 pairs**) assumed the decision is early and is now falsified,
    while the flat (339) and late-only (246) weightings are the ones the measurement
    supports. So `n` sits at the **cheap** end, not the expensive one. Their per-stage
    `pi` row survives untouched (geometry plus corrected `rho`, no `reach` term); it is
    only the weighting that moves, and it moves favourably.
Corrections to my own record, in one place: D70's framing reversed (D102); D99's "neither
  reading can be asserted" is resolved in favour of the late reading; the caveat I added
  to the report -- that the early window's primacy was unsupported -- is replaced by a
  measurement showing it is wrong.
Caveats: `aggressive` arm only, one injection parameter set, N = 40, and the alive
  threshold is calibrated from the same runs (robust across 5-100, but not independent).
  A `quadrature` profile is the obvious next measurement and is blocked by the same
  buffer non-convergence that blocked the campaign (302 s/run at 2^21).

## 2026-09-17 (ai) — the profile hardened: independent criterion, second parameter set
Done:
  - **D104 — the threshold circularity is removed.** The campaign session's first caution
    on D101 was that the alive threshold was calibrated from the same runs. Re-derived the
    profile under an **independent** criterion instead: the metric mismatch `m <= 1.0`
    that Phase 3's recovery test already used (`m_recover`), evaluated in each level's own
    accumulated baseline. Nothing calibrated from the data.
  - **D105 — and the second caution, one parameter set, is addressed.** Two injections,
    N = 30 each: set A `(accel 1.0, jerk 0.05, snap 0.001)`, set B
    `(accel -3.2, jerk 0.18, snap -0.004)`.

    | | set A | set B | pooled |
    |---|---|---|---|
    | alive at final level | 15/30 | 17/30 | — |
    | losses **by level 10** | **0/15** | **0/13** | **0/28** |
    | losses by level 20 | 9/15 | 11/13 | 20/28 |
    | median loss level | 17 | 15 | 15 |
    | earliest loss | 13 | 13 | 13 |

    **Zero losses before level 13 under both criteria and both parameter sets.** So D101's
    decisive claim -- 0% by level 10 against the model's 99% -- is robust to the criterion
    and to the signal. What does move is the *centre* of the window: median loss level 15
    under the metric criterion against 27 under the excursion one, because `m <= 1.0` is
    stricter than end-to-end recovery (final alive 0.50-0.57 against the 0.725 that
    matched their 38/50). So the excursion threshold was better calibrated to recovery,
    and the metric criterion is better anchored; both are reported.
  - **D106 — the gain survives every decision-window definition, and the 100% does too.**
    Re-measured on windows drawn from the metric criterion's losses:

    | window | Taylor | Chebyshev |
    |---|---|---|
    | 13-47 (excursion crit.) | 2.47x (54/54) | 1.57x (54/54) |
    | 13-37 (metric crit.) | 2.09x (60/60) | 1.50x (60/60) |
    | 13-16 (earliest quartile) | 1.74x (24/24) | 1.54x (24/24) |
    | 1-10 (the falsified premise) | 1.71x (47/60) | 1.15x (47/60) |

    `quadrature` is closer in **100% of cells in every decision window**, against 78% in
    the falsified early window. Taylor's gain is window-dependent (1.74-2.47x) so it must
    be quoted as a range; Chebyshev's is stable at ~1.5x. D102's conclusion stands under
    all of them.
  - **Deliberately NOT done: recomputing their `DECISION_STAGES`.** Their window
    `[2, 6, 10, 14, 20, 28]` comes from the premise my profile reverses, and four of six
    stages sit at or below 14 while none exceeds 28 -- so their stratification and strata
    boundaries are built over stages where little is being decided. That is their model
    and their doc, and recomputing it here would put a new headline on someone else's
    unverified inputs, which is the failure mode this branch has repeated. Their decision
    to leave it recorded-but-not-recomputed is the right one and I am matching it.
Conventions fixed: when an inherited measurement blocks someone else's work, harden the
  measurement rather than applying it for them.
Open questions: a `quadrature` profile, still blocked by the buffer non-convergence
  (302 s/run at 2^21); their strata, if the configuration question is reopened.

## 2026-09-17 (aj) — the validation fails on more data; the shape survives
Done:
  - **D107 — their reconciliation query, answered: a denominator mismatch in my
    presentation, not a non-monotone predicate.** "Losses by level 20 = 20/28" was a
    fraction of the **28 runs that ever dip**, while the alive fractions were over all
    **60** runs. Alive at level 20 = (60 - 20)/60 = 0.667, which matches the computed
    instantaneous 0.667. The figures reconcile; the presentation put two denominators in
    one table. Their hypothesis (a) was reasonable and wrong; the query was right.
  - **D108 — and the predicate is *essentially* monotone, so it can serve as `reach(s)`.**
    Checked explicitly: the `m <= 1.0` predicate is re-satisfied after a dip in **3** of
    3 720 level-observations, the excursion predicate in **1**. So the never-dipped and
    instantaneous curves agree to three decimals and either is usable. Non-monotonicity is
    physically possible -- refinement can bring a child inside tolerance where no parent
    was -- but it is negligible here.
  - **D109 — PARTLY OVER-CORRECTED; see D111. The substance holds (one batch was never a
    validation) but the "discrepant" label was a statistical error.** Original entry:

    | | alive at final level |
    |---|---|
    | set A, batch 1 (N=40) | 29/40 = 0.725 |
    | set A, batch 2 (N=30) | 15/30 = 0.500 |
    | **combined set A (N=70)** | **44/70 = 0.629** |

    My own two batches are only marginally consistent with each other (Fisher exact
    `p = 0.080`), and combined against their independently measured 38/50 = 0.760 the
    two-sided `p = 0.019` -- **discrepant**, with 0.760 lying outside the combined 95%
    interval [0.514, 0.743]. I reported the batch-1 agreement as the validation; on more
    data it fails. That validation was, in their words, "the only thing tying your
    instrument to an independently measured quantity", so its loss matters.
    Candidate causes, none yet established: my injected parameters may differ from theirs
    (I have asked -- if the signals differ the rates were never comparable and the
    batch-1 agreement was coincidence); their recovery criterion may not be my alive
    predicate; or the instrument is biased. **Until resolved, the profile's absolute
    level is not validated and should not be used to set `reach(s)` values.**
  - **D110 — the SHAPE survives, and that is what the conclusions rest on.** "No losses
    before level 13" holds across **two criteria** (excursion, and metric `m <= 1.0`),
    **two parameter sets** (A and B), and **two batches**: 0/11, 0/15, 0/13 by level 10
    in each. The falsified model put 99% there. So D101/D102 -- the decision is not made
    early, and the gain sits in the decision window -- do not depend on the absolute
    calibration that D109 withdraws. What does depend on it is any numeric `reach(s)`, and
    hence their `pi`.
  - Their pushback on "better anchored" is **accepted and withdrawn**: the metric
    criterion is independently *defined* but its final alive fraction (0.50-0.57) is the
    one that sits furthest from their measured recovery, while the calibrated-threshold
    excursion version was closer. Being independent of the data is not the same as
    agreeing with external data. The two measure different things, `m <= 1.0` is stricter
    than end-to-end recovery, and per D109 neither is currently validated as `reach(s)`.
  - Also noting their 368 -> 368-691: they checked my window shift against their model
    rather than assuming stability, and it moves by 1.9x. Their defensible claim is now
    "the early-decision corner is excluded and `n` is in the high hundreds" -- which is
    exactly the shape-not-level split D110 describes, arrived at independently.
Conventions fixed: never put two denominators in one table. And an agreement on one batch
  is not a validation -- run the second batch before calling it one.

## 2026-09-17 (ak) — "discrepant" was my error: it is consistent, and merely unvalidated
Done:
  - **D111 — D109's "discrepant" is withdrawn. I tested an estimate as though it were a
    known parameter.** Their 38/50 is an *estimate* from n = 50, not a fixed rate, so the
    correct comparison is two-sample rather than one-sample:

    | test | result |
    |---|---|
    | one-sample, 44/70 against `p = 0.760` fixed (**what I did**) | `p = 0.0192` |
    | two-sample, 44/70 against 38/50 (**correct**) | Fisher exact `p = 0.1641` |

    So there is **no evidence my instrument disagrees with theirs**. Pooled, the best
    estimate is **82/120 = 0.683**, 95% interval [0.600, 0.767], which contains both my
    0.629 and their 0.760. The honest label is **unvalidated** -- because one batch
    agreeing was never evidence -- and *not* discrepant. D109's substance stands and its
    headline does not, which is the same shape of error as the framings this branch has
    repeatedly got wrong, this time in the pessimistic direction against my own work.
  - **D112 — their non-exchangeability point is verified in the code, and it widens the
    intervals further.** `calibrate_scale_on_folds` (`simulation/pulse.py:72-110`) solves
    iteratively for the injected amplitude so that the *measured* boxcar S/N on the
    *realised* folded noise equals `snr_target`, with `noise_scale` taken from that
    realisation's own off-pulse std (`pulse.py:381`). So the injected signal is tuned per
    draw and runs are not exchangeable Bernoulli trials in the way a two-proportion test
    assumes. Their documented spread at identical settings supports it in one place:
    `quadrature` 4/20 against 11/20 gives Fisher `p = 0.048`, while `aggressive` 7/20
    against 11/20 gives 0.341. Either way the effect is overdispersion, which makes every
    nominal `p` above **anti-conservative** -- so the corrected 0.164 is if anything an
    understatement of the consistency.
  - **D113 — one real residual difference remains, and it is theirs to name: end-to-end
    versus final-level.** Their 38/50 is a candidate in the final `ScatteredPeriodogram`
    within `m <= 1.0` of truth in the full-baseline metric; my predicate is evaluated at
    the last prune level, before the final ascend and the report/resolve steps. So the two
    need not agree exactly even on identical data. Cheap to test if it ever matters -- run
    both on the same realisations and compare per-run -- and currently not worth it, since
    nothing uses the absolute level.
  - Noted from their side, and it is the better catch: every line-number citation in their
    upstream note pointed at *this branch's* `prune.py`, which carries 166 lines this
    branch added -- a note whose purpose is to point a maintainer at specific code, in the
    wrong file's numbering. Now machine-checked against `upstream/main` and the diff
    copied from `git diff` rather than transcribed (the transcription had drifted, claiming
    "three lines" against an actual +16/-2). Their framing, which I am adopting: **my
    denominator mismatch and their wrong-file numbering are the same mistake -- a number
    correct in its own frame, presented in a frame where it is not.** Theirs would have
    reached a maintainer; mine stayed between us.
Conventions fixed: when comparing against someone else's measured rate, test two-sample
  -- their number has an interval too. And check whether the trials are exchangeable
  before using any binomial test on this pipeline; `calibrate_scale_on_folds` means they
  are not.
```
