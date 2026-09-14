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

## Still to do in Phase 2

- **Step 2 (`transforms.py`) was never implemented, and it — not step 6b — is now
  what keeps `"metric"` from running end-to-end.** `shift_taylor_errors` and
  `shift_taylor_full` still `raise ValueError(f"Invalid tiling strategy: metric")`.
  Measured after step 6b landed, with `tiling_strategy="metric"`, `poly_order=3` on the
  `tests/test_prune.py` fixture, a full `prune_dyp_tree` run now gets through **branch,
  validate, resolve, shift_add and score** on metric children and dies at the
  `transform` step, `prune.py:651 -> transforms.py:136`. `generate_branching_pattern`
  fails the same way, so the threshold scheme cannot be built for `"metric"` either.

  Step 2 is not mechanical, and the difficulty is worth stating before the session that
  does it. The plan says to "transform the Cholesky factor (or `g`) exactly with
  `transform_metric` and return the new axis extents for column 1. No inflation, no
  diagonal truncation." But under D12 a leaf stores only `ellipsoid_axis_extents`, the
  ellipsoid's **bounding box**, and a bounding box does not determine the ellipsoid — so
  the exact transform cannot be done from leaf state alone. The options look like:
  (a) rebuild `g` for the stage inside the transform, which needs `m_max` and the
  duty cycle on the structref and the interval from `coord_next`/`coord_cur`;
  (b) accept that column 1 is a diagnostic under `"metric"` (D12 already says nothing
  derives spacing from it) and transform it by the `"conservative"` AABB rule, which
  the plan's "no inflation" explicitly rules out; or (c) store the stage's `g` or its
  Cholesky factor somewhere the transform can reach. This wants its own session and its
  own decision entry.

- **Steps 4 and 5 are also still open**: the `resolve` diagnostic against the base grid
  (step 4) and the `validate`/`report`/`ascend`/`io/cands.py` audit (step 5).

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
```
