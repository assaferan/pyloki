# Phase 0 — Pipeline notes

Base commit: `6b11aba` (identical to `upstream/main`; includes merged PRs #1 and #2).
Branch `metric-gridding`, worktree `worktrees/metric-gridding` (a child of the repo root,
excluded from git via `.git/info/exclude` rather than a tracked `.gitignore` rule).

Every claim below was checked by reading the code at this commit. Plan facts that were
wrong or incomplete are flagged **CORRECTION** / **ADDITION**.

## 1. Verification of the plan's "Repository facts"

| plan fact | verdict |
|---|---|
| Leaf is `(n_leaves, poly_order + 2, 2)`; col 0 value, col 1 half-width | **confirmed** |
| Rows `[:-2]` = `[d_kmax..d_1]`, row `[-2]` = `d_0`, row `[-1]` = `f_0` / basis flag | **confirmed** |
| Branching per-axis via `branch_param_padded` + `cartesian_prod_padded` | **confirmed** |
| Steps from `poly_taylor_step_d_vec` -> `poly_taylor_step_f`, `use_cheby=True` | **confirmed** |
| Three tiling strategies in `transforms.py::shift_taylor_errors` | **confirmed** |
| `poly_taylor_transform_batch` calls `shift_taylor_full` | **confirmed** |
| DP plumbing in `dynamic/dyn_poly_taylor.py`; buffer in `utils/world_tree.py` | **confirmed** |
| `B(s)` / thresholds in `detection/schemes.py`, `detection/thresholding.py` | **confirmed** |
| Existing tests: `tests/cheby.py`, `mat_inv_tests.py`, `test_maths.py` | **confirmed at this commit** |

**ADDITION (leaf row layout, exact).** From `poly_taylor_seed` (`core/taylor.py:52-64`),
counting from the end: row `[-1]` = `f_0` (col 0) and basis flag (col 1); row `[-2]` =
`d_0`, never branched on (`leaves[:, -2, 0] = 0`, comment "we never branch on d0"); row
`[-3]` = `d_1` (velocity), seeded value 0 with `dparam = df * C_VAL / f0`; rows `[:-3]` =
`[d_kmax .. d_2]`. So `leaves[:, :-1, 0]` is the full Taylor vector
`[d_kmax, ..., d_1, d_0]` (`poly_order + 1` entries) and `leaves[:, :-2, *]` is the
`poly_order` branchable axes. Branching uses `n_params = poly_order`
(`core/taylor.py:106`).

**ADDITION (tests on other branches).** `tests/test_transforms.py` and
`tests/test_prune.py` exist on sibling branches (`Chebyshev`, `fix-prune-segfault`) but
not here. `test_prune.py` is the only end-to-end pruning coverage anywhere and it
depends on the `prune.py` SIGSEGV fix (PR #3, open upstream). **Phase 3 cannot run at
all on this base** — every `prune_dyp_tree` call segfaults until PR #3 lands or
`fix-prune-segfault` is merged in. Phases 0-2 are unaffected.

## 2. Consumers of leaf column 1 (the half-width)

| consumer | file:line | what it does with col 1 |
|---|---|---|
| `common.get_leaves` | `core/common.py:50` | writes initial `dparams[j]` into col 1 |
| `poly_taylor_seed` | `core/taylor.py:54,60,64` | reads `df` from `[-1,1]`; writes `d_1` half-width `df*C/f0`; sets basis flag `[-1,1]=0` |
| `poly_taylor_branch_batch` | `core/taylor.py:109` | reads `dparam_cur_batch = leaves[:, :-2, 1]` — the branching input |
| " | `core/taylor.py:112,165,168` | carries basis flag; writes `branched_dparams` back to col 1 |
| `poly_taylor_transform_batch` | `core/taylor.py:354` (via `shift_taylor_full`) | shifts values and errors together; **this is where tiling strategy acts** |
| `poly_taylor_report_batch` | `core/taylor.py:346,348,353,359` | gauge-transforms col 1 into reported parameter uncertainties |
| `periodogram.add_run` | `periodogram.py:249` | writes col 1 out as the `d<name>` columns of the candidate table |
| `io/cands.py` | `io/cands.py:343` | persists `leaves_report` (values + col 1) to HDF5 |
| `circular` branch/validate | `core/circular.py:152-156` | reads `dsnap/djerk/daccel` from col 1 (circular basis only) |
| `chebyshev` equivalents | `core/chebyshev.py:54,60,68,114,124,173,176,383-397` | same roles in the Chebyshev basis |

**CORRECTION to the plan's step-2 list.** `validate` is a **no-op** in the Taylor path:
`dyn_poly_taylor.py::validate_func` (line 465) returns `(leaves_batch, leaves_origins)`
unchanged, so it consumes column 1 not at all. Nothing to audit there for
`tiling_strategy="metric"`.

**`world_tree.py` does not interpret column 1.** It only moves whole leaf rows
(`leaves[pos] = leaf`, slice copies at lines 473/499/517/608) and reads **column 0** for
parameter values (lines 577, 589) and for its de-duplication keys (919, 937). So the
candidate buffer is agnostic to what column 1 means — good news for the plan's
Phase 2 decision (a) vs (b).

**`ascend` does not read column 1.** `ascend_func` (`dyn_poly_taylor.py:588`) re-resolves
using `leaves` values and rewrites `folds`/`scores` only.

So the binding consumers of column-1 *semantics* are exactly three: `branch` (reads it as
a per-axis half-width), `transform`/`shift_taylor_errors` (propagates it), and
`report`/`periodogram`/`cands` (publishes it as an uncertainty).

## 3. How `B(s)` and the thresholds are computed

- `PulsarSearchConfig.generate_branching_pattern(kind, ref_seg, use_cheby_coarsening=True)`
  (`config.py:758`) returns the per-stage branching factors. `kind` selects
  `generate_bp_poly_taylor{,_approx}` / `generate_bp_poly_chebyshev{,_approx}`.
- It reaches the step sizes through `get_dparams(niters_ffa, use_cheby_coarsening)` and
  `get_dparams_actual(...)` (`config.py:789-791`), i.e. `eta`, `nbins` and the
  `2**k` coarsening all enter here.
- `tiling_strategy` enters `generate_bp_*` via `self.tiling_strategy` (`config.py:753`,
  `817`) — so **`B(s)` already depends on the tiling strategy**, and adding `"metric"`
  will change the branching pattern and therefore require threshold recalibration
  (the plan's Phase 3 step 2).
- Thresholds: `detection/thresholding.py::determine_scheme(survive_probs, branching_pattern, ...)`
  derives a ladder from target per-stage survival; `evaluate_scheme(thresholds, ...)`
  only *scores* a given ladder. `StatesInfo.thresholds` extracts the numbers.
  `DynamicThresholdScheme` is the Viterbi optimiser referred to in Phase 3 step 2.

## 4. Current gridding invariant (as implemented)

For any true parameter vector `λ*` inside a leaf's axis-aligned box, the intended
guarantee is
`max_t |Φ(t; λ*) − Φ(t; λ_centre)| ≤ η / N_b`
over the current accumulated interval — a **sup-norm, worst-instant** criterion. It is
enforced per axis and independently: `branch_param_padded` subdivides each axis until its
own step is below the target, and the box is the Cartesian product of those intervals.

**Two measured caveats on that invariant, from work on the sibling `Chebyshev` branch.**
Both bear directly on this project's premise and are worth carrying in:

1. The `2**k` Chebyshev coarsening in `poly_taylor_step_f` means the box does **not**
   actually honour `η`. Worst-case per-cell sup-norm phase error is `1.5x` the nominal
   `η/N_b` at `poly_order=2`, rising to `7.5x` at 4 and `15.5x` at 5; typical (random
   in-cell position) error is `1.5x/3.1x/6.1x` at `poly_order` 3/4/5. With the factor
   removed the numbers land exactly on the Chebyshev-basis values. So the existing
   invariant is loose by a factor growing exponentially in `poly_order`.
2. The tiling dilemma the plan sets out to dissolve is **quantified**: at `poly_order=4`,
   5 paired replicates give `aggressive` costing `+0.584` of recovered score in the
   Taylor basis (95% CI [+0.290, +0.878], p=0.005, 5/5) against `conservative`, for
   roughly 9-11x less pruning time and ~90-135x fewer leaves. In the Chebyshev basis the
   cost is `+2.028`. So "redundant coverage vs coverage gaps" is not hypothetical, and a
   metric-based criterion has a real ~0.6-2.0 of score to recover.

Mechanism note, relevant to Phase 2 step 2: in `shift_taylor_errors` the shift matrix
`t_mat = delta_t**powers / fact(powers)` has `powers` **zero on the diagonal**, so
`diag(t_mat)` is all ones and `"aggressive"` returns the error vector **completely
unchanged**. It does not shrink errors; it simply ignores the off-diagonal shear. That is
why it leaves gaps, and why the gaps are worse in a basis whose transform has non-unit
diagonal.

## 5. Target invariant

For any `λ*` inside a leaf's ellipsoid,
`m(λ*, λ) = δλᵀ g δλ ≤ m_max`,
where `m` is the fractional S/N loss to second order. This is a **mean-square** (time
averaged) criterion, not a worst-instant one, so it is **not equivalent** to the current
sup-norm invariant: a template can satisfy one and violate the other. Any before/after
comparison must therefore go through an explicit bridge (`m_max_from_eta`) and should be
reported as "matched along axis X" or "matched by volume", never as "same tolerance".

For a pure quadratic phase difference over an interval, the sup-norm and RMS of the same
polynomial differ by an `O(1)` factor that depends on the polynomial's shape, so the
bridge constant is not universal across `poly_order`. Recorded as an open question.

## 6. Conventions (see DECISIONS.md for the fixed list)

- **Ordering**: reverse `k`, `[d_kmax, ..., d_2, d_1]` for the branchable axes, with
  `d_0` and `f_0` in the last two rows. Any metric matrix must use the same order.
- **`Δt` sign**: `shift_taylor_params(vec, delta_t)` is called with
  `delta_t = t_new − t_old` (e.g. `t0_add − t0_cur` in `poly_taylor_resolve_batch:222`,
  `coord_next[0] − coord_cur[0]` in `poly_taylor_transform_batch:352`). `t_mat` is lower
  triangular with unit diagonal.
- **Units**: `d_k` in m·s⁻ᵏ (distance derivatives; `d_0` metres, `d_1` m/s, `d_2` m/s²),
  `f_0` in Hz, `C_VAL` in m/s (`utils/misc.py`). Phase is in **cycles**, and the
  code's own delay-to-phase conversion is `phase = f0 * d / C_VAL`
  (`poly_taylor_resolve_batch:228-234`), i.e. **no factor of 2π**.

## 7. Open questions for the human (Phase 0 exit gate)

1. **Phase 3 is blocked on PR #3.** Pruning segfaults on this base. Merge
   `fix-prune-segfault` into `metric-gridding` when Phase 3 starts, or wait for the PR?
2. **`m_max` bridge.** Match along the best-conditioned axis, or by matched ellipsoid
   volume? §5 above argues neither is canonical. Recommendation: report **both** in
   Phase 1 test 5 and pick after seeing the ratios.
3. **Harmonic weighting.** The plan asks for a `(2π)^2` factor and a harmonic weighting
   choice. The codebase works in cycles with no `2π` (see §6). Recommendation: define
   `g` so that `m` is dimensionless fractional S/N loss with the `(2π)^2` folded *into*
   `g`, and state it once in `DECISIONS.md` rather than carrying a loose factor.
4. **Does the `2**k` coarsening stay?** The metric strategy replaces the box, but
   `poly_taylor_step_f`'s coarsening also feeds `B(s)` and the base grid `G0`. Should
   `"metric"` bypass the coarsening entirely (my reading of the intent) or inherit it?
