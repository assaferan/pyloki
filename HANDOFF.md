# Handoff: pyloki fork — Chebyshev-economized resolve, its parity limit, and the prune segfault

Updated 2026-09-09 (second session). Working note, not part of the codebase.

## 1. Repo / remotes

- `origin` = `git@github.com:assaferan/pyloki.git` (fork, pushable)
- `upstream` = `git@github.com:pravirkr/pyloki.git` (read-only)
- `main` tracks `origin/main`. Local venv `.venv/` must be built from a
  tkinter-capable Python (`/opt/homebrew/bin/python3.13`); see §1 of the git history
  for the original note. 35 tests pass as of this writing (32 pre-existing + 2 parity + 1 prune smoke).

## 2. UNCOMMITTED WORK — read this first

`git status` currently shows these modified, **none committed**:

| file | change |
|---|---|
| `src/pyloki/prune.py` | **Segfault fix** (§4). Independent of the economization work and by far the most important change here. |
| `src/pyloki/utils/transforms.py` | `Notes` section on `economize_taylor_params` documenting the parity limitation (§3). |
| `src/pyloki/core/taylor.py` | Comments at the three resolve call sites pointing at that limitation; two 89-char lines wrapped to satisfy the configured 88-char limit. |
| `tests/test_transforms.py` | Two regression tests for the parity limitation. |
| `tests/test_prune.py` | **New.** First test to touch the pruning path at all (§4). |
| `paper/QUESTIONS.md` | Pre-existing §3.2 answer from the previous session, still uncommitted. |
| `.gitignore` | Pre-existing change from the previous session. |

The `prune.py` fix is worth committing and pushing on its own merits, separately from
the economization thread.

## 3. The parity limitation (settled, now documented in-tree)

Confirmed and generally proven (`scratch/econ_experiment/parity_check.py`, plus the
two new tests): economization perturbs only the retained coefficients that **share the
parity** of a dropped order, because a Chebyshev polynomial of even/odd order expands
into even/odd monomials only. When exactly one order is dropped, the result is
bit-identical to naive truncation for every opposite-parity coefficient.

For the `resolve` step (`n_keep=3` always):

- **`poly_order=3`** (jerk search, the common case): only jerk (order 3, odd) is
  dropped, so **accel (order 2) and delay (order 0) are untouched** — the resolved
  accel cell is exactly what naive truncation picks. Only velocity/frequency moves.
- **`poly_order>=4`**: two-plus orders dropped, both parities spanned, accel does move.

**The sting:** the truncation error is dominated by the *jerk* term (lowest dropped
order, entering at `t_s**3` against snap's `t_s**4`). A `size_experiment.py` bisection
found that at `tseg = 100 s` the injected jerk of 6 alone exceeds one phase bin with
snap set to zero. So the configuration where truncation error matters most is precisely
the one where the fix does least for the accel grid.

## 4. The prune segfault: ROOT-CAUSED AND FIXED

The previous session recorded an unexplained `Segmentation fault: 11` and guessed it was
an edge case of extreme small scale. **That guess was wrong.** The crash was
configuration-independent — it reproduced at every `nsegments`, `max_sugg`,
`branch_max`, `batch_size` (including `batch_size=1`), and at both `poly_order=3` and
`poly_order=4`. The previous session's own `run_search.py` also died immediately.

**Cause.** `pruning_iteration_batched` (`src/pyloki/prune.py`) returned

```python
stats = {"n_leaves": n_leaves,          # int64
         "n_leaves_phy": n_leaves_phy,  # int64
         "score_min": ...,              # float64
         "score_max": ...}              # float64
```

Numba cannot box a heterogeneous `LiteralStrKey` dict back to Python. It fails with
`TypeError: cannot convert native LiteralStrKey[Dict](...) to Python object`, returns
NULL with an exception set, and the caller's three-way tuple unpack at
`prune.py:634` then dereferences NULL. The macOS crash report confirms it:
`EXC_BAD_ACCESS (SIGSEGV) KERN_INVALID_ADDRESS at 0x0` with the faulting frame in
CPython's `_PyEval_UnpackIterable` and **no numba frames on the stack**.

Minimal reproducer (three lines, no pyloki): an njit function returning
`{"a": 1, "b": 2.5}` fails; `{"a": 1.0, "b": 2.5}` and `{"a": 1, "b": 2}` both box fine.

**Fix.** Cast the two counts to `float` so the dict is homogeneous, and restore `int`
at the single Python consumer. Note the now-commented-out older implementation at
`prune.py:156` already wrote `float(n_leaves)` — the batched rewrite dropped that cast.
This is a regression, not a novel constraint.

Environment where it bites: numba 0.67.0, numpy 2.5.3, Python 3.13.15.

**Why CI missed it:** the suite had *zero* coverage of `prune_dyp_tree`. `tests/test_prune.py`
now runs a minimal prune end to end in ~11 s. A segfault can't be caught in-process, so
the test guards the path by crashing the worker if the boxing regresses.

## 5. Corrected: the sup-norm probe bug was a SHAPE bug, not units

The previous handoff called this "an unresolved units/scaling bug." Wrong, and
misleading to whoever picks it up. `phase_of_exact_model` was called with a `(1,1)`
array as `delta_t`; `shift_taylor_params` broadcasts that silently, so `phi_exact` came
out `(401,1)` instead of `(401,)` and `phi_naive - phi_exact` broadcast to a
`(401,401)` **outer difference**. `np.max(np.abs(...))` then returned the full phase
*spread* across the segment — `freq*tseg/2` after midpoint alignment — which is exactly
the nonsensical `149.796` reported for both methods.

Fixed in `resolve_probe2.py` (old version kept as `.py.bak`) by evaluating the exact
polynomial directly, as `tests/test_transforms.py` already did.

**With shapes fixed, the metric works and the ratio is clean: economization reduces
sup-norm phase error to 0.2438 of naive truncation (spread 0.2429–0.2447) across six
octaves of segment duration.** Rows where the error exceeds 1 cycle are excluded — there
the constant-accel model has broken down entirely and the ratio drifts meaninglessly
above 1.

## 6. Experiment infrastructure in `scratch/econ_experiment/` (gitignored)

Rebuilt this session. The previous session's `data.npz` + `run_search.py` **cannot
answer the question** — that data has no snap at all and `tobs = 2.1 s`, giving a
resolve error of ~1e-7 cycles, and `run_search.py` uses `poly_order=3`, the known null
case. Both are kept only for reference.

- `econ_metrics.py` — shared, importable metrics: `resolve()` (mirrors
  `poly_taylor_resolve_batch`), `sup_norm_phase_error()`, `phase_error_pair()`,
  `snap_for_phase_budget()` (bisection), `max_v_over_c()`.
- `size_experiment.py` — **screens a candidate config against five preconditions
  without paying for pruning** (~6 s). Use this before any expensive run; it is what
  made the working config findable instead of guessed.
- `make_data_snap.py` — generates `data_snap.npz` at the validated config, solving snap
  so the naive error is exactly one phase bin.
- `build_variant_naive.py` — rebuilds `variant_naive/` from the **current** `src/pyloki`
  with the three economize calls stripped. **Re-run after any `src/` change**: the first
  naive run silently used a pre-fix copy and segfaulted while econ succeeded.
- `run_snr.py <econ|naive>` — one end-to-end run. Asserts which pyloki tree it loaded,
  so a PYTHONPATH slip cannot fake a null result.
- `replicate.py N` / `analyze_replicates.py` — paired replicates and their analysis.
- `resolve_probe2.py` — fixed sup-norm probe plus the scale sweep (§5).
- `parity_check.py` — fast, self-contained proof of the parity limitation (§3).

### The five preconditions, and why they fight each other

1. `poly_order >= 4` — else accel is provably untouched (§3).
2. accel grid must have > 1 cell, or the resolved cell cannot differ.
3. naive resolve error must be a real fraction of a phase bin at
   `half_width_add = 0.5 * tseg_ffa` (`snail.py:249`, fed from `dyp.tseg`).
4. `nsegments` moderate.
5. **Neglected 2nd-order Doppler must stay well below the 1st-order error under test.**
   pyloki models `f = f0*(1 - v/c)` to first order only. Raising snap to satisfy (3)
   also raises `v/c`, and the ratio of neglected-to-measured error goes as
   `nseg**7 / (freq * tseg)` — brutal in `nseg`. The first config that satisfied (2)+(3)
   drove `v/c` to **0.19**, with a 2nd-order error of 172 cycles swamping the 1-bin
   effect; its SNR numbers were meaningless. This is the trap to avoid.

(2) and (5) are coupled through snap, because `ParamLimits.from_upper` derives the accel
range from the snap bracket. `--accel-widen` in `size_experiment.py` decouples them by
scaling the accel search range about its centre.

Also note `_check_init_param_arr` (`ffa.py:429`) requires every non-frequency parameter
to occupy exactly **one** initial FFA cell, so the snap bracket must be narrow and
centred on the true value, not spanning zero.

### The validated config

`nsamps=2**20`, `dt=64e-6` (`tobs=67.1 s`), `period=0.007`, `nseg=4`
(`tseg_ffa=16.78 s`), `brute_div=128`, `accel_widen=20`, `nbins=64`,
`poly_order=4`, `accel=500`, `jerk=6`, **`snap=17.336`** (solved for one phase bin).
Gives `grid=[1, 1, 15, 1143]`, 4 pruning stages, `v/c=4.3e-4`,
2nd-order/naive = 0.114. All five preconditions pass. One pair of runs ≈ 35 s.

Clarification the previous handoff got half right: `nsegments` is `nsamps//bseg_brute`
**pre**-execute, but the pruning loop iterates the **post**-execute value,
`tobs/tseg_ffa`. Size the experiment with the latter.

## 7. End-to-end SNR at `poly_order=4`: NO measurable difference

With the segfault fixed and a physically valid config (§6), the pipeline runs to
completion for both variants. 8 **paired** replicates (each regenerates the data —
`PulseSignalConfig` uses an unseeded `default_rng` — and runs both variants on that same
realization, so realization variance cancels):

| metric | naive | econ | mean delta | sd(delta) | p (paired t) | econ wins |
|---|---|---|---|---|---|---|
| closest-to-true score | 7.1673 | 7.1358 | **-0.0315** | 0.3281 | **0.79** | 4/8 |
| closest score_ep | 7.1731 | 7.2783 | +0.1052 | 0.2739 | 0.31 | 4/8 |
| best score | 11.3475 | 11.3999 | +0.0524 | 0.2641 | 0.59 | 4/8 |
| best score_ep | 11.2668 | 10.8755 | -0.3913 | 0.5249 | 0.073 | 1/8 |
| n candidates | 4708.6 | 4677.4 | -31.3 | 144.9 | 0.56 | 3/8 |

Per-replicate deltas in closest-to-true score:
`+0.078 +0.153 -0.165 -0.250 +0.507 -0.104 -0.598 +0.128` — signs mixed, no trend.

**95% CI on the mean delta: [-0.306, +0.243], i.e. [-4.3%, +3.4%] of score.**

This is not merely "underpowered". A crude smearing estimate (a `ducy=0.1` pulse is
~6.4 bins wide at `nbins=64`; a 1-bin sup-norm error costs ~7% amplitude, a 0.244-bin
error ~1.9%) predicts econ should win by **+0.37 score**. The measurement is
`-0.03 +/- 0.12`, which **excludes that prediction at >3 sigma**. So the effect is not
hiding below the noise — the naive model of how resolve phase error turns into lost
sensitivity is wrong.

**Leading hypothesis for the null.** `poly_taylor_resolve_batch` does not only pick a
grid cell; it also returns `relative_phase`, computed from the (economized) `delay`, and
the pipeline applies that shift explicitly. Much of the phase error the sup-norm metric
charges against naive truncation may therefore be absorbed by the tree's own phase
bookkeeping, never reaching the folded profile. If so, the sup-norm-over-the-segment
metric — including the corrected one in `resolve_probe2.py` and the assertion in
`tests/test_transforms.py` — systematically **overstates** the sensitivity relevance of
economization. Checking that is the single most valuable next step: instrument
`relative_phase` and the recovered profile directly, rather than the coefficients.

**Bottom line.** Economization provably reduces sup-norm phase error to 0.244x naive
(§5), is a no-op for accel at `poly_order=3` (§3), and produces **no detectable
end-to-end SNR gain at `poly_order=4`** in the one regime where the experiment is
physically valid. It is defensible as a correctness/robustness improvement, not as a
sensitivity improvement. Nothing in the commit message for `ef20ac0` or the plan file
says this.

## 8. Suggested next steps

1. **Commit and push the `prune.py` segfault fix + `tests/test_prune.py`** (§4),
   separately from the economization thread. This is a real bug that made every
   `prune_dyp_tree` call die; it is worth an upstream issue/PR to `pravirkr/pyloki`.
2. Test the §7 hypothesis: does `relative_phase` absorb the truncation error? If yes,
   soften the sensitivity claim in `test_transforms.py`'s docstring and reconsider
   whether economization is worth its cost at all.
3. If a sharper end-to-end bound is still wanted, note that `nseg` is capped hard by
   precondition 5 (`nseg**7`). Raising `freq` (shorter period, finer `dt`) buys
   headroom; going past `nseg=4` at `freq=142.86 Hz` requires `tobs` in the thousands
   of seconds.
4. Decide whether the `poly_order=3` null (§3) plus the `poly_order=4` null (§7)
   together warrant reverting `ef20ac0`, or keeping it with the documented caveats now
   in `economize_taylor_params`.
5. Still open from the previous session: ~10 unanswered questions in `paper/QUESTIONS.md`,
   and its uncommitted §3.2 answer.
