# The `max_sugg` overflow ratchet is silent in the logs

**Status: draft for review. Not posted, not pushed.** One defect, one patch, one
verification. Separate from `04_upstream_report.md`, which another session owns.

## Summary

When the candidate buffer overflows, `pyloki` raises the effective pruning threshold
above the value the threshold scheme asked for. Nothing in `PruneStats` records that it
happened or what the realised cut was, so a run whose thresholds were silently tightened
is indistinguishable in its log from one that ran at the scheme's thresholds. Users
comparing configurations — which is exactly what a `threshold_scheme` is for — can be
comparing buffer pressure instead, with no indication in the output.

## Mechanism

The nominal threshold is read from the scheme and passed down (`prune.py:640`):

```python
threshold = self.threshold_scheme[self.prune_level - 1]
```

Inside `pruning_iteration_batched` a second, mutable threshold is maintained
(`prune.py:225`), and it is the one that actually cuts (`prune.py:294`):

```python
current_threshold = threshold
...
passing_mask = batch_scores >= max(threshold, current_threshold)
```

It is raised whenever the buffer overflows. `add_batch` returns the effective threshold
and it is fed back into `current_threshold` (`prune.py:325`); the raise itself is
`prune_on_overload_func` (`utils/world_tree.py:547`):

```python
effective_threshold = max(current_threshold, topk_threshold, median_threshold)
```

so once the buffer fills, the cut ratchets up to the top-K or median score and **never
comes back down** for the rest of that level.

The stats dictionary returned by the iteration carries `n_leaves`, `n_leaves_phy`,
`score_min` and `score_max` (`prune.py:342-346`) — `current_threshold` is not among
them. The caller then records the nominal value (`prune.py:682-685`):

```python
pstats_cur = PruneStats(
    level=self.prune_level,
    seg_idx=seg_idx_cur,
    threshold=threshold,          # <- the scheme value, not the realised cut
```

and that is what `get_summary` prints as `score thresh:` (`io/cands.py:98`).

**The realised cut cannot be recovered from the log.** `score_min` does not help: it is
taken over every scored leaf, before thresholding (`prune.py:288`), not over the
survivors. `n_leaves_surv` equal to `max_sugg` does reveal that the buffer was full at
the end of a level, but not the value of the cut, and it misses the case where the
ratchet fired mid-level and survivors later fell below the buffer size.

## Why it is worth fixing

Measured on this branch, on **one fixed noise realisation**, changing only `max_sugg`
from 2^14 to 2^18 changed the surviving candidate count from 0 to 110 (`aggressive`) and
from 8 553 to 134 418 (`quadrature`). Over a fixed set of 20 realisations `aggressive`
recovered 11/20 at 2^14 and 13/20 at 2^18. Same data, same scheme, same everything else.

A direct probe showed the *final*-stage cut was the nominal one in all four runs tested
(minimum surviving score below `thresholds[-1]`), so the ratchet is firing at
**intermediate** stages — where it changes what survives to be scored later, and where
there is nothing at all in the output to show it.

This is a reporting defect, not a correctness defect. Ratcheting on overflow is a
reasonable way to stay inside a fixed buffer; the argument is only that it should be
visible.

### For a high-branching configuration the ratchet is not occasional — it is permanent

The case above would be weak if raising `max_sugg` reliably escaped the ratchet: a user
who noticed could just raise it. Measured on one fixed set of 10 realisations, sweeping
the buffer at a 64-segment, `poly_order = 4`, `branch_max = 16` Chebyshev search:

| branching-pattern product | `max_sugg` | median `ncand` | median saturation |
|---|---|---|---|
| 1.5 × 10¹² (low) | 2^18 | 16 778 | 0.064 |
| | 2^19 | **16 778** | 0.032 |
| | 2^20 | **16 778** | 0.016 |
| 5.4 × 10²⁰ (high) | 2^18 | 152 007 | 0.580 |
| | 2^19 | 390 575 | 0.745 |
| | 2^20 | 808 116 | **0.771** |

The low-branching configuration converges: identical candidate counts at three buffers,
saturation falling, the ratchet gone. The high-branching one **never** converges —
`d log2(ncand) / d log2(max_sugg)` is 0.95, 1.36 and 1.05 across three successive
increases, and saturation *rises*. Doubling the buffer doubles the count, because the
ratchet relaxes the cut to keep the buffer full at whatever size it is given.

**So a user running a high-branching search is running against the ratchet at every
buffer they can afford, permanently, and the log tells them their threshold scheme is
being applied when it is not.** There is no value of `max_sugg` at which they would find
out, and no field in the output that would tell them. That is the case for this patch,
independent of anything about tiling strategies: it is what makes the difference between
an occasional tightening a careful user could detect and a silent, permanent
substitution of a different cut for the one they configured.

## Patch

Three lines. `current_threshold` is already a float, so adding it to the stats dict does
not reintroduce the heterogeneous-dict problem the comment at `prune.py:336-341`
describes.

```diff
--- a/src/pyloki/prune.py
+++ b/src/pyloki/prune.py
@@ stats = {
         "n_leaves": float(n_leaves),
         "n_leaves_phy": float(n_leaves_phy),
         "score_min": score_min if np.isfinite(score_min) else 0.0,
         "score_max": score_max if np.isfinite(score_max) else 0.0,
+        "threshold_eff": float(current_threshold),
     }
@@ pstats_cur = PruneStats(
         threshold=threshold,
+        threshold_eff=stats_dict["threshold_eff"],
```

with a matching `threshold_eff: float` attribute on `PruneStats` (`io/cands.py:17`) and
one more field in `get_summary`, printed next to the nominal one so a tightened run is
obvious at a glance.

## What this note does not claim

It does **not** claim the ratchet is arm-dependent — that it penalises a
higher-branching `tiling_strategy` more than a lower-branching one. That is a separate
question, it is being measured separately, and it is not needed for this fix: the case
for recording the realised threshold stands on the fact that it differs from the logged
one at all.

## Verification

**Done.** The patch is implemented and verified on branch `upstream-max-sugg-logging`,
off `upstream/main` at `18d04b3` — a separate worktree, deliberately not this one, so
`src/` here was never touched while a measurement was running.

Measured with the shipped scaled-down Extreme-Pruning pipeline
(`tests/test_example_ep_accel.py`'s configuration) on **one fixed time series**, run at
two buffers:

| `max_sugg` | levels logged | levels with `eff > thresh` |
|---|---|---|
| 2^20 (non-binding) | 8 | **0** — every level logged `eff == thresh` |
| 2^10 (binding) | 8 | **4** |

and in the binding case the tightening is large:

    thresh  3.00  ->  eff  3.96   (+0.96)
    thresh  3.63  ->  eff  5.41   (+1.78)
    thresh  3.42  ->  eff  6.24   (+2.82)
    thresh  3.90  ->  eff  7.34   (+3.44)

**Nearly a factor of two on the realised cut at the last stage, and before this patch
none of it appeared anywhere in the output.**

The final implementation is slightly larger than the three lines sketched above: the
stats dict gains `threshold_eff`, `PruneStats` gains the field (defaulting to `NaN` so
existing direct constructions keep working) and prints it beside the nominal value, and
the initial pre-pruning record reports 0 rather than `NaN`.

*Checks run:* `src/` is ruff-clean; the upstream suite passes (79 tests, everything bar
the long `test_example_ep_circular.py`); and
`tests/test_prune_threshold_eff.py` asserts **both** directions — the log must show the
ratchet when the buffer binds and must not when it does not. The second half matters:
a patch that always reported a raised threshold would pass a one-sided test and be
useless.

## Provenance of the motivating numbers

The saturation figures quoted above under "Why it is worth fixing" come from this
branch's own campaign work (`05_injection_design.md` §6.5) and are reproduced there with
all 200 runs. They are what prompted the patch; they are not needed to justify it. The
verification above stands on the shipped example alone and depends on nothing from the
`metric-gridding` or `injection-design` branches.
