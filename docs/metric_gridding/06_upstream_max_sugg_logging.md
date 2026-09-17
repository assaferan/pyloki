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

To be filled in once the patch is applied and a saturating run is logged. Deliberately
not done while the paired experiment is in flight: editing `prune.py` would change the
library under a running measurement.
