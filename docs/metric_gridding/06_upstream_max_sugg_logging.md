# The `max_sugg` overflow ratchet is silent in the logs

**Status: draft for review. Not posted, not pushed.** One defect, one patch, one
verification. Separate from `04_upstream_report.md`, which another session owns.

*All line numbers refer to `upstream/main` at `18d04b3`.*

## Summary

When the candidate buffer overflows, `pyloki` raises the effective pruning threshold
above the value the threshold scheme asked for. Nothing in `PruneStats` records that it
happened or what the realised cut was, so a run whose thresholds were silently tightened
is indistinguishable in its log from one that ran at the scheme's thresholds. Users
comparing configurations — which is exactly what a `threshold_scheme` is for — can be
comparing buffer pressure instead, with no indication in the output.

## Mechanism

The nominal threshold is read from the scheme and passed down (`prune.py:626`):

```python
threshold = self.threshold_scheme[self.prune_level - 1]
```

Inside `pruning_iteration_batched` a second, mutable threshold is maintained
(`prune.py:218`), and it is the one that actually cuts (`prune.py:285`):

```python
current_threshold = threshold
...
passing_mask = batch_scores >= max(threshold, current_threshold)
```

It is raised whenever the buffer overflows. `add_batch` returns the effective threshold
and it is fed back into `current_threshold` (`prune.py:315`); the raise itself is
`prune_on_overload_func` (`utils/world_tree.py:547`):

```python
effective_threshold = max(current_threshold, topk_threshold, median_threshold)
```

so once the buffer fills, the cut ratchets up to the top-K or median score and **never
comes back down** for the rest of that level.

The stats dictionary returned by the iteration carries `n_leaves`, `n_leaves_phy`,
`score_min` and `score_max` (`prune.py:331-335`) — `current_threshold` is not among
them. The caller then records the nominal value (`prune.py:653-656`):

```python
pstats_cur = PruneStats(
    level=self.prune_level,
    seg_idx=seg_idx_cur,
    threshold=threshold,          # <- the scheme value, not the realised cut
```

and that is what `get_summary` prints as `score thresh:` (`io/cands.py:98`).

**The realised cut cannot be recovered from the log.** `score_min` does not help: it is
taken over every scored leaf, before thresholding (`prune.py:279`), not over the
survivors. `n_leaves_surv` equal to `max_sugg` does reveal that the buffer was full at
the end of a level, but not the value of the cut, and it misses the case where the
ratchet fired mid-level and survivors later fell below the buffer size.

## Why it is worth fixing

Measured on this branch, on **one fixed noise realisation**, changing only `max_sugg`
from 2^14 to 2^18 changed the surviving candidate count from 0 to 110 (`aggressive`) and
from 8 553 to 134 418 (`quadrature`). Over a fixed set of **50** realisations, raising
the buffer from 2^14 to 2^18 moved `aggressive` from 30 to 38 recoveries and
`quadrature` from 28 to 40 — **20 flips, every one toward recovery and none away**
(McNemar `p` = 0.008 and `p` < 0.001). Same data, same scheme, same everything else.

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
the buffer at a single 64-segment, `poly_order = 4`, `branch_max = 16` Chebyshev search — the two rows are the two `tiling_strategy` values of **one** configuration, not two configurations:

| tiling strategy (branching product) | `max_sugg` | median `ncand` | median saturation |
|---|---|---|---|
| `aggressive` (1.5 × 10¹²) | 2^18 | 16 778 | 0.064 |
| | 2^19 | **16 778** | 0.032 |
| | 2^20 | **16 778** | 0.016 |
| | 2^21 | **16 778** | **0.008** |
| `quadrature` (5.4 × 10²⁰) | 2^18 | 152 007 | 0.580 |
| | 2^19 | 390 575 | 0.745 |
| | 2^20 | 808 116 | 0.771 |
| | 2^21 | **1 610 052** | **0.768** |

The low-branching strategy converges: **identical candidate counts at four
buffers spanning a factor of 8**, saturation falling to 0.008, the ratchet gone. The
high-branching one **never** converges — over the same factor of 8 its count grows by a
factor of 10.6, with `d log2(ncand) / d log2(max_sugg)` = 1.135 over the full span and
1.022 and 0.994 over the two widest sub-spans, while saturation stays pinned near 0.77.
Doubling the buffer doubles the count, because the ratchet relaxes the cut to keep the
buffer full at whatever size it is given. (Span exponents are quoted rather than
step-to-step ones, which are noisy enough to invite reading a knee into noise — a single
realisation gave 1.54 then 0.49 on consecutive steps.)

**So a user running a high-branching search is running against the ratchet at every
buffer they can afford, permanently, and the log tells them their threshold scheme is
being applied when it is not.** There is no value of `max_sugg` at which they would find
out, and no field in the output that would tell them. That is the case for this patch,
independent of anything about tiling strategies: it is what makes the difference between
an occasional tightening a careful user could detect and a silent, permanent
substitution of a different cut for the one they configured.

## Patch

Implemented on branch `upstream-max-sugg-logging` off `upstream/main` at `18d04b3`.
`current_threshold` is already a float, so adding it to the stats dict does not
reintroduce the heterogeneous-dict problem the comment at `prune.py:326-330` describes.
The `NaN` default on the new field keeps existing direct constructions of `PruneStats`
working. **+16 / −2 lines across two files**, plus a regression test.

This diff is the committed one, copied from `git diff`, not transcribed:

```diff
diff --git a/src/pyloki/io/cands.py b/src/pyloki/io/cands.py
index ec979fb..84ae89e 100644
--- a/src/pyloki/io/cands.py
+++ b/src/pyloki/io/cands.py
@@ -24,7 +24,11 @@ class PruneStats:
     seg_idx : int
         The segment index being added.
     threshold : float
-        The threshold value.
+        The nominal threshold from the threshold scheme.
+    threshold_eff : float
+        The threshold actually applied. Equal to `threshold` unless the candidate
+        buffer overflowed, in which case pruning ratchets the cut up to
+        `max(threshold, top-K, median)` and never lowers it again for that level.
     score_min : float
         The minimum leaf score.
     score_max : float
@@ -50,6 +54,7 @@ class PruneStats:
     level: int
     seg_idx: int
     threshold: float
+    threshold_eff: float = float("nan")
     score_min: float = 0.0
     score_max: float = 0.0
     n_branches: int = 1
@@ -95,7 +100,8 @@ class PruneStats:
             f"branch_frac: {self.branch_frac:5.2f},",
         )
         summary.append(
-            f"score thresh: {self.threshold:5.2f}, max: {self.score_max:5.2f}, "
+            f"score thresh: {self.threshold:5.2f}, eff: {self.threshold_eff:5.2f}, "
+            f"max: {self.score_max:5.2f}, "
             f"min: {self.score_min:5.2f}, P(surv): {self.surv_frac:4.2f}",
         )
         return "".join(summary) + "\n"
diff --git a/src/pyloki/prune.py b/src/pyloki/prune.py
index 22499de..34209c6 100644
--- a/src/pyloki/prune.py
+++ b/src/pyloki/prune.py
@@ -333,6 +333,12 @@ def pruning_iteration_batched(
         "n_leaves_phy": float(n_leaves_phy),
         "score_min": score_min if np.isfinite(score_min) else 0.0,
         "score_max": score_max if np.isfinite(score_max) else 0.0,
+        # The cut actually applied. Differs from the caller's nominal threshold
+        # whenever the buffer overflowed (utils/world_tree.py: prune_on_overload_func
+        # ratchets it to max(threshold, top-K, median)). Without this the tightening
+        # is invisible: PruneStats logs the scheme value, and `score_min` is taken
+        # over every scored leaf before thresholding, not over the survivors.
+        "threshold_eff": float(current_threshold),
     }
     return tree_new, stats, timers
 
@@ -517,6 +523,7 @@ class Pruning:
             level=self.prune_level,
             seg_idx=self.scheme.get_segment_idx(self.prune_level),
             threshold=0,
+            threshold_eff=0,   # initial record: nothing has been pruned yet
             score_min=self.world_tree.score_min,
             score_max=self.world_tree.score_max,
             n_branches=self.world_tree.size,
@@ -654,6 +661,7 @@ class Pruning:
             level=self.prune_level,
             seg_idx=seg_idx_cur,
             threshold=threshold,
+            threshold_eff=stats_dict["threshold_eff"],
             n_branches=self.world_tree.valid_size,
             n_leaves_surv=world_tree.valid_size,
             # stats_dict values are all float (see pruning_iteration_batched);
```

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
