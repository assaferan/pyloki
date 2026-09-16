# 04_upstream_report.md — draft report for pravirkr/pyloki

Status: **draft, not yet posted.** Upstream has Discussions disabled, so this would be an
issue. Numbers re-verified 2026-09-16; see DECISIONS.md session (p).

**Title:** `tiling_strategy` cannot affect the final grid cell, so `quadrature`/`conservative` buy no sensitivity for their cost

Section 5.2.4 of the paper closes by asking for the sensitivity loss of aggressive tiling to be quantified. I measured it, and the answer is more specific than expected: the choice of `tiling_strategy` has no effect on sensitivity at all.

Sampling parameter offsets uniformly inside each leaf's own cell and evaluating the incurred phase error against the promised `eta / N_b` (268.4 s, 64 segments, `poly_order=4`, `N_b=64`, `eta=1.0`, 7 ms spin period), in units of that tolerance:

| `tiling_strategy` | median over stages | worst corner | `prod B(s)` |
|---|---|---|---|
| `aggressive` | 3.35 | 16.80 | 1.5e12 |
| `quadrature` | 3.58 | 14.38 | 1.1e17 |
| `conservative` | 3.42 | 14.20 | 1.5e32 |

Twenty orders of magnitude of cost, and the phase error is the same to within the sampling noise.

### Why

`branch_param_padded` (`utils/psr_utils.py:360-368`) sets

```python
num_points = max(1, math.ceil(dparam_cur / dparam_new - FLOAT_EPSILON))
dparam_new_actual = dparam_cur / num_points
```

so a branched child cell is `dparam_cur / ceil(dparam_cur / dparam_new)`, which is always `<= dparam_new` and depends on `dparam_cur` only through the integer `ceil`. The complementary case is in `core/taylor.py:148-154`: when the accumulated shift on an axis is below `eta`, the axis is not branched at all and keeps its `dparam_cur`.

Either way the criterion, not the transport, fixes the scale: a branched axis is pulled to just under `dparam_new`, and an unbranched one is by construction still within tolerance. The transported width survives only as the `ceil` remainder and as the slack in that guard, which is why the measured cells sit within a factor of two of the criterion on both sides.

`shift_taylor_errors` is what the three strategies differ in, and it sets only `dparam_cur`. So the strategy controls how much redundant subdivision happens on the way to the cell, never the cell itself. Measured final half-widths against the criterion step `[9.7e-3, 0.163, 3.64, 122]`:

| `aggressive` | `[1.10e-2, 0.123, 5.00, 96.4]` |
| `quadrature` | `[1.10e-2, 0.288, 4.02, 133]` |
| `conservative` | `[1.10e-2, 0.344, 6.99, 239]` |

All within a factor of two of the criterion, on both sides — cells can exceed it where the shift guard skipped branching. This also means the geometric inflation that Figure `moving_reference` costs out is pure overhead in the current branching implementation — it cannot translate into coverage, because the cell it produces is re-refined to the criterion at the next stage regardless.

### Separately: what `eta` delivers

The same measurement pins what `eta` means for the shipped default, and it is not what eq. `grid_criteria` suggests. Two self-checks:

- the naive per-axis grid (`use_cheby=False`) costs exactly one tolerance at a full step — the promise is kept per axis, to 1e-6;
- the shipped grid is coarser by exactly `2**(k-1)` on the order-`k` axis, i.e. eq. `dk_optimal` / appendix `app:optimal_gridding`, as intended.

The coarsening is per-coefficient and derived in the Chebyshev basis, where `|T_k| <= 1` bounds each coefficient independently. Applied as independent bounds on monomial `d_k` cells, a signal offset in several coefficients at once accumulates: hence `3.4x` for a typical signal and `~15x` at a corner. So a sensitivity estimate taken from `eta` alone is optimistic by about a factor of three under the recommended `aggressive` Taylor default. Tightening `eta` to recover the nominal budget is the remedy 5.2.4 already suggests; this just puts a number on how much.

### Caveats

Deterministic geometry only. Pruning is not modelled, and it is the one place a wider claimed region could still help — a signal near a cell edge may sit in a leaf that was thresholded away, and a strategy claiming more territory keeps more such leaves alive. That is a threshold-scheme interaction and needs injections to settle; nothing here rules it out. Taylor basis only: the Chebyshev and circular transforms are not unit-diagonal triangular, so the argument above does not carry over to them.

Reproducer: `docs/metric_gridding/sensitivity_loss.py` on
https://github.com/assaferan/pyloki/tree/metric-gridding — self-contained, pure NumPy, carries the two self-checks above. Happy to open a PR adding it as a test if useful.
