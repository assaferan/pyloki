# Phase 3 results — Taylor analogue

Deliverable for `metric_PLAN.md` Phase 3. **Provisional**: steps 1, 2 and 3 are done,
steps 4 and 5 are not, and step 3 is at low statistics.

## Configuration

Not the plan's config, for two measured reasons (D29, D30):

- the plan specifies a **circular-orbit** search, which is `prune_poly_order = 5`, and
  `"metric"` is refused on the circular basis (D16) because `core/metric.py` builds `g`
  on the Taylor kinematic basis;
- the plan's **18-minute / 128-segment** config needs a **3.8 TB** FFA fold array and
  cannot run on this machine.

So: Taylor basis, `poly_order = 4`, **268.4 s over 64 segments**, `eta = 1`, `N_b = 64`,
`ducy_max = 0.2`, spin period 7 ms. That is the smallest config that both fits in memory
and has a non-degenerate base grid — at 67 s the grid is `[1, 1, 1, 12]`, only frequency
is searched, and all three strategies return an identical answer. Pinned in
`phase3_config.py`. The metric runs at `(m_max, R) = (0.00208, 16)`, its cost-optimal
point for matching `conservative`'s worst-case mismatch.

The paper's Figures 8/11/12 are circular-orbit panels, so what follows tests the same
*mechanism* — corner gaps versus metric coverage — and is not a reproduction of them.

## Step 1 — branching factor

| strategy | `prod B(s)` | levels that branch | max `B` |
|---|---|---|---|
| `aggressive` | 1.51e12 | 19 / 63 | 27 |
| `metric` | 1.35e19 | 5 / 63 | 17 549 |
| `conservative` | 5.29e32 | 56 / 63 | 96 |

The two box strategies bracket the tiling dilemma and differ by twenty orders of
magnitude. The metric sits between them, and branches in a handful of enormous bursts
rather than steadily.

## Step 2 — recalibrated thresholds

Viterbi-optimised per strategy at `P_d = 0.1`. Required fixing two upstream bugs
(issues #8 and #9, PR #10) that made the optimiser return no scheme at all.

| strategy | `P_d` achieved | log2 complexity | log2 cost | thresholds |
|---|---|---|---|---|
| `aggressive` | 0.1031 | **8.16** | 11.44 | 2.30 – 7.70 |
| `metric` | 0.1031 | **18.22** | 21.50 | 1.40 – 8.40 |
| `conservative` | 0.1031 | **49.99** | 53.27 | 1.70 – 9.10 |

At equal detection probability the metric costs `2^10.1` = **1070x** `aggressive`, and
`2^31.8` = **3.7e9 x less** than `conservative`.

## Step 3 — injection-recovery

Injected `accel = 1.0`, `jerk = 0.05`, `snap = 0.001` — a signal the base grid cannot
resolve on its own, so pruning has to find all four axes. Recovery is judged by mismatch
to the truth in the full-baseline metric, not parameter distance, so the criterion means
the same thing for every strategy. Three noise realisations per point,
`max_sugg = 2**14` for all.

| injected S/N | `aggressive` | `metric` | `conservative` |
|---|---|---|---|
| 6 | 0/3 | 0/3 | 0/3 |
| 8 | 0/3 | 0/3 | 0/3 |
| 10 | 0/3 | 0/3 | 0/3 |
| 12 | 1/3 | 0/3 | 1/3 |
| 15 | **3/3** | **0/3** | 0/3 |
| 20 | 3/3 | 2/3 | 3/3 |

Wall-clock per injection: `aggressive` 2-4 s, `conservative` ~11 s, `metric` 100-260 s.

