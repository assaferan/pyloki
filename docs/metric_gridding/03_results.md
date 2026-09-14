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

### The buffer, not the covering

`max_sugg` is **not a neutral knob** for the metric strategy. Its branching is spiky —
a single parent emits **17 549** children — so at `max_sugg = 2**14 = 16 384` the buffer
cannot hold even one parent's offspring, `trim_threshold` fires mid-burst, and the true
track is discarded. Re-run at `S/N = 15` with the buffer as the only change:

| `max_sugg` | recovered | median mismatch | median score | seconds (3 reps) |
|---|---|---|---|---|
| 2^14 | 0/3 | 8.3e5 | 3.79 | 130 |
| **2^17** | **3/3** | **0.00025** | **14.06** | 1694 |

So the `0/3` in the table above is an artefact of a shared buffer, not a detection
deficit — the second time in this phase that a shared resource setting made the metric
look broken (D31 was the first, with a shared threshold ramp). **Given room to run, the
metric matches `aggressive` at `S/N = 15` (3/3 each) and localises the signal better:
median mismatch 2.5e-4 against `aggressive`'s 1.5e-3, at an equal score (14.06 vs
14.02).**

The price is the whole story:

| | `aggressive` | `metric` |
|---|---|---|
| recovered at S/N 15 | 3/3 | 3/3 |
| candidate buffer needed | 2^14 | **2^17** (8x) |
| wall-clock per injection | ~1.2 s | **~565 s** (~450x) |
| log2 complexity at equal `P_d` (step 2) | 8.16 | **18.22** (~1070x) |

The three independent cost measures agree to within a factor of a few, which is
reassuring: step 2's predicted ~1070x shows up as ~450x wall-clock and 8x memory.

## Verdict

**Inconclusive on the central question, and negative on cost.**

The metric strategy works. It covers correctly, it recovers what `aggressive` recovers,
and it localises better. It costs about three orders of magnitude more to do so.

But this config cannot answer the question the project exists to ask. The premise is
that `aggressive` loses signals in the corners its diagonal error propagation does not
track — and here `aggressive` **under-reports its own region by only 1.02x to 1.33x**
(measured, see DECISIONS "Regime dependence"). There is essentially no gap to fix, so
the metric's guarantee buys nothing and only its cost is visible. A fair test needs a
regime where `aggressive` actually drops signals, and the plan named one: the 18-minute
circular-orbit search. That config needs 3.8 TB here and `"metric"` does not support the
circular basis.

So the recommendation is **not** "abandon". It is that the decisive experiment has not
been run, and cannot be run on this machine as the code stands. The two things that
would change that, in order:

1. **A regime where `aggressive` visibly fails.** Either a larger machine for the
   plan's config, or find a Taylor-basis config on this machine where `aggressive`'s
   under-report is large. The under-report grew with baseline and order in the small
   scan (1.001x at 17 s to 1.33x at 67 s), so a longer or higher-order Taylor run is
   worth scanning for before reaching for more hardware.
2. **The circular-orbit extension** (Phase 4), without which the plan's own config is
   unreachable at any scale.

Cheaper things that would sharpen, but not settle, the picture: more noise realisations
(3 is too few — `conservative` goes 1/3, 0/3, 3/3 across S/N 12, 15, 20, which is
noise), per-strategy buffers rather than a shared one, and settling `(m_max, R)` jointly
with O6 and O7 instead of using the cost-optimal point.

## Not done

Step 4 (`P_d[q]` versus anchor segment) and step 5 (wall-clock and per-stage candidate
counts at equal recalibrated `P_d`). Step 5 is partly covered by the table above.
