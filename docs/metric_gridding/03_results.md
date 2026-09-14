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

## The benefit side — and it appears to be zero

Every measurement above is a **cost**. This is the missing half: how much does
`aggressive` actually lose to the coverage gaps the project exists to close?

**Measured answer: nothing, geometrically.**

`shift_taylor_errors` builds `T` lower-triangular with **unit diagonal** (C4). So row `i`
of `T @ ((idx - n) * spacing)` depends on `n_i` only through `-n_i * spacing_i`, and
choosing `n_0, n_1, ...` in order brings every row's residual within `spacing_i / 2` by
rounding. An axis-aligned box of half-width `spacing / 2`, placed at the **sheared**
lattice of leaf centres, therefore tiles the space exactly — for any shear.

And `aggressive` records exactly `spacing * |diag(T)| == spacing`. Its widths are
unchanged because they *should* be: the lattice shears underneath them but each axis
keeps its period.

Uncovered fraction of a sheared cell, `poly_order = 4`, real Phase 3 spacings:

| `delta_t` | `aggressive` | `quadrature` | `conservative` |
|---|---|---|---|
| 0.5 | 0.0000 | 0.0000 | 0.0000 |
| 2.0 | 0.0000 | 0.0000 | 0.0000 |
| 10.0 | 0.0000 | 0.0000 | 0.0000 |
| 100.0 | 0.0000 | 0.0000 | 0.0000 |

Verified two independent ways: forward substitution, and brute force over an explicit
±400 lattice in 2D. They agree to the fourth decimal.

**This took three attempts and the first two were wrong**, both in the same direction —
searching a fixed block of neighbouring lattice points instead of solving for the
covering one. Because `T` is triangular, the covering index can be arbitrarily far from
the nearest one on the lower-order axes, so a too-small window reports spurious holes
(the second attempt claimed `aggressive` leaves 100% uncovered). Recorded because the
error is easy to repeat.

### What this means

If it holds, the premise of the project does not: `aggressive` is not gappy, it is
exactly sufficient, and `quadrature` and `conservative` record wider boxes for **no
coverage benefit at all** — `quadrature` costs 1.10e17 against `aggressive`'s 1.51e12,
and `conservative` 1.48e32, all covering identically.

The metric strategy is then paying ~1070x at equal `P_d` to guarantee something
`aggressive` already guarantees for free. That is consistent with step 3, where
`aggressive` matched or beat it.

### Caveats, in order of how much they could overturn this

1. **Pruning is not modelled.** This is pure geometry, and assumes every lattice point
   exists. In a real run only survivors continue, so a signal migrating into a leaf that
   was thresholded away is still lost. That is a threshold-scheme interaction, not a
   tiling one — and a better tiling does not obviously fix it, since metric leaves get
   pruned too.
2. **Taylor basis only.** C4's unit-diagonal triangular `T` is the Taylor kinematic
   transform. `core/chebyshev.py` and the circular path transform differently, and the
   argument above does not carry over. If the paper's gap is real, this is where to
   look for it.
3. It contradicts the stated motivation of Kumar & Zackay (2026) §5.2.4 as this plan
   reads it, so it deserves an independent check before anyone acts on it.

## The sensitivity loss, quantified — the paper's actual open question

§5.2.4 closes with *"Further work is required to quantify the sensitivity loss and
compare basis strategies."* That is a different measurement from anything above, and it
is cheap. For a signal anywhere inside a leaf's cell, what phase error does it actually
incur, against the promised `eta / N_b`?

Two facts the self-checks pin first, both about what `eta` means:

- the **naive** per-axis grid is calibrated so a full-step offset costs exactly
  `eta / N_b` — the promise is kept, per axis;
- the **shipped** grid is coarser by exactly `2**(k-1)` on the order-`k` axis. That is
  the paper's own metric diagonalisation (§3.4, D38), deliberate, and it means `eta`
  does *not* mean "phase error below `eta / N_b`".

Sampling offsets uniformly inside each leaf, 268 s / 64 seg / `po=4`, in units of
`eta / N_b`:

| strategy | median stage | worst corner | median mismatch |
|---|---|---|---|
| `aggressive` | 3.35 | 16.80 | 0.0240 |
| `quadrature` | 3.58 | 14.38 | 0.0258 |
| `conservative` | 3.42 | 14.20 | 0.0251 |

**The three strategies are indistinguishable.** A typical signal costs ~3.4x the nominal
phase budget and ~2.5% in amplitude; a corner signal costs ~15x. Identical across
strategies that differ by **twenty orders of magnitude** in cost.

### Why, and it settles the question

`branch_param_padded` refines each parent to the criterion step `dparam_new` at every
stage: `num_points = ceil(dparam_cur / dparam_new)`, so the child cell is
`dparam_cur / num_points`, which lands just under `dparam_new` whatever
`dparam_cur` was. Measured cell half-widths at the final stage, against the criterion
step of `[9.7e-3, 0.163, 3.64, 122]`:

| strategy | final cell half-widths | `prod B(s)` |
|---|---|---|
| `aggressive` | `[1.10e-2, 0.123, 5.00, 96.4]` | 1.51e12 |
| `quadrature` | `[1.10e-2, 0.288, 4.02, 133]` | 1.10e17 |
| `conservative` | `[1.10e-2, 0.344, 6.99, 239]` | 1.48e32 |

All within a factor of ~2 of the criterion and of each other; the scatter is the integer
`ceil`. **The transported width does not set the cell — the criterion does. It sets only
how much redundant subdivision is done on the way.** So it can only buy cost, never
sensitivity.

That is the general statement, and it applies to the metric strategy too: a covering is
just another transport rule. It cannot improve sensitivity, and it measurably worsens
cost. Which closes the question this branch was opened to answer.

### What is worth reporting upstream

The `~3.4x` typical and `~15x` worst-corner excursion is a real, quantified answer to
the paper's request, and it is actionable: `eta` is optimistic by about a factor of
three for a typical signal under the recommended strategy, most of it the `2**(k-1)`
coarsening compounded across axes, the rest multi-axis corner addition. Anyone choosing
`eta` from equation `grid_criteria` alone will over-estimate their sensitivity. That
holds for the shipped default and does not depend on any of the metric work.

### Caveats

Deterministic geometry only. Pruning is not modelled, and it is the one place a wider
claimed region could still help — a signal near a cell edge may fall in a leaf that was
thresholded away, and a strategy claiming more territory keeps more such signals alive.
That is a threshold-scheme interaction, measurable only by injection, and nothing here
rules it out. Taylor basis only: the Chebyshev transform is not unit-diagonal
triangular, so neither D36 nor this section carries over to it.
