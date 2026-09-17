# 04_upstream_report.md — draft report for pravirkr/pyloki

Status: **draft, unposted, and awaiting human review.** Rewritten again on 2026-09-17
after session (s) replaced the withdrawn measurement with an exact one. History matters
here: an earlier version of this file claimed `tiling_strategy` buys no sensitivity
(refuted), then that `eta` is optimistic by 7.5x (retracted — that figure is a distance
to the containing cell's centre, not a covering radius). What follows is only what a
validated method now supports.

## What the exact search establishes

`nearest_template.py` computes the true minimum phase error to *any* leaf in the tree.
It rests on a structural fact — `branch_param_padded` derives a child's offset from
`(dparam_cur, dparam_new)` alone, which every leaf at a stage shares, so child offsets
are parent-independent and the leaf set is *exactly* a Minkowski sum of small per-stage
offset sets. 1e34 leaves live in a few hundred points per stage, and branch and bound
over that sum, with an admissible bound, never discards the optimum. It agrees with
brute-force enumeration to 1e-12 on every case small enough to enumerate, across all
three strategies, and its minimum is monotone in the search width.

Config: 268.4 s / 64 segments / `poly_order=4` / `N_b=64` / `eta=1`, 6 signal positions
x 15 stages = 90 cells.

1. **`aggressive`'s nearest template sits at 1.069 `eta/N_b`** (median; range 0.240 to
   2.115), exact in all 90 cells. So under the shipped default the closest template to a
   signal is a little over one tolerance away.

2. **`tiling_strategy` does affect sensitivity in the Taylor basis.** Seeding the
   incumbent with `aggressive`'s exact value makes the comparison a decision problem:

   | vs `aggressive` | strictly closer (proved) | not closer (proved) | unresolved |
   |---|---|---|---|
   | `quadrature` | **89** | 1 | 0 |
   | `conservative` | 28 | 0 | 62 |

   Median proven gain where closer: **>= 2.30x** for `quadrature` (a lower bound, since
   it comes from upper bounds on the other side). The mechanism is redundancy: the box
   strategies over-claim under transport, sibling cells overlap, and the signal is
   covered many times, so the template that scores is the nearest rather than the one
   whose cell contains it.

3. **The corner figure is not a covering radius.** The corner of the optimal grid costs
   exactly `(2^(k_max-1) - 1/2) * eta/N_b` — verified for orders 2-8, independent of
   `t_s` and `f_max` provided the corner is evaluated at the span that enters the step
   formula. That is 7.5x at `poly_order=4`. But it measures the distance to the
   *containing cell's* centre, and a search enumerates every grid point: minimising over
   neighbouring centres gives ~0.5. The `2^(k-1)` coarsening is economization working,
   not a defect, and `eta` is not optimistic by 7.5x. Nor is the corner a good proxy for
   sensitivity: the Chebyshev grid's corner is `k_max/2` (2.0 at `poly_order=4`, 3.75x
   better than Taylor's 7.5) and it still costs slightly *more* in S/N — see the basis
   section. Quote these closed forms as grid geometry, not as a sensitivity ordering.

4. **The amplitude cost is small.** Converting the phase geometry into S/N *without* a
   metric -- the folded profile is `P(k) * S(k)` with `S(k)` the characteristic function
   of the phase residual, scored with the shipped boxcar bank and cross-checked against a
   matched filter. Median fractional S/N loss, Taylor basis:

   | pulse duty | `aggressive` | `quadrature` | `quadrature`'s advantage (paired) |
   |---|---|---|---|
   | 0.05 | 5.34% | 2.26% | +1.56% |
   | 0.10 | 1.73% | 0.75% | **+0.47%** |
   | 0.20 | 0.40% | 0.17% | +0.11% |

   The last column is the **median over cells of the per-cell ratio**, which is the
   paired comparison and the one to use. Note it is *not* what dividing the two medians
   in the same row gives (+1.00% at 10% duty): the two losses are not co-monotone across
   cells, so the paired statistic is the smaller and the more honest of the two. Duty
   cycle is the dominant uncertainty -- an order of magnitude across the range -- and
   0.05 is marginal at `N_b = 64`, so the 0.10 row is the one to quote.

   So `aggressive` loses 1.73% at 10% duty to grid coarseness, of which about half a
   point is the tiling choice and the rest is irreducible at `eta = 1`. Both loss figures
   are upper bounds, since the template comes from a sup-norm search rather than a loss
   search, so the advantage is indicative; the proven part is the >= 2.30x above. Worth
   stating separately: a sup-norm phase error *overestimates* the loss by ~3x, because
   the residual attains its peak only briefly.

## What this does NOT establish

- **No proven amplitude *difference*.** The two losses in (4) are upper bounds, since the
  template comes from a sup-norm search rather than a loss search. The proven claim is
  the `>= 2.30x` phase gap in (2); the half-percent that follows from it is indicative.
- **Pruning is not modelled, and at this configuration it cannot be.** These are
  geometric distances to the nearest leaf; whether that leaf survives thresholding is a
  separate question. An injection campaign to settle it was designed, powered and then
  blocked: in the `quadrature` arm the candidate buffer never converges (doubling
  `max_sugg` from 2^18 to 2^20 doubles the candidate count each time, with saturation
  *rising*), because `prune_on_overload` relaxes the cut to keep the buffer full. So the
  scheme's thresholds are not the operative cut in that arm at any buffer tested, and the
  search cannot be run in a regime where the geometry is what decides survival. The
  geometry below stands; its conversion into detections is untested and, here,
  untestable. (Measured by another session; the buffer ratchet is theirs to report.)
- **The circular basis is untested.** Chebyshev is now measured (below); the circular
  transform is not.
- `conservative`'s 62 unresolved cells are absence of proof, not evidence against it.

## Where the gain sits — measured, not modelled

The gains in (2) are medians over stages 4 to 60. Which stages actually decide detection
is a question that was first answered with a model and then measured, and the two
disagree completely.

`survival_profile.py` drives the real pruning loop and records, at every prune level, the
minimum phase excursion from the injected signal to any surviving candidate — so
"is a covering leaf still alive?" is observed rather than inferred. Over 40 real runs (`aggressive`, S/N 14, `max_sugg` = 2^18):

| level | 1 | 10 | 20 | 30 | 40 | 50 | 63 |
|---|---|---|---|---|---|---|---|
| fraction with a covering leaf alive | 1.000 | **1.000** | 0.900 | 0.850 | 0.825 | 0.725 | 0.725 |

First-loss levels are 13, 13, 17, 19, 24, 27, 36, 41, 41, 44, 47 — **none before level
13**. A single-leaf survival model had put 99% of first losses in stages 1-10; the
measurement puts 0% there. The decision is made in the middle and late stages.

**What is and is not validated here.** The *shape* is robust: "no losses before level 13"
holds across two criteria (the excursion threshold above, and the metric mismatch
`m <= 1.0` this project already uses for recovery), two injected parameter sets and two
independent batches — 0/11, 0/15 and 0/13 losses by level 10. The *absolute level* is merely
unvalidated, which is weaker: the final-level alive fraction is 29/40 in one batch and
15/30 in another, and combined at 44/70 = 0.629 it is **consistent** with an
independently measured recovery rate of 38/50 = 0.760 (Fisher exact `p` = 0.16; pooled
82/120 = 0.683). It is unvalidated because one batch agreeing was never evidence, not
because the two disagree. Note also that this pipeline calibrates the injected amplitude
against each realisation's own noise, so runs are not exchangeable and binomial `p`-values
here are anti-conservative. The conclusions below rest on *where* losses occur, not on how
many.

Measuring the gain in that window:

| basis | stages 1-10 | **stages 13-47 (where leaves are lost)** | stages 32-60 |
|---|---|---|---|
| Taylor | 1.71x (47/60 cells) | **2.47x (54/54)** | 2.61x (48/48) |
| Chebyshev | 1.15x (47/60) | **1.57x (54/54)** | 1.72x (48/48) |

In the decision window `quadrature` is closer in **every one of 54 cells**, and the gain
is essentially the headline figure of (2). So the advantage is not concentrated where
nothing is at stake — it sits where the decision is made.

Checked under an independent criterion and a second signal, since the threshold above is
calibrated from the same runs: re-deriving "alive" as the metric mismatch `m <= 1.0` that
this project's recovery test already uses — nothing calibrated from the data — over two
injected parameter sets, N = 30 each, gives **0 of 28 losses by level 10** and an
earliest loss at level 13 in both. The centre of the window moves (median loss level 15
under the metric criterion, 27 under the excursion one, the former being the stricter
test) but the decisive claim does not.

The gain likewise survives every window definition, and so does the 100%:

| window | Taylor | Chebyshev |
|---|---|---|
| 13-47 (excursion criterion) | 2.47x (54/54 cells) | 1.57x (54/54) |
| 13-37 (metric criterion) | 2.09x (60/60) | 1.50x (60/60) |
| 13-16 (earliest quartile) | 1.74x (24/24) | 1.54x (24/24) |

Taylor's figure is window-dependent and should be read as **2.1-2.5x**; Chebyshev's is
stable at **~1.5x**. Remaining caveats: `aggressive` arm only (a `quadrature` profile is
blocked by the buffer non-convergence above), and N = 40 / 2x30.

## The Chebyshev basis

The same search ports to `poly_basis="chebyshev"` -- the branch there also goes through
`branch_param_padded`, and `shift_cheby_errors` ignores the values, so child offsets stay
parent-independent and the Minkowski-sum structure holds. Same 90 cells:

| basis | `aggressive` nearest template | `quadrature` closer | proven gain | `aggressive` loss | `quadrature` loss |
|---|---|---|---|---|---|
| Taylor | 1.069 | 89/90 | >= 2.30x | 1.73% | 0.75% |
| Chebyshev | **0.904** | 87/90 | >= 1.56x | **2.03%** | 1.14% |

(losses at 10% duty cycle; paired advantage +0.47% Taylor, +0.52% Chebyshev)

Chebyshev puts a ~15% closer template near the signal, and the tiling choice matters
*less* there, not more. Its nominal cell corner is also far better behaved --
`k_max/2 * eta/N_b` exactly, i.e. 2.0 at `poly_order=4` against the Taylor grid's 7.5,
linear in the order rather than geometric.

But note the nominal corner predicts the amplitude ordering **wrongly**: Chebyshev is
3.75x better on that measure and 15% better on sup-norm, yet costs slightly *more* in
S/N (2.03% against 1.73%). The Chebyshev residual oscillates across the whole domain
while the Taylor residual peaks only briefly at the window edge, and what smears the
profile is the residual's distribution, not its peak. On amplitude -- the thing that
matters -- the two bases are within ~0.3 points of each other.

## One practical point about `branch_max`

Max per-axis child count, against the shipped `branch_max = 16`:

| | `aggressive` | `quadrature` | `conservative` |
|---|---|---|---|
| Taylor | 7 | **28 — raises** | **27 — raises** |
| Chebyshev | 9 | 16 — at the limit | **20 — raises** |

The guard is strict (`num_points > branch_max`, `psr_utils.branch_param_padded:363`), so
16 builds and 17 does not. Three of the four non-`aggressive` configurations are
therefore unreachable at the default, and the Taylor gain above needs `branch_max`
raised. **Chebyshev + `quadrature` is the exception: it builds as shipped**, which makes
it the only configuration here with a proven gain that costs no config change.

That said, 16 of 16 is *at* the limit with no headroom. `num_points` is
`ceil(dparam_cur/dparam_new)`, which moves with `poly_order`, `tobs`, `eta` and the
segment count, so a nearby configuration turns it into a raised `ValueError` rather than
a slower search. It fits for this configuration; it is not robust.

## Why the advantage is not realisable

The geometry above is threshold-independent: no threshold, score or buffer quantity
enters it, and the set it enumerates is the full branching tree rather than the surviving
set. So the gains hold under any cut — and equally, they say nothing about which templates
survive one.

That distinction is decisive here, because `quadrature` cannot be run at its own
threshold scheme. `prune_on_overload` raises the effective cut to
`max(configured, top-K, median)` and never lowers it, and `quadrature`'s candidate count
is proportional to `max_sugg` over a factor of 8 in buffer size (no convergence), so the
ratchet is always active for it. Putting the two on the same scale, at the ladder
thresholds spanning the measured decision window:

| | score units |
|---|---|
| `quadrature`'s signal-template advantage (from the ~1% amplitude edge) | **+0.05 to +0.07** |
| cut elevation observed from the ratchet | **+3.4** (nominal 3.90 → effective 7.34) |

The advantage is real, exactly quantified, and about two orders of magnitude smaller than
the elevation it would have to overcome. Retention under the ratchet *is* score-ordered,
so a closer template is genuinely more likely to be kept — the sign favours `quadrature`
— but `K` is the buffer capacity and is the same for both arms while `quadrature` presents
9.1x more candidates at 2^18, so it meets the same rank cut over a far larger pool.

(The +3.4 figure is from a shipped-example configuration at `max_sugg = 2^10`, not from
the configuration above, so this is an order-of-magnitude comparison rather than a
measurement.)

So `tiling_strategy` is a real but small sensitivity knob, and `quadrature`'s ~1e5x
extra branching (and its higher recalibrated threshold ladder, 9.10 against 7.70 at
equal `P_d`) buys about half a percent of amplitude. That does not look like a trade
worth making, but it is a different statement from "the strategies are equivalent",
which is what an earlier draft of this file wrongly said.

Reproducers on https://github.com/assaferan/pyloki/tree/metric-gridding :
`docs/metric_gridding/nearest_template.py` and `amplitude_loss.py`, with
`tests/test_nearest_template.py`.

---

# ⚠ EVERYTHING BELOW THIS LINE IS THE WITHDRAWN EARLIER DRAFT

Kept only so the retractions in DECISIONS.md sessions (r) and (s) have something to
point at. **Its title and its central claim are both retracted** — the corner figure is
not a covering radius (D48), so `eta` is *not* optimistic by 7.5x. Do not quote, post,
or excerpt any of it. The live report is above.



**Title:** Taylor grid: the cell corner costs `(2^(k_max-1) - 1/2) * eta/N_b`, so `eta` is not the phase bound it looks like

Everything below is the **Taylor basis** (`poly_basis="taylor"`, the shipped default).
The Chebyshev path transforms differently and I have not measured it; see Scope.

Section 5.2.4 asks for the sensitivity loss of aggressive tiling to be quantified. The
largest term turns out not to be the tiling at all — it is in the grid criterion itself,
and it has a closed form.

### The corner of the optimal grid

Eq. `dk_criteria` promises `|dPhi| <= eta/N_b`. Eq. `dk_optimal` sets

    Delta d_k^opt = 2^(k-1) * (c/f_max) * (eta/N_b) * k!/t_s^k

Each axis is bounded *independently*, so a signal at the corner of the cell — offset by
half a step on every axis at once — incurs the sum. Axis `k` contributes

    (Delta d_k^opt / 2) * (f_max/c) * t_s^k/k!  =  2^(k-2) * eta/N_b

and summing `k = 1..k_max` telescopes:

    corner phase error = (2^(k_max-1) - 1/2) * eta/N_b

One hypothesis is doing work there and is worth stating: the corner is evaluated at the
**same span `t_s` that enters the step formula**. In general axis `k` contributes
`2^(k-2) * (eta/N_b) * (t/t_s)^k`, so it is only at `t = t_s` that the `t_s` and `f_max`
dependence cancels and the sum telescopes. That is the relevant case for a grid quoted
over its own span, and it is also why the figure is larger under a moving reference
frame, where a leaf's validity window is displaced from its own epoch and `t/t_s > 1`.

Verified against `psr_utils.poly_taylor_step_d_vec` for `k_max = 2..8`, exact to machine
precision and **independent of `t_s` and `f_max`** (both cancel):

| `poly_order` | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|
| corner, in `eta/N_b` | 1.5 | 3.5 | 7.5 | 15.5 | 31.5 | 63.5 | 127.5 |

At the default `poly_order = 4` the corner costs **7.5x** the nominal tolerance; at
`poly_order = 5`, which is what a circular-orbit search runs, **15.5x**. The `2^(k-1)`
coarsening is sound per coefficient — it comes from `|T_k| <= 1` in the Chebyshev basis
(appendix `app:optimal_gridding`) — but applied as independent bounds on monomial `d_k`
cells it is the same factor that makes the corner sum grow geometrically.

So `eta` regulates grid density, not the worst-case phase error, and the gap widens by
`2x` per polynomial order. Anyone reading `eta = 1` as "at most one bin of drift" is
optimistic by `7.5x` at `poly_order=4`. The amplitude cost is second order and therefore
much gentler (below), but the phase bound is not what the symbol suggests.

### What the tiling strategies actually buy

Two things that were not obvious to me, both measured at 268.4 s / 64 segments /
`poly_order=4` / `N_b=64` / `eta=1`:

The final cell is set by the **criterion**, not by what transport delivered.
`branch_param_padded` (`utils/psr_utils.py:360-368`) uses
`num_points = ceil(dparam_cur/dparam_new)` and returns `dparam_cur/num_points`, so a
branched axis lands just under `dparam_new` whatever `dparam_cur` was; the complementary
case in `core/taylor.py:148-154` leaves an axis unbranched while its shift is below
`eta`. Measured final half-widths agree within the factor of two that the `ceil` and
that guard allow, across strategies differing by `10^20` in `prod B(s)`.

But the strategies still differ in sensitivity, through **redundancy** rather than cell
size. `quadrature` and `conservative` over-claim under transport, so sibling cells
overlap and the signal is covered many times; the template that scores is the *nearest*,
not the one whose cell contains it. Following the true signal down the tree (12 random
signal positions, seeding an 81-cell base block as the real search does):

| strategy | nearest-template phase error | worst signal | leaves covering the signal | mismatch |
|---|---|---|---|---|
| `aggressive` | 1.89 | 3.09 | 1 | 0.0088 |
| `quadrature` | 0.617 | 1.07 | ~3 700 | 0.0021 |
| `conservative` | 0.615 | 1.35 | ~13 000 | 0.0021 |

`aggressive` partitions exactly (multiplicity 1), and is the one strategy that does not
deliver `eta/N_b` even in this averaged sense. The redundant strategies do, and the
phase-error gain is `3.1x` (per-signal 1.8-4.9x). In amplitude, though, that gain is
about **0.34% in S/N** (0.44% loss against 0.105%) — second order, so three orders of
magnitude of extra branching buy a third of a percent.

### The one actionable item

**`conservative` is strictly dominated by `quadrature`.** Same nearest-template error
(0.615 vs 0.617), same mismatch (0.0021), for `2^15` times the cost
(`prod B(s)` 1.48e32 vs 1.10e17). If `conservative` is kept as an option it is worth
saying in the docstring that it buys nothing over `quadrature`.

### Scope

- **Taylor basis only.** `T(delta_t)` is lower-triangular with unit diagonal here, which
  is what makes the box strategies' cells tile and the corner analysis clean. The
  Chebyshev interval-change matrix mixes orders off-diagonally, so none of the tiling
  results above transfer to `poly_basis="chebyshev"`, and I have not measured that path.
  The closed form for the corner is basis-independent in derivation but is stated for the
  Taylor kinematic grid.
- **Deterministic geometry.** The nearest-template numbers are geometric; they do not
  model pruning. The 0.34% amplitude edge is far smaller than the difference between the
  recalibrated threshold ladders (top threshold 7.70 for `aggressive` vs 9.10 for
  `conservative` at equal `P_d`), so I would not expect it to survive a real search, and
  no injection run at reachable statistics could resolve it either way.
- Mismatch figures are second-order estimates from a parameter-space metric built for
  this work, not from folded S/N.

Reproducers on https://github.com/assaferan/pyloki/tree/metric-gridding :
`docs/metric_gridding/sensitivity_loss.py` (corner closed form, self-checks) and
`docs/metric_gridding/pruning_multiplicity.py` (nearest-template and multiplicity).
Happy to open a PR adding the closed-form check as a test — it is three lines and it
pins eq. `dk_optimal` against the code.
