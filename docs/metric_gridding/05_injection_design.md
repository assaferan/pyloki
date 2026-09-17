# 05_injection_design.md — pre-registered design for the injection campaign

Status: **design, pre-registration. No production run has been made.** Written before
looking at any outcome, because three results on this branch were withdrawn after the
fact and a pre-committed design is the cheapest available protection against a fourth.
Approved by assaferan 2026-09-17.

## Hypothesis

Tiling redundancy keeps signals alive through pruning. A signal away from a leaf centre
scores lower; if it falls below the threshold at any stage the leaf is discarded and the
signal is lost outright. A strategy that over-claims under transport puts a nearer
template next to the signal (D52: `quadrature` strictly closer in 89/90 cells, gain
`>= 2.30x`), so it should lose fewer signals near the survival boundary.

**This is not the 0.5% amplitude effect of D55.** A thresholded-away leaf loses the
signal entirely, so the effect is bimodal: ~0 for signals comfortably above the ladder,
total for signals near it. Powering against a mean amplitude shift would understate it,
and that error is why an earlier session wrongly called the campaign futile.

## Arms, and which question they answer

`branch_max` defaults to 16 and the guard is strict, so (D62):

| basis | `aggressive` | `quadrature` |
|---|---|---|
| Taylor | 7 — builds | 28 — **raises** |
| Chebyshev | 9 — builds | 16 — builds, at the limit |

- **Arm A (shipped behaviour): Chebyshev, `aggressive` vs `quadrature`.** The only pair
  that runs at the default. This is the question a maintainer would ask.
- **Arm B (raised limit): Taylor, `aggressive` vs `quadrature` at `branch_max >= 28`.**
  Larger measured effect, but explicitly not a shipped configuration.

Run A first. Note A is the *harder* test: in Chebyshev the tiling gap is smaller
(`>= 1.56x` against `>= 2.30x`), so if only one arm can be afforded, A risks a null that
B would not have produced. Say which question is being answered, in the write-up.

## Pairing

Identical injected signal and identical noise realisation across arms, with
`tiling_strategy` the only difference. The test is then McNemar on discordant pairs
rather than two independent `P_d` estimates, which is the difference between ~300
injections and tens of thousands. **Prerequisite:** the pipeline must reproduce a noise
realisation across arms; if it cannot be seeded to do so, fixing that is task one.

## Stratification, and the built-in null

Inject at controlled positions relative to the leaf lattice, using
`nearest_template.py` to place signals by *measured* distance to the nearest template
rather than by nominal cell coordinate:

- **centre** — at the nearest template. **This is the null stratum: the arms must agree.**

  **REVISED with the corrected table (inherited, not verified here). The effect stratum
  does not exist.** Decision-weighted `pi` over 24 positions, corrected against the
  profile-overlap version that was in use:

  | `rho_AB` from | min | q25 | median | q75 | max | pi<0.55 | pi>0.70 |
  |---|---|---|---|---|---|---|---|
  | profile (was in use) | 0.440 | 0.595 | 0.645 | 0.704 | 0.741 | 12% | 29% |
  | corrected | 0.470 | 0.556 | **0.590** | 0.619 | **0.643** | 21% | **0%** |

  So `pi > 0.70` is unreachable: **three strata collapse to two**, and the primary test
  runs on the mid stratum alone. The null stratum roughly doubles in incidence (~1 in 5
  rather than ~1 in 8), which makes the design's most important control cheaper to
  screen for. Their control for the substitution: the profile row reproduces the cached
  per-position values to 5e-5 on all 24 positions.
- **face** — half way to the nearest template.
- **corner** — at the maximum distance the strategy permits.

If the centre stratum shows a strategy difference, the measurement is picking up
threshold-calibration mismatch rather than coverage, and that is known *before*
publication. A single aggregate `P_d` comparison has no such check, which is how the
0.34% and the 1.89 got through.

## Controls

- **positive** — signal exactly on a base-grid point: all strategies must recover it, so
  a null here means the pipeline cannot detect agreement.
- **negative** — signal placed outside any claimed cell: all strategies must fail, so a
  null here means the pipeline cannot detect failure.

Without both, a null result is uninterpretable.

## Thresholds

Recalibrate per strategy at equal **on-grid** `P_d`, or the comparison measures the
ladder rather than the tiling: the cached ladders top out at 7.70 (`aggressive`) against
9.10 (`conservative`), which would swamp any coverage effect. Only `aggressive`,
`conservative` and `metric` ladders are cached; **`quadrature` and both Chebyshev
variants must be calibrated first**, at whatever `branch_max` the arm uses.

## Power, and the n it implies

Model (`injection_power.py`): the leaf accumulates segments, its score at stage `s` is
`a (1 - loss_s) sqrt(s+1) + (1/sqrt(s+1)) sum_{i<=s} n_i` with per-segment noise shared
between arms, and the signal is recovered iff it clears the calibrated threshold at every
stage. `loss_s` is the measured per-stage S/N loss of the nearest template. Sanity check:
the model gives `P_d` 0.05 / 0.18 / ~0.45 at S/N 10 / 12 / 15, against the Phase 3
injections' 0/3, 1/3, 3/3 — consistent within small-sample noise.

> **The table below is retired.** Its `P_d` column comes from `injection_power.survival`,
> which has since been falsified against a converged measurement: at S/N 14 and 2^18
> (where `aggressive` converges, saturation 0.064) the real search recovers 38/50 = 0.760
> against the model's 0.368, `P(X >= 38 | 0.368) = 1.8e-08`. The model is ~4.5 S/N
> **pessimistic**, so its absolute `P_d`, its operating point and its pair counts are all
> unreliable. Use the 954 pairs computed against a *measured* `p_disc` instead. Kept for
> the record, and because the `r`-parameterised structure is still the right shape.

Pairs for 80% power at alpha = 0.05, against the loss ratio `r = loss_quad / loss_agg`
(`r = 1` is no effect; measured `r ~ 0.34`, IQR 0.21-0.45):

| S/N | `P_d` (`aggressive`) | r=0.21 | r=0.34 | r=0.45 | r=0.70 |
|---|---|---|---|---|---|
| 8 | 0.005 | 3 536 | 4 333 | 5 145 | 8 731 |
| 10 | 0.050 | 524 | 630 | 748 | 1 445 |
| 12 | 0.178 | 273 | 323 | 383 | 681 |
| **14** | **0.369** | **243** | **288** | **343** | **596** |
| 16 | 0.567 | 282 | 336 | 410 | 713 |
| 20 | 0.848 | 610 | 725 | 867 | 1 588 |

**REVISED — the affordable figure was wrong, for the same reason the report's headline
gain was.** `r = 0.34` is a median over all stages, but 99% of first threshold failures
happen in **stages 1-10**, and there `r = 0.675`: the advantage is much weaker where
detection is actually decided, and `aggressive`'s own loss there is 1.03% rather than
2.05%. Applying the measured ratio only in the window that matters:

| S/N | `dP_d` | discordance | pairs for 80% power |
|---|---|---|---|
| 10 | +0.08 pt | 0.08% | 5 011 |
| 12 | +0.27 pt | 0.27% | 1 430 |
| **14** | **+0.50 pt** | **0.50%** | **772** |
| 16 | +0.53 pt | 0.53% | 723 |

So **~770 pairs for Arm B**, not 290. The **S/N 15-17 recommendation is withdrawn.**

It came from fitting a constant offset to a single pilot point, and a second measurement
falsified it: against the real search my model is *high* at S/N 12 (0.177 vs 0.125) and
*low* at S/N 14 (0.368 vs 0.600 at 2^14), so the discrepancy changes sign and the
measured curve is 2.0x steeper rather than shifted. An offset fitted to one point cannot
tell those apart. **The operating point must come from the measured on-grid `P_d` ~ 0.5
criterion at the final buffer, which brackets below S/N 14**, not from this model's
discordance peak. Arm A (Chebyshev) is weaker again —
its early gain is 1.15x against Taylor's 1.71x — so it needs materially more than that,
and the figure should be recomputed on the Chebyshev ladder before committing.

And ~770 is a **floor**, not an estimate: this model treats the arms as differing only by
a deterministic `loss_s`. They also score *different templates*, which adds a stochastic
term that a peer session measures at **sd 0.230** against a deterministic amplitude
difference of 0.044 — noise roughly **5x** the signal. (The 0.133 first reported came
from a profile-overlap correlation that has since been corrected to a score correlation;
`sqrt(2 x 0.0265) = 0.230`.) Pairing on the noise realisation does
not remove it. Carried into McNemar that shrinks the discordance asymmetry and raises `n`,
plausibly by an order of magnitude. **That term must be in the model before any n is
committed to.** Their `injection_power.py:power()` already implements the aggregation,
and `rho_check.py:propagate()` reproduces the uncorrected `pi` exactly as a control, so
this is a substitution rather than a re-derivation.

One caveat inherited with the deterministic term is now **partly retired**: it assumes
`L_A` and `L_B` are upper bounds whose looseness need not be symmetric. Over stages 1-10
-- the only window that decides detection -- both are *exact* (identical at 4M and 16M
nodes, searches completing in 0-2 s), so no asymmetry can enter where it would matter.
What that check did reveal is that `quadrature` is genuinely **worse** at stages 2 and 8
(1.667 vs 1.027, 2.067 vs 1.657), so the deterministic term is not sign-definite early
and the design must expect two-sided discordance there.

Two caveats, both one-sided in the safe direction. `quadrature`'s per-stage losses are
*upper* bounds from a budget-limited search, so `r` is overestimated, the effect
underestimated and `n` overestimated. And the numbers above use the Taylor ladder and
Taylor losses; Chebyshev's loss ratio is nearer 0.7, so **Arm A needs ~600 pairs, not
~290**.

## Affordability (inherited, unverified here)

At the operating point the design selects -- on-grid `P_d ~ 0.5` in both arms, where
discordance is maximal and the pilot measured `p_disc = 0.30` -- the campaign needs
**340 to 2 840 pairs**, quoted as 954 under one stage weighting. The `p_disc` is measured;
the **stage weighting is not**, and it rests on a modelled per-stage survival curve from
the same family as the model falsified in section "Power". Re-running the aggregation
under flat / `sqrt` / late-only / early-only weightings spans a factor of 8 in `n`. Treat
954 as a point in that range, not the requirement.
Powering instead against the `p_disc = 0.08` seen at S/N 12 would demand ~3600 pairs,
but that is a regime the design avoids by construction, so it is not the requirement.

`n = 2000` covers `p_disc >= 0.145`, which is no longer the whole measured range, so the
3-point operating-point pilot at the final `max_sugg` is a **precondition**: below 0.145,
`n` rises to ~3600 or the design is re-scoped.

## Pre-committed analysis

- **Test**: McNemar (exact binomial on discordant pairs), two-sided, alpha = 0.05.
- **Primary endpoint**: recovery of the injected signal, corner stratum, Arm A.
- **n**: 600 pairs for Arm A, 300 for Arm B, at S/N 14.
- **Null stratum**: centre, same test; a significant result there invalidates the primary.
- **Stopping rule**: fixed n. No interim looks, no extension on a near-miss.
- **Effect of interest**: `dP_d >= 1` point. Below that the campaign cannot distinguish
  the effect from the model's own simplifications.

## What the model ignores

Competing candidates and the beam (this is survival of the true leaf, not the whole
search), a single duty cycle, Gaussian noise, and the assumption that both arms are
calibrated to equal on-grid `P_d`. Each would need arguing separately; none is a reason
to defer the run, but all belong in the write-up.

## A null is a good outcome

If coverage does not matter once pruning is in play, that closes the question this branch
was opened on, and it is publishable. The campaign does not have to find an effect to be
worth its budget — it has to be able to detect one if it is there, and be believed either
way.
