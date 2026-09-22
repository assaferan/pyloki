# 05_injection_design.md — the `tiling_strategy` injection campaign: design, and why it cannot be run here

> ## ⚠ The campaign was NOT run, and at this configuration it CANNOT be
>
> Written across two sessions. The second one completed the two experiments the first
> left open and, in doing so, established that the campaign this document designs is not
> executable at Phase 3's configuration. **Read §1, then §6.6.**
>
> One thing did ship: the silent-ratchet defect found along the way is
> [PR #14](https://github.com/pravirkr/pyloki/pull/14) upstream, written up in
> `06_upstream_max_sugg_logging.md`.
>
> **This document contains retracted claims that are kept visible rather than deleted**,
> because the branch's failure mode is a correct measurement carrying a framing that
> does not hold. Struck-through text and ⚠/⚑ boxes mark them. §8 lists the three
> conclusions that were withdrawn; §10 is the assumption register and now has **no**
> independently verified row.
>
> No inherited library code was modified by the design work. The one library change,
> the logging patch, was built on a separate branch off `upstream/main`.

Arm analysed: **Chebyshev basis, `aggressive` vs `quadrature`, at the shipped
`branch_max = 16`** — the only pairing with a proven geometric gain (D60) that builds at
the default (D62). Config throughout is Phase 3's: 268.4 s / 64 segments /
`poly_order = 4` / `N_b = 64` / `eta = 1`.

---

## 1. Verdict

**The campaign cannot be run at this configuration, and the obstacle is not cost.** The
comparison requires both arms to run against the threshold scheme rather than the
candidate buffer. `aggressive` can; `quadrature` cannot, at any buffer size, and raising
the buffer does not help because the candidate count rises with it.

### What blocks it

**`quadrature`'s candidate set never converges (§6.6).** On one fixed set of
realisations swept across four buffers:

| arm | 2^18 | 2^19 | 2^20 | 2^21 | saturation at 2^21 |
|---|---|---|---|---|---|
| `aggressive` | 16 778.5 | **16 778.5** | **16 778.5** | **16 778.5** | 0.008 |
| `quadrature` | 152 007 | 390 575 | 808 116 | **1 610 052** | 0.768 |

`aggressive` is done by 2^19 — its **per-realisation** counts are bit-identical across
2^19–2^21, so the buffer is provably irrelevant there, not merely indistinguishable.
(The 2^18 median coincides but one realisation is still clipped; see §6.6.) Confirmed
directly in §6.7: `aggressive` ratchets on **0 of 192** levels while `quadrature`
ratchets on 60, with the effective cut roughly double the nominal. `quadrature`'s count
grows by a factor of **10.6** over a factor of 8 in buffer — span exponents 1.135, 1.022,
0.994 — with saturation pinned near 0.77. The overflow ratchet relaxes the cut to keep
the buffer full at whatever size it is given, so there is no escape by spending more.

**Consequence: in the `quadrature` arm the scheme's thresholds have never been the
operative cut**, at any buffer used on this branch. The equal-`P_d` recalibration (§6.3)
is not doing what §6.3 says in that arm.

### What was measured before that became clear

- **The paired `max_sugg` experiment (§6.5, 50 realisations, all four cells,
  pre-registered before the first cell ran).** Two answers that must not be collapsed:
  - *Buffer **pressure** is arm-dependent, significantly* — at 2^18, 10/50 `quadrature`
    runs exceed 0.9 saturation against 1/50 for `aggressive`, paired McNemar **p = 0.012**.
  - *Whether it biases **outcomes** is unresolved, and the point estimate is worse than
    the effect being measured* — the pre-registered sign test gives p = 0.34, a failure
    to reject and **not** a null, with mean Δ = **+0.080 [−0.043, +0.203]** against a
    campaign effect of interest of **0.06**. Closing it needs n ≈ 211 per cell.
  - Raising the buffer produced **20 recovery flips, all toward recovery, none away.**
- **`rho_AB` was the wrong quantity (§5.2).** The code used the *profile* overlap where
  the model needs the *score* correlation — the same smearing factor weighted by the
  pulse's spectrum rather than the filter's, and the selected boxcar has 4.10x the
  pulse's `<k²>`. Median `1 − rho` 0.0088 → **0.0265**; `pi` 0.697 → **0.628**.
- **The stratification inherited that error (§5.1).** Corrected, the **`effect` stratum
  does not exist** — no sampled position reaches `pi` > 0.70 — so §9's three strata
  collapse to two. The `null` stratum, the design's main control, roughly doubles to
  ~1 position in 5.

### What was withdrawn

- **§4.3's "the decision is made where the geometry isn't" is backwards.** The measured
  per-stage survival profile puts the first loss at level 13 and none before it, against
  a modelled 99% of losses by stage 10. This was *the* structural reason the effort
  expected a small effect. The per-stage `pi` profile survives (it is geometry, with no
  survival model in it); the weighting over it does not.
- **The equal-`P_d` ladder check (§10 row 8).** Batch A's 3/24 agreement ran at 2^14
  where both arms saturate, so it cannot be credited to the ladders. The register now
  has no verified row.
- **The "defensible upstream sentence" about tiling (§8).** Nothing in this document
  currently supports one.

### `n`, and why it is a range

`n` is **368–691** and the precise value is not established. The early-decision corner
(2 842 pairs) is excluded by the measured survival profile; the spread that remains is
the profile's own, whose two criteria disagree on where the survival curve falls. The
survival profile is **inherited from the `metric-gridding` session, uncontradicted but
not validated** (§4.3 ⚑) — its shape is robust across two criteria, two parameter sets
and three batches; its absolute levels are not. Earlier figures of 954 and 368 in this
document are superseded and marked.

### What would unblock it

1. ~~**A stricter ladder.**~~ **Tested and it fails — §6.7.** A tenfold reduction in the
   detection target buys 4 ratcheted levels out of 192 and leaves zero clean runs,
   because the effective cut is set by `top-K`, not by the scheme: raising the ladder by
   0.87 moved the cut that actually ran by **−0.06**. No `P_d` target fixes a cut that
   is not set by `P_d`.
2. **The user-facing question instead.** "What does a user get at the shipped default,
   ratchet included" is answerable now with §6.5's apparatus — but it is a claim about
   pyloki-as-shipped, not about tiling, and must be reported as one.
3. **A smaller configuration.** Tests conversion in a regime nobody deploys; the
   geometry is already proven at the deployed one.
4. **A different buffer policy (§6.7, out of scope here).** The ratchet exists because
   `max_sugg` is a hard cap enforced by discarding the lowest-scoring candidates. A
   policy that refused to proceed, or reported the shortfall as a first-class output,
   would let `quadrature` fail honestly instead of silently substituting a different
   cut. A library change, and the natural successor to
   [PR #14](https://github.com/pravirkr/pyloki/pull/14) — that made the substitution
   visible; this would make it optional. **The only remaining route to the idealised
   question at the deployed configuration.**

---

## 2. What is being powered against, and why not a 0.5% mean shift

Everything established so far (D50–D62) is deterministic geometry: the distance from a
signal to the nearest leaf, and the S/N that distance costs. None of it models whether
that leaf **survives thresholding**, which is the open question.

The mechanism is bimodal, not a mean shift. A signal whose covering leaf is cut at some
stage is lost entirely; one that clears every threshold is essentially unaffected. So
the quantity to power against is the discordance structure of a **paired** experiment:

    p_disc = P(recovered_A ≠ recovered_B),        the discordance rate
    pi     = P(quadrature wins | discordant),     the effect size, null = 1/2

and the test is **McNemar**, which conditions on the discordant pairs and reduces to a
binomial test of `pi = 1/2`. Powering against the 0.5% mean amplitude difference would
be the wrong calculation twice over: it would assume the effect is a shift when it is a
flip, and it would ignore that the two arms' scores are correlated (§3, M4), which is
what actually sets the effect size.

---

## 3. The score model, stated in full

This is where a hidden assumption would do the most damage, so it is written out rather
than left in code. The implementation is `injection_power.py`; these are its M1–M6.

**M1 — the score is a matched-filter S/N in units of the noise sigma.**
`scoring.snr_score_batch_func` whitens per bin (`ts_e/√ts_v`, `scoring.py:287`) and
correlates with a boxcar that is exactly zero-mean and unit-L2-norm
(`scoring.py:136-139`: `width·h = size_w·b` and `width·h² + size_w·b² = 1`). Under H0 a
per-leaf score is therefore N(0, 1) and the `threshold_scheme` values are directly
comparable to it. Read off the implementation; not an assumption about the physics.

**M2 — accumulation.** `common.shift_add_batch` adds `ts_e` and `ts_v` linearly across
segments, so after `s+1` of `nseg` segments the on-grid mean score is

    mu(s) = snr_final · sqrt((s + 1) / nseg).

This is exactly the model the shipped ladder is built on — `schemes.bound_scheme` sets
S/N² linear in the stage index — so the power calculation stands on the same footing as
the thresholds it is powered against.

**M3 — THE PHASE-ERROR-TO-SCORE MAPPING.** A phase residual does not attenuate the
pulse, it smears it. `amplitude_loss.snr_ratio` computes the smeared folded profile
exactly, as `P(k)·S(k)` with `S(k)` the characteristic function of the residual
distribution, and scores it with the pipeline's own filter bank (D53–D54). Writing
`L_X(s)` for the resulting fractional S/N loss of arm X's best covering leaf,

    score_X(s) = mu(s) · (1 − L_X(s)) + N_X(s),     N_X ~ N(0, 1).

Multiplicative in amplitude, additive in noise: the loss is a loss of *recovered
amplitude*, while the noise level is unchanged because folding with a different
ephemeris averages the same samples. **A sup-norm phase error must not be substituted
for `L`** — D56 measured that doing so overestimates the loss by ~3x, because the
residual attains its peak only briefly.

**M4 — the two arms' noise is correlated, and that correlation is the effect's
denominator.** This is the point a naive paired power calculation gets wrong. Pairing
fixes the data; it does **not** make the comparison noise-free, because the two arms
score *different templates*:

    D(s) = score_B(s) − score_A(s) = mu(s)·(L_A − L_B) + eps,
    eps ~ N(0, 2·(1 − rho_AB(s))),

with `rho_AB` the overlap of the two arms' templates. `rho_AB` is computed with the same
smearing machinery applied to the phase residual **between** the templates,
`delta_A − delta_B`; that is legitimate because the residual is linear in the
coefficient vector. Measured (§4), the stochastic term is **three times larger than the
deterministic one**, which is the single most consequential number in this document.

**M5 — single crossing.** A pair counts as discordant when, at some stage, the threshold
separates the two arms' scores. The model sums the per-stage marginal probability
weighted by the probability of having reached that stage. That double-counts pairs
marginal at more than one stage, so the modelled `p_disc` is an **upper** bound and is
anti-conservative for n. Measured against the pilot it is high by ~4x (§6).

**M6 — one covering leaf per arm per stage.** Survival is really a max over all leaves
covering the signal, and `quadrature` has far more of them. Ignoring the max
**understates** `quadrature`. It is deliberate: more leaves is a *more trials* channel,
not a *closer template* channel, and it is exactly what the equal-`P_d` recalibration
exists to null out. It is also why the recalibration is not optional (§7).

---

## 4. Measured inputs

### 4.1 Geometry, from the exact search

`injection_power.py --geometry`, Chebyshev basis, 6 signal positions × 15 stages, 83 of
90 cells (7 dropped because branch-and-bound hit its node budget and a budget-limited
value is only an upper bound, per the convention fixed in D52).

| ducy | median `L_aggressive` | median `L_quadrature` | paired advantage | median `1 − rho_AB` | `sigma` of `eps` |
|---|---|---|---|---|---|
| 0.05 | 7.19% | 4.27% | +2.79% | 0.0346 | **0.263** |
| 0.10 | 1.90% | 1.10% | **+0.73%** | 0.0088 | **0.133** |
| 0.20 | 0.48% | 0.28% | +0.18% | 0.0022 | **0.066** |

(Matched-filter scoring throughout; `04_upstream_report.md` quotes the boxcar figures,
2.03% / 1.14% / +0.52% at ducy 0.10, over a different stage set. The two agree to
within the filter difference and neither is quoted in place of the other.)

The last column is the new quantity. At ducy 0.10 the deterministic part of `D` is
`mu · 0.0073` — about **0.044 at mu = 6** — against a stochastic part of sd **0.133**.
The template-to-template noise difference is 3x the amplitude difference. That is why
`pi` comes out near 0.6–0.7 rather than near 1, and it is why pairing on the noise
realisation buys much less than it looks like it should.

> **The `1 − rho_AB` and `sigma` columns above are the uncorrected ones (§5.2).** At
> ducy 0.10 the score correlation gives `1 − rho_AB` = **0.0265**, so `sigma` is
> **0.230**, not 0.133, and the ratio to the deterministic term is **~5x, not 3x**. The
> ducy 0.05 row barely moves (`<k²>` ratio 1.03) and 0.20 moves more than 0.10 did; only
> 0.10 has been recomputed in full. The qualitative point — that the template-to-template
> noise dominates and pairing buys less than it looks like — is unchanged and
> strengthened.

### 4.2 Thresholds, recalibrated per strategy on the Chebyshev pattern

The cached schemes under `schemes/` are Taylor (`kind="poly_taylor_moving"`) and contain
no `quadrature`, so they are **not** valid here and were regenerated:
`schemes/cheby_aggressive.npz`, `schemes/cheby_quadrature.npz`, Viterbi-optimised at
`P_d = 0.1`, `snr_final = 10`, `kind="poly_chebyshev_moving"`.

| | branching product | log₂ complexity | top threshold | achieved `P_d` |
|---|---|---|---|---|
| `aggressive` | 1.51e12 | 9.16 | 7.10 | 0.1031 |
| `quadrature` | 5.41e20 | 27.03 | 7.00 | 0.1031 |

Two things worth noting. The ladders are **nearly identical** (7.10 vs 7.00), unlike the
Taylor pair (7.70 vs 9.10) — one fewer confound, and a point in this arm's favour.
`quadrature` pays 2^17.9 ≈ 2.4e5 times the work for it.

### 4.3 ~~Where the decision is actually made — and it is not where the geometry is~~ — **REVERSED, see the box below**

The cumulative H1 success of the Viterbi ladder drops from 1.0 to 0.1031 over stages
1–30 and is **flat thereafter**. So the survival decision is made early, where
`mu(s) = 10·√((s+1)/64)` is still below the threshold:

| stage | 2 | 6 | 10 | 14 | 20 | 28 | 40 | 50 |
|---|---|---|---|---|---|---|---|---|
| threshold | 3.10 | 4.00 | 4.30 | 5.10 | 5.80 | 6.20 | 6.60 | 6.70 |
| `mu(s)` | 2.17 | 3.31 | 4.15 | 4.84 | 5.73 | 6.73 | 8.00 | 8.93 |
| reach `P(survived)` | 0.742 | 0.238 | 0.177 | 0.149 | 0.118 | 0.103 | 0.103 | 0.103 |
| `Delta = mu·(L_A−L_B)` | −0.000 | 0.000 | 0.007 | 0.044 | 0.032 | 0.064 | 0.129 | 0.123 |
| `sd` of `eps` | 0.000 | 0.003 | 0.158 | 0.158 | 0.165 | 0.155 | 0.190 | 0.161 |
| `pi` at this stage | 0.515 | 0.533 | 0.558 | 0.721 | 0.677 | 0.748 | 0.928 | 0.866 |

**This is the structural reason the effect is small.** The stages where the tiling
geometry is largest (`pi ≈ 0.87–0.93` at stages 40–50) are stages where nothing is being
decided any more — `reach` has already flattened at 0.103. The stages that decide
survival (2–14) are the ones where the two strategies have barely diverged under
transport, and `Delta` there is 0.000–0.044 against an `eps` of 0.158. The effect is
suppressed twice at the same time: `mu` is small *and* `L_A − L_B` is small, and for the
same reason — few segments have been accumulated.

> ⚠ **The `reach` row — and therefore the whole "where the decision is made" claim — is
> a MODELLED survival curve, and the model family is now measured to be badly wrong.**
> `reach` is the ladder's own cumulative H1 success. The `metric-gridding` session's
> survival model, of the same single-leaf Gaussian family, predicts `P_d` = 0.368 at
> S/N 14 where §6.5 measures **0.760** in the converged `aggressive` arm — pessimistic
> by ~4.5 S/N, `p` = 1.8 × 10⁻⁸ — and it has traced the gap to M6: thresholds
> calibrated for the whole branching pattern applied to a signal given exactly one leaf,
> when `aggressive` had 16 778 candidates alive. If survival is really a max over many
> near-covering leaves, `reach` decays far more slowly than the row above, and the
> flat-after-stage-30 premise fails.
>
> What survives and what does not:
> - **The per-stage `pi` row survives.** It is geometry plus §5.2's corrected `rho`, and
>   `reach` does not enter it. `pi` really is 0.515–0.558 early and 0.87–0.93 late.
> - **The claim that the early stages are the ones that matter does not survive.** That
>   is `reach`, and `reach` is the modelled quantity.
> - So the section's title is not established. The effect may be small for the reason
>   given, or the decision may extend into the stages where `pi` is large — which would
>   make it *bigger*, not smaller.
>
> This propagates: `DECISION_STAGES = [2, 6, 10, 14, 20, 28]` in `injection_power.py`
> is chosen from this premise, so §5.1's stratification and §9's strata are conditioned
> on it, and §9's `n` moves by a factor of 8 across plausible reweightings (§10 row 8b).
>
> ### ⚑ Now measured, and the premise is not merely wrong — it is backwards
>
> The `metric-gridding` session measured the real profile at assaferan's request
> (its `1f8a7e8`): 40 runs, `aggressive`, S/N 14, 2^18, recording the minimum phase
> excursion from the signal to any survivor at each prune level. **Inherited and not
> re-run here**; §11.1's rule applies.
>
> | level | 1 | 10 | 20 | 30 | 40 | 50 | 63 |
> |---|---|---|---|---|---|---|---|
> | fraction alive | 1.000 | 1.000 | 0.900 | 0.850 | 0.825 | 0.725 | 0.725 |
>
> First-loss levels are 13, 13, 17, 19, 24, 27, 36, 41, 41, 44, 47 — median 27, **none
> before 13**. So **0% of losses occur by stage 10**, where the modelled curve put 99%.
>
> **Status of the instrument: uncontradicted, not validated.** An early claim that its
> first batch validated it (29/40 = 0.725 against my 38/50 = 0.760, `p` = 0.81) was
> **retracted by its author** — one batch agreeing is not evidence, and a second batch
> gave 15/30 = 0.500. Across all its batches, 44/70 = 0.629 against my 38/50 = 0.760 is
> **consistent**, two-sample Fisher exact `p` = 0.164 (reproduced here), pooled
> 82/120 = 0.683 with a 95% interval containing both. An initial report of `p` = 0.019
> was a one-sample test treating my 38/50 as a known rate rather than an estimate.
>
> Two reasons not to read even the consistency as strong. The runs are **not
> exchangeable** — `calibrate_scale_on_folds` (`pulse.py:72-110`, `:381`) tunes the
> injected amplitude against each realisation's own noise, so the signal varies per draw
> and every nominal `p` here is anti-conservative. And my 38/50 is **end-to-end**
> recovery from the final periodogram while its predicate is evaluated at the last prune
> level, before the final ascend and resolve; the two need not agree even on identical
> data.
>
> **So: the shape is solid and the absolute levels are not.** "No losses before level
> 13" holds across two criteria, two injected parameter sets and three batches. Nothing
> below uses the absolute levels.
>
> **The decision is late, not early — so §4.3's title is not just unsupported, it is
> the wrong way round**, and the "structural reason the effect is small" is withdrawn.
> The stages where `pi` is 0.87–0.93 are inside the window where leaves are actually
> lost, not outside it. The consequences are favourable in every direction: `n` falls
> (below), and the `metric-gridding` gain measured in the real decision window 13–47 is
> 2.47x (Taylor) and 1.57x (Chebyshev) with `quadrature` closer in **54 of 54** cells,
> against the 1.71x / 1.15x its early window gave.
>
> **`DECISION_STAGES` is therefore misplaced**: [2, 6, 10, 14, 20, 28] sits almost
> entirely *before* the measured decision window, so §5.1's stratification and §9's
> strata are built over stages where nothing is being decided. Recomputing them over
> 13–47 has **not** been done — it would put a new headline number on an inherited,
> unverified N = 40 measurement from one arm and one parameter set, which is the exact
> failure mode this branch has repeated. It is the first thing to do if the
> configuration question is ever reopened.

---

## 5. The power calculation

Aggregating §4 over stages, weighting each by `reach(s) · phi(T(s) − mu(s)(1−L_A))`
(the probability of being in the race times the density of arm A's score at the
threshold), and taking `pi = E[w·D⁺] / E[w·|D|]`:

| ducy | `pi` | `p_disc` (modelled, M5 upper bound) | discordant pairs needed | **pairs needed** | core-hours at 50 s/pair |
|---|---|---|---|---|---|
| 0.05 | 0.787 | 0.699 | 22 | **31** | 0.4 |
| 0.10 | ~~0.697~~ **0.628** | 0.370 | 49 | ~~131~~ | 1.8 |
| 0.20 | 0.607 | 0.168 | 168 | **997** | 13.8 |

At the pilot's *measured* `p_disc = 0.30` (§6.3) rather than the modelled one, the same
three rows are n = 74 / 163 / 560, i.e. 1.0 / 2.3 / 7.8 core-hours.

> **The ducy 0.10 row of both tables above is superseded by §5.2**, which corrects
> `rho_AB`: `pi` = 0.628 and n = 391 at `p_disc` = 0.30, not 0.697 and 163 — and §9
> supersedes even that with 954, the count for the test actually run. The 0.05 and
> 0.20 rows have not been recomputed. They are left in place because the *shape* of the
> argument — that `p_disc` enters linearly and `pi` quadratically — is what §5 is for,
> and that is unchanged.

`m = (z_{α/2}·√¼ + z_β·√(pi(1−pi)))² / (pi − ½)²`, `n = m / p_disc`, 80% power,
two-sided α = 0.05.

### n across a range, since neither parameter is known to better than a factor of a few

| `pi` \ `p_disc` | 0.01 | 0.03 | 0.10 | 0.30 |
|---|---|---|---|---|
| 0.52 | 490 320 | 163 440 | 49 032 | 16 344 |
| 0.55 | 78 253 | 26 085 | 7 826 | 2 609 |
| 0.60 | 19 385 | 6 462 | 1 939 | 647 |
| 0.65 | 8 482 | 2 828 | 849 | 283 |
| 0.70 | 4 663 | 1 555 | 467 | 156 |
| 0.80 | 1 927 | 643 | 193 | 65 |
| 0.90 | 950 | 317 | 95 | 32 |

**Read this table as the answer, not the point estimate.** Two things fall out of it.

1. **`p_disc` barely matters, and it is the one parameter that is now measured.** It
   enters linearly, so moving it by a factor of 30 (0.30 → 0.01) moves `n` by 30 — from
   hours to a fortnight on one core. The pilot puts it at 0.08–0.30 (§6.3), the bottom
   two columns of the table, so the verdict does not turn on it. That is fortunate,
   because it is the parameter the model estimates worst (M5).

2. **`pi` is the whole question.** `n ∝ (pi − ½)^-2`, so `pi = 0.52` instead of 0.70
   costs a factor of 105. **`pi` is the single assumption the conclusion is most
   sensitive to**, and §4.1 says what drives it: the ratio
   `mu·(L_A − L_B) / √(2(1−rho_AB))`.

   The reassuring part is that this ratio is **more robust than it looks to the duty
   cycle**, which D55 identified as the dominant uncertainty in the amplitude numbers.
   Both the numerator and the denominator scale with the loss, so the ratio goes like
   `√L` rather than `L`: a factor of 15 in the losses between ducy 0.05 and 0.20 becomes
   a factor of 1.3 in `pi` (0.787 → 0.607) and a factor of 32 in `n` (31 → 997). All
   three remain affordable at 50 s/pair: 0.4 to 14 core-hours.

### 5.1 Stratification: measured, and it does not stratify (corrected in §5.2)

The brief asks for strata by within-cell position, measured rather than assumed.
`injection_power.py --strata 24` samples 24 signal positions and computes, for each, the
decision-weighted `pi` over stages 2–28:

    decision-weighted pi:  min 0.440   q25 0.595   median 0.645   q75 0.704   max 0.741
    positions with pi < 0.55 ("null"  stratum):   3/24
    positions with pi > 0.75 ("effect" stratum):  0/24

**There are no clean centre / face / corner strata in this configuration.** A signal is
placed once in 4-D parameter space; its per-stage geometry then follows and cannot be
dialled separately, and the resulting distribution is narrow and unimodal around 0.645.
A genuine null stratum exists but is rare — about 1 position in 8 — so it has to be
**pre-screened with `nearest_template_cheby.py`**, not obtained by random placement, and
a "corner" stratum with a large effect cannot be obtained at all. This is a measured
constraint on the design, and it is the opposite of what the sketch in the brief
assumed.

**Corrected for §5.2, and it changes the design, not just the numbers.**
`injection_power.strata` builds each position's `pi` from the same wrong `rho_AB`
(`injection_power.py:337`), so everything above inherits the error. Recomputed on the
same 24 positions with the corrected score correlation (`rho_check.py --strata 24`,
`rho_strata.json`):

| `rho_AB` from | min | q25 | median | q75 | max | `pi` < 0.55 | `pi` > 0.70 |
|---|---|---|---|---|---|---|---|
| profile overlap (in use) | 0.440 | 0.595 | 0.645 | 0.704 | 0.741 | **12%** | **29%** |
| **corrected** | 0.470 | 0.556 | **0.590** | 0.619 | **0.643** | **21%** | **0%** |

*Control:* the profile row reproduces the cached per-position values to 5 × 10⁻⁵ on all
24 positions, so the shift is the correction and nothing else.

Two consequences, in opposite directions:

- **The `effect` stratum of §9 does not exist.** It was defined as `pi > 0.70`; the
  corrected maximum over 24 positions is **0.643**. No position reaches it, so the
  three-stratum design collapses to two. The earlier `pi > 0.75` finding ("not
  reachable") was right for the wrong threshold — it is `pi > 0.70` that is unreachable.
- **The `null` stratum is roughly twice as easy to find**: about 1 position in 5 rather
  than 1 in 8. That is the control §9 leans on hardest, and pre-screening for it gets
  cheaper.

The distribution is also *narrower* (0.470–0.643 against 0.440–0.741), which is what a
larger common noise term does: it pulls every position toward ½.

---

### 5.2 `rho_AB` corrected: it was the profile overlap, not the score correlation

§10 row 3 and §11.4 both flagged this as unchecked and as the thing `pi` depends on most.
It is now checked (`rho_check.py`), and it was wrong, in the direction that costs pairs.

`injection_power.py` computes `rho_AB = snr_ratio(delta_A − delta_B, filt="matched")`,
which is `<p, p∗s>/<p, p>` — the overlap of two **signals**. M4 needs the correlation of
two **scores** under noise. Writing the score as a linear filter `h` on the folded
profile, and using that folding on two ephemerides assigns the same samples to bins
differing by `dPhi(t)`, both quantities are the *same* smearing factor `S(k)` under a
different weight:

    rho_score = Σ_k |H(k)|² S(k) / Σ_k |H(k)|²          the score correlation
    rho_prof  = Σ_k |P(k)|² S(k) / Σ_k |P(k)|²          what is in the code

They agree only if the filter **is** the pulse. It is not: at ducy 0.10 with `N_b = 64`
the bank selects a boxcar of **width 3**, whose `<k²>` is **4.10x** the pulse's, so it
weights exactly the high `k` where `S(k)` is small. The direction — `rho` overstated,
`pi` overstated, `n` understated — was predicted from the two spectra before anything
was simulated.

Measured over the same 83 exact cells the power calculation uses, at ducy 0.10:

| `rho_AB` from | median `1 − rho` | vs. profile | `pi` | `n` at `p_disc` = 0.30 | at 0.08 |
|---|---|---|---|---|---|
| profile overlap (in use) | 0.0088 | — | 0.697 | 161 | 603 |
| **with a signal present** | **0.0265** | **2.95x** | **0.628** | **391** | **1 464** |
| fixed selected boxcar | 0.0375 | 4.23x | 0.611 | 527 | 1 976 |
| max over bank, pure H0 | 0.176 | 19.3x | 0.560 | 1 805 | 6 768 |

**Use the signal-present row.** With a signal the two arms' bank maxima are anchored to
the same filter, which is the campaign's situation; the closed-form fixed-boxcar row
brackets it from below and agrees with it to within the re-optimisation the bank is free
to do. The pure-H0 row is only a **lower bound** on the correlation — with no signal the
two maxima land wherever they like — and it drives the model's `p_disc` to 1.0, which is
M5 breaking down, not a result. It is quoted to show where the bound is, not as a
candidate value.

*Controls.* The sampler draws the two folded profiles from their exact joint law in
Fourier rather than Monte-Carloing the time series, and is checked against the closed
form on a one-width bank before anything else runs. Substituting the *uncorrected*
values back through the same path reproduces the cached `pi = 0.69684` and
`p_disc = 0.36994` on all 83 cells exactly, so the correction is the only thing that
moved.

**What it does and does not change.** `pi` falls from 0.697 to **0.63**, and the pairs
needed at the measured `p_disc` rise by ~2.4x, to 391 at `p_disc` = 0.30 and 1 464 at
the pilot's worst 0.08. **These are aggregate-`pi` figures and are not what the campaign
needs** — §9 carries the operative count (954), because the primary test runs on the mid
stratum and the null stratum does not enter it. The pre-committed `n = 2 000` covers the
operating point either way, so the
verdict moves from comfortable to **adequate**, not from reachable to unreachable.

**Across duty cycle, and it is not uniform.** Only ducy 0.10 was recomputed in full, but
the leading-order factor — the `<k²>` ratio that sets the correction in the
small-smearing limit — is a pure spectral calculation and costs nothing:

| ducy | selected boxcar | `<k²>` boxcar | `<k²>` pulse | ratio |
|---|---|---|---|---|
| 0.05 | 3 | 104.2 | 101.6 | **1.03** |
| 0.10 | 3 | 104.2 | 25.4 | 4.10 |
| 0.20 | 9 | 37.4 | 7.0 | **5.35** |

**At ducy 0.05 the substitution is very nearly right, and for a comprehensible reason:**
a 5% pulse in 64 bins occupies ~3.2 bins, which is the width the bank selects, so the
filter *is* the pulse and the two weightings coincide. That is the one regime where the
original code was doing the right thing. At ducy 0.20 the correction is larger than at
0.10. So §5's table moves hardly at all in its top row and more than measured here in
its bottom one — the opposite of assuming "a similar factor" throughout.

These are leading-order indicators, not corrections: the measured inflation at ducy 0.10
was 2.95x against the 4.10 this predicts, because the linearisation overstates it once
the smearing is large and the bank can re-optimise. The 0.05 and 0.20 rows of §5's table
remain formally uncorrected.

## 6. Prerequisites and cost, measured

### 6.1 Pairing is possible — verified, and it needs no library change

The library **never seeds its noise**: `simulation/pulse.py:338` is a bare
`np.random.default_rng()`, and there is no `seed`/`rng` parameter on
`PulseSignalConfig`, `PulseSignalConfig.generate`, `DynamicProgramming`, `Pruning` or
`prune_dyp_tree`. `np.random.seed` does not help, because `default_rng` ignores the
legacy global seed. Worse, the injected amplitude itself varies run to run, because
`calibrate_scale_on_folds` (`pulse.py:392`) calibrates against the realised noise.

The route that works is to generate the `TimeSeries` once and persist `(ts_e, ts_v, dt)`.
Verified end to end in `injection_pilot.py`: the saved series round-trips to an identical
SHA-256, and two separate processes running the *same* arm on it produced **bit-identical
sorted score vectors** (`cc66644060fdc8f2` twice for `aggressive`, `9ccd6c6400eb642a`
twice for `quadrature`). Everything downstream of the time series is deterministic — a
grep of `src/pyloki` for `np.random|default_rng|shuffle|randint|permutation` hits only
`simulation/pulse.py`, `detection/thresholding.py`, `sensitivity/sim_ffa.py` and
`kepler.py`, none of which is on the search path.

One stochastic input sits *beside* the pipeline and must be frozen:
`DynamicThresholdScheme.__init__` (`thresholding.py:740`) is also unseeded, so a scheme
re-derived per run would differ between arms. **Generate each ladder once, commit the
array, reuse it** — which is what `schemes/cheby_*.npz` are for.

Storage is the one real constraint: a realisation is **32 MB** (2^22 samples × 2 arrays
× float32), so 5 000 pairs is 160 GB if they are all kept. The campaign must therefore
**stream** — generate one series, run both arms, delete — which is a driver change, not
a library one.

### 6.2 Per-pair cost

Numba JIT dominates the first run in a process (~25 s), so the number that matters is
marginal. Measured over 20–24 repetitions per arm:

| injected S/N | `max_sugg` | `aggressive` s/run | `quadrature` s/run | **s/pair** |
|---|---|---|---|---|
| 12 | 2^14 | 0.67 | 3.11 | 3.78 |
| 14 | 2^14 | 1.01 | 3.54 | 4.55 |
| 16 | 2^14 | 1.38 | 3.54 | 4.92 |
| 14 | **2^18** | 2.31 | **47.47** | **49.78** |

Plus ~0.15 s to generate each realisation. **The `max_sugg` row is the one to budget
from**: §6.4 shows the 2^14 runs are not measuring the thresholds they claim to. At
2^18, `quadrature` costs 13x what it did and the pair costs **≈ 50 s single-threaded**,
so 2 000 pairs is **28 core-hours** and 10 000 is 138. Still not an expensive
experiment, but two orders of magnitude from the naive figure, and it will rise again at
the `max_sugg` that finally stops binding.

### 6.3 What the pilot found

Five independent batches of 20–24 pairs with the per-strategy recalibrated Chebyshev
ladders. Every run is in `injection_pilot_results.json`.

| batch | S/N | `max_sugg` | `aggressive` rec. | `quadrature` rec. | disc. | quad-only | agg-only | max sat. |
|---|---|---|---|---|---|---|---|---|
| A | 12 | 2^14 | 3/24 | 3/24 | 2/24 | 1 | 1 | 0.880 |
| B | 14 | 2^14 | 7/20 | 4/20 | 5/20 | 1 | 4 | 0.978 |
| C | 16 | 2^14 | 16/20 | 17/20 | 3/20 | 2 | 1 | 0.986 |
| D | 14 | 2^18 | 13/20 | 15/20 | 6/20 | 4 | 2 | 0.979 |
| E | 14 | 2^14 | 11/20 | 11/20 | 4/20 | 2 | 2 | 0.972 |

*The ladders are honestly calibrated.* Batch A: both arms recovered exactly 3 of 24
(12.5%) against a nominal `P_d = 0.1031`. The equal-`P_d` recalibration does what it
claims, in both arms, on real data. That is the check all three retracted results lacked,
and it is the one assumption in §10 that has been independently verified.

*`p_disc` is real and measurable:* **0.08 to 0.30** across the five batches, against the
model's 0.37. M5 is high by ~2–4x, as expected. The campaign should be powered off the
measured value.

*The recovery criterion is not the weak link.* Over all 153 finite pilot runs the
mismatch to the truth is sharply bimodal: the largest among recovered candidates is
**0.60** and the smallest among non-recovered is **45.9**, a factor of 76 with
`m_recover = 1.0` inside the gap. So "recovered" means the search found the signal, not
that one of `quadrature`'s ~10⁵ extra candidates strayed close by luck — the obvious way
a candidate-count difference could fake a sensitivity difference, ruled out rather than
assumed away.

*Batches B and E say n = 20 is useless for anything but cost.* Same S/N, same
`max_sugg`, same ladders, different noise: `aggressive` 7/20 vs 11/20 and `quadrature`
**4/20 vs 11/20**. Part of that is ordinary binomial scatter and part is that
`calibrate_scale_on_folds` (`pulse.py:392`) calibrates the injected amplitude against
the realised noise, so the *signal* varies too. Either way, no 20-pair comparison on
this branch should be read as evidence of anything.

### 6.4 The `max_sugg` confound: what is established and what is not

**Established — the buffer changes the outcome on data the search has already seen.**

- *Single fixed realisation, all four combinations.* Changing only `max_sugg`:
  `aggressive` went from **0 surviving candidates at 2^14 to 110 at 2^18**;
  `quadrature` from 8 553 to 134 418. Same time series, same ladder, same everything.
- *Fixed set of 20 realisations, `aggressive` at both buffers.* 11/20 recovered at
  2^14, **13/20 at 2^18**, with 2 runs flipping, both from lost to recovered.
- *The ratchet fires at intermediate stages, not the last one.* In all four runs of the
  single-realisation test the minimum surviving score was **below** `thresholds[-1]`
  (2.87 vs 7.10 for `aggressive` at 2^18; 1.50 vs 7.00 for `quadrature`), so the final
  cut was the nominal one. The buffer is acting earlier in the tree, where nothing in
  the output reveals it.
- *It is silent.* `world_tree.py:527-551` sets the effective cut to
  `max(stage threshold, top-K, median)` on overflow; `PruneStats` records only the
  nominal `threshold` (`prune.py:685`). Worth reporting upstream on its own.

**NOT established — that the confound is arm-dependent.** The hypothesis is that
`quadrature` is penalised more, because its branching-pattern product is 2^28.4 times
`aggressive`'s (5.4e20 vs 1.5e12) so it overflows earlier and harder. The experiment
that would settle it is both arms at both buffers on **one** fixed set of realisations;
it was launched and **the `quadrature`-at-2^18 arm was killed before it finished**. The
partial result is in `injection_pilot_results.json` under `paired_max_sugg_test`.

**A retraction of my own, in this document, before it left the building.** An earlier
draft of §1 read the cross-batch pair B (2^14: `aggressive` 7-4, discordances 4:1 its
way) against D (2^18: `quadrature` 15-13, discordances 4:2 its way) as "the buffer sets
the direction of the answer". That comparison is between **independent batches of 20**
and cannot carry the claim: batch E, at B's exact settings, gives 11-11 with
discordances 2:2. The claim was removed. The reason it is recorded rather than quietly
deleted is that it is the same failure mode as the three retractions in `DECISIONS.md` —
a real mechanism, a measurement too small to see it, and a conclusion drawn anyway.

*And 2^18 is still not obviously enough.* Saturation reached 0.94–0.98 somewhere in
every batch. The buffer has to go higher, and `quadrature`'s cost goes with it.

---

### 6.5 The paired `max_sugg` experiment: §6.4 is NOT closed, and the confound may be bigger than the effect

§11.6's experiment, run as specified: **one** fixed set of **50** realisations at S/N 14,
all four `(arm, max_sugg)` cells, 2^14 and 2^18. The analysis was pre-registered in
`paired_max_sugg.py:plan()` and **committed before the first cell ran** (`d238142`);
`paired_max_sugg_results.json` carries the pre-registration alongside the result and all
200 runs.

**Recovery, all four cells (out of 50):**

| | 2^14 | 2^18 |
|---|---|---|
| `aggressive` | 30 | 38 |
| `quadrature` | 28 | 40 |

**1. Raising the buffer only ever helps, in both arms.** `aggressive` gained 8 and lost
**0** (McNemar p = 0.008); `quadrature` gained 12 and lost **0** (p < 0.001). Twenty
flips, all in the same direction. The ratchet is purely destructive to recovery, and
this confirms §6.4's central claim on 2.5x the realisations.

**2. The buffer pressure IS arm-dependent — measured, and significant.** At 2^18,
runs with `saturation > 0.9`:

| buffer | `aggressive` | `quadrature` | paired `n01`/`n10` | McNemar `p` |
|---|---|---|---|---|
| 2^14 | 5/50 | 8/50 | 8 / 5 | 0.58 |
| **2^18** | **1/50** | **10/50** | **10 / 1** | **0.012** |

with median saturation 0.165 vs **0.724** and median candidate counts 43 314 vs
**189 828**. At 2^14 both arms are slammed against the buffer and the asymmetry is
invisible; at 2^18 `aggressive` has largely escaped and `quadrature` has not. **The
mechanism for an arm-dependent bias is present and is now measured**, which is more than
§6.4 could say.

**3. Whether that becomes an arm-dependent OUTCOME bias is NOT resolved — and the point
estimate is worse than the effect being measured.** The pre-registered primary, a sign
test on `Delta_i = d_i(2^18) − d_i(2^14)`:

    nonzero Delta   10/50        k positive   7        exact two-sided p = 0.34
    mean Delta      +0.080       95% CI  [−0.043, +0.203]
    campaign's effect of interest                       0.06

**This is a failure to reject, not a null.** Read the interval, not the p-value: the
best estimate of the arm-dependent buffer bias is **+0.080, which is larger than the
0.06 the campaign is powered to detect**, and it points *toward* `quadrature` — the same
direction as the effect the campaign would be trying to claim. The interval is
consistent with no interaction, and equally consistent with one three times the size of
the signal. n = 50 simply cannot tell them apart.

**This was pre-committed, in exactly these words, before the data existed**: *"NOT a
green light, and must not be reported as one … it excludes a gross arm-dependence and
nothing more."* That is the reading, and the point estimate makes it stronger than
anticipated, not weaker.

**4. What it would take to close it.** From the measured `se`, bounding the interaction
below 0.06 needs **211 realisations per cell**, 4.2x this experiment — roughly 10–14
core-hours for all four cells at the measured per-pair costs. That is the honest price
of an experiment that could *clear* the campaign rather than merely fail to condemn it.

**5. The discordance direction flips with the buffer — and this is exactly the claim
§6.4 retracted, so it is not being made again.** `n01`/`n10` is 4/6 at 2^14 and 5/3 at
2^18. On *one fixed set* this time rather than independent batches, so it is not the
same error — but 10 and 8 discordances cannot carry a direction, the sign test on
`Delta` **is** the correct test of precisely this flip, and it gives p = 0.34. Recorded
as an observation, tested, and not concluded from.

**6. S/N 14 is the wrong operating point at 2^18 — and in the opposite direction to
what the `metric-gridding` model predicts.** §9 wants on-grid `P_d ≈ 0.5` in both arms,
because discordance is maximised at the steepest part of the recovery curve. Measured
here: **0.76 and 0.80**. The operating point is well *above* the optimum, so the 3-point
pilot should bracket **below** S/N 14 — while D78 from the `metric-gridding` session
recommends 15–17, which is further above it. The two disagree on direction; this one is
a direct measurement of recovery at the final buffer and that is the criterion §9 states.

#### Verdict on §6.4

**Not closed, and the campaign must not run at 2^18.** §7.1's remedy is no longer
advisory: `quadrature`'s 99th-percentile saturation at 2^18 is **0.975**, against the
< 0.9 the design requires, and it is 20% of runs above 0.9 versus `aggressive`'s 2%. The
buffer must go higher — and the measured asymmetry says the *reason* to raise it is
specifically that `quadrature` is still inside it when `aggressive` has left.

### 6.6 Stage 1 of the remedy fails: `quadrature`'s candidate set never converges

§7.1 requires the buffer to stop binding in **both** arms. Swept over four buffers on
**one fixed set of 10 realisations** (the first 10 of §6.5's, so every row below is the
same data):

| arm | `max_sugg` | median `ncand` | median saturation | growth exponent |
|---|---|---|---|---|
| `aggressive` | 2^14 | 1 446 | 0.088 | — |
| `aggressive` | 2^18 | 16 778.5 | 0.064 | 0.88 |
| `aggressive` | 2^19 | **16 778.5** | 0.032 | **0.00** |
| `aggressive` | 2^20 | **16 778.5** | 0.016 | **0.00** |
| `aggressive` | 2^21 | **16 778.5** | 0.008 | **0.00** |
| `quadrature` | 2^14 | 10 906 | 0.666 | — |
| `quadrature` | 2^18 | 152 007 | 0.580 | 0.95 |
| `quadrature` | 2^19 | 390 575 | 0.745 | 1.36 |
| `quadrature` | 2^20 | 808 116 | 0.771 | 1.05 |
| `quadrature` | 2^21 | **1 610 052** | 0.768 | **0.99** |

`growth exponent` = `d log2(median ncand) / d log2(max_sugg)` against the previous
buffer. Zero means the candidate set has converged and the buffer is irrelevant; one
means the buffer is the only thing setting the answer.

**`aggressive` converges at 2^19 and stays converged.** The claim is stronger than the
medians: the **per-realisation** candidate counts are *bit-identical* across 2^19, 2^20
and 2^21 — all ten of them, [9, 44, 628, 1482, 14 921, 18 636, 32 949, 33 184, 54 011,
198 037] — so the buffer is provably irrelevant over a factor of 4, not merely
statistically indistinguishable. Median 16 778.5, saturation falling 0.032 → 0.016 →
0.008. That is what a non-binding buffer looks like.

*It is 2^19, not 2^18, and the median hides why.* The medians coincide at 2^18, but the
per-realisation counts do not: `tim_0004` returns 192 661 at 2^18 (saturation 0.735) and
198 037 at every larger buffer — still clipped at 2^18. Over all 50 of §6.5's
realisations one more (`tim_0045`) exceeds 0.9 saturation at 2^18, giving p99 = 0.915,
which fails §7's criterion. **So `aggressive`'s clearing point is 2^19**, and a
median-only reading of this table would have put it one buffer too low. This does not
affect the conclusion below, which turns entirely on `quadrature`.

**`quadrature` never converges.** Its exponent is 0.95, 1.36, 1.05 across three
successive increases, and its median saturation does not fall — it *rises*, 0.580 →
0.745 → 0.771. Doubling the buffer simply doubles the candidate count. The mechanism is
the ratchet itself: the effective cut is `max(scheme threshold, top-K, median)`, so as
the buffer grows the cut relaxes to keep it full, and the buffer stays full at every
size.

#### What this means

**There is no feasible `max_sugg` at which `quadrature` clears §7.1's criterion, so
stage 1 of the remedy cannot be completed and the campaign in §9 cannot be run as
designed.** This is not a cost problem that a bigger machine solves. The count is
tracking the buffer with exponent ≈ 1 at 2^20 with no sign of a knee; the true
unconstrained count is far above 10⁶, consistent with `quadrature`'s branching product
being 2^28.4 times `aggressive`'s (5.4 × 10²⁰ vs 1.5 × 10¹²). Cost at 2^20 is already
152 s/run against `aggressive`'s 1.3.

**The consequence for everything measured in the `quadrature` arm on this branch is
that the scheme's thresholds were never the operative cut.** At every buffer ever used
here — 2^14, 2^18, 2^19, 2^20 — `quadrature` ran against the ratchet, not against
`cheby_quadrature.npz`. The equal-`P_d` recalibration (§6.3) is therefore not doing in
the `quadrature` arm what §6.3 says it does, and the one check that appeared to validate
it — batch A, 3/24 in both arms at S/N 12 — was run at 2^14, where **both** arms were
saturating and the agreement cannot be attributed to the ladders.

**What is still answerable.** §7.1 already separated two questions, and only one of them
dies here:

- *Does tiling matter for the idealised thresholded search?* Needs a non-binding buffer
  in both arms. **Not answerable in this configuration.**
- *Does tiling matter for what a user actually runs?* Runs at the shipped default and
  treats the ratchet as part of the system under test. **Still answerable**, and §6.5's
  apparatus answers it directly — but it is a statement about pyloki-as-shipped, not
  about tiling, and it must be reported as such.

A third option is to shrink the configuration — fewer segments, smaller
`prune_poly_order`, smaller `branch_max` — until `quadrature`'s true candidate count is
affordable, and ask the idealised question there. That is a different experiment with a
different external validity, and it has not been costed.

**A fourth option, which nothing above has considered and which keeps both the
configuration and the idealised question: make the ladder stricter.** The candidate
count is set by how much the threshold scheme admits, and the schemes here are
calibrated to `P_d` = 0.1031. A ladder calibrated to a smaller `P_d` has higher
thresholds at every stage and admits proportionally fewer leaves, so `quadrature`'s true
count could fall by orders of magnitude while the comparison stays a fair equal-`P_d`
one — both arms simply move to a stricter operating point together. That changes *which*
operating point the campaign speaks about, which is a real cost and must be stated, but
it does not change the configuration or make the question hypothetical.

It is also cheap to test before committing to anything: generate one stricter
`quadrature` ladder, run 10 realisations at 2^20, and see whether the count converges.
If it does, stage 1 is back. **Untested** — noted here because it is the only route
found so far that preserves both the deployed configuration and the idealised question,
not because there is evidence it works. The obvious risk is that a stricter ladder
lowers `P_d` below the ~0.5 operating point §9 wants, and the two constraints may not be
simultaneously satisfiable; that is exactly what the 10-realisation test would show.

#### The 2^21 confirmation: no knee

The conclusion was held provisional on one more buffer, because it kills the campaign
and consecutive exponents are noisy. 2^21 is now in and it confirms rather than softens
it. Span exponents, which are less noisy than step-to-step ones:

| span | median `ncand` | exponent |
|---|---|---|
| 2^18 → 2^19 | 152 007 → 390 575 | 1.361 |
| 2^18 → 2^20 | 152 007 → 808 116 | 1.205 |
| **2^18 → 2^21** | **152 007 → 1 610 052** | **1.135** |
| 2^19 → 2^21 | 390 575 → 1 610 052 | 1.022 |
| 2^20 → 2^21 | 808 116 → 1 610 052 | 0.994 |

Over a factor of **8** in buffer the candidate count grows by a factor of **10.6**, and
median saturation is pinned at 0.745 / 0.771 / 0.768 — flat to three buffers. Meanwhile
`aggressive` returns **16 778 candidates at 2^18, 2^19, 2^20 and 2^21**, identical to
the unit across the same factor of 8, with saturation falling 0.064 → 0.008. The two
behaviours could not be more cleanly separated on the same data. There is no knee, and
cost is now 302 s/run against `aggressive`'s 1.4.

*Scope.* Measured at Phase 3's configuration only (268.4 s, 64 segments,
`poly_order = 4`, `N_b = 64`, `branch_max = 16`, Chebyshev, `eta = 1`) and on 10
realisations. The `aggressive` convergence is exact. The `quadrature` non-convergence
now rests on four buffers spanning a factor of 8, with the two widest spans giving
exponents 1.022 and 0.994. **No longer provisional.**

### 6.7 The stricter ladder fails too, and the ratchet is not set by the ladder at all

§6.6 left four routes open. Option 4 — recalibrate the ladder to a smaller `P_d` so it
admits fewer leaves, keeping both the deployed configuration and the idealised question
— is the only one that preserved everything. **It does not work**, and the reason rules
out any variant of it.

**A better criterion than saturation, and a correction to §6.5 and §6.6.** Both sections
used `saturation = ncand / max_sugg` as the evidence for ratcheting. It is only a
correlate: saturation is measured at the *end* of a level, while the ratchet fires on
overflow *during* one and never relaxes, so a run can finish well below the buffer
having spent most of the level ratcheted. The `threshold_eff` field added by
[PR #14](https://github.com/pravirkr/pyloki/pull/14) measures it directly — the cut
actually applied, per level — so the question becomes exactly the one that matters:
**did the search run the ladder it was given?**

**Measured at `max_sugg` = 2^20 on three fixed realisations, 64 levels each:**

| arm | ladder | realised `P_d` | thresholds | levels ratcheted | clean runs | recovered |
|---|---|---|---|---|---|---|
| `aggressive` | shipped | 0.1031 | 2.10–7.70 | **0 / 192** | **3/3** | 3/3 |
| `quadrature` | shipped | 0.1031 | 2.00–7.90 | 60 / 192 (20.0/run) | **0/3** | 3/3 |
| `quadrature` | matched baseline | 0.1083 | 1.40–8.10 | 60 / 192 (20.0/run) | **0/3** | 3/3 |
| `quadrature` | **strict** | **0.0100** | 2.40–8.50 | 56 / 192 (18.7/run) | **0/3** | 3/3 |

`aggressive` is the control and it is perfectly clean — 0 of 192 levels, three runs out
of three running exactly the ladder they were given. That demonstrates what §6.6's
convergence only predicted, and it means the `quadrature` numbers are not an artefact of
the instrument.

**A ten-fold reduction in the detection target buys 4 levels out of 192.** The strict
and matched-baseline ladders come from the *same* `DynamicThresholdScheme` run,
backtracked at different `P_d`, so they differ only in the target and not in the
unseeded draw (§6.1). Going from `P_d` = 0.108 to 0.010 — thresholds raised by ~1.0 at
the low end — moves the ratcheted count from 20.0 to 18.7 per run and leaves **zero**
clean runs. The candidate count does not fall either: median 730 462 against 949 921,
i.e. it went *up*, because the count is pinned to the buffer rather than set by the
scheme.

**Why no ladder can fix this.** The effective cut is
`max(nominal, top-K, median)` (`utils/world_tree.py:547`). Measured max excess over the
nominal is **4.19 / 4.98 / 5.32** on thresholds of order 2–8, so the operative cut is
roughly **double** what the scheme asked for, and it is `top-K` that is setting it —
a function of the buffer size and the score distribution, **not of the ladder**. Raising
the nominal threshold by 1.0 changes nothing when top-K already exceeds it by 5.

**Measured directly, per level, on one realisation under both ladders.** Recording the
nominal and effective cut at every level (`ratchet_probe.py` stores `levels_detail`):

| level | `thresh` base | `thresh` strict | Δ nominal | `eff` base | `eff` strict | **Δ effective** |
|---|---|---|---|---|---|---|
| 4 | 1.80 | 3.30 | +1.50 | 2.93 | 3.30 | +0.37 |
| 5 | 2.30 | 3.50 | +1.20 | 4.35 | 4.40 | +0.05 |
| 7 | 2.90 | 4.20 | +1.30 | 5.39 | 5.38 | −0.01 |
| 10 | 3.90 | 5.00 | +1.10 | 5.92 | 5.92 | **0.00** |
| 13 | 4.60 | 5.80 | +1.20 | 6.36 | 6.37 | +0.01 |
| 19 | 6.00 | 6.70 | +0.70 | 7.29 | 7.31 | +0.02 |

Over all 20 ratcheted levels:

    mean change in NOMINAL threshold (strict - baseline):  +0.870   sd 0.344
    mean change in EFFECTIVE cut                        :  -0.063   sd 0.175
    |d eff| / |d thresh| = 0.07

**Raising the ladder by 0.87 moved the cut that actually ran by −0.06.** The effective
cut is independent of the scheme to within 7%, and at level 10 it is identical to two
decimal places while the nominal rose by 1.10. This is not a weak effect that a stricter
ladder might overcome; the ladder is simply not the operative term.

**Level 4 is the exception, and it explains the one level the strict ladder gained.**
There the strict nominal (3.30) rises *above* the baseline's effective cut (2.93), so
the nominal becomes binding again and `eff = thresh` exactly. That is `max(nominal,
top-K, median)` behaving as written, and it is why the ratcheted count went 20 → 19: the
single level recovered is precisely the one where the nominal crossed above top-K. To
recover all of them the ladder would need nominal > top-K at *every* level — thresholds
of 6–8 from the first stage — which is not a stricter operating point but no search at
all.

#### Verdict on §6.6's options

**Option 4 is dead, structurally rather than by tuning.** No `P_d` target fixes a cut
that is not being set by `P_d`. What remains:

- **Option 2 — the user-facing question** at the shipped default, ratchet included.
  Answerable now. A claim about pyloki-as-shipped, not about tiling.
- **Option 3 — a smaller configuration.** Untested, different external validity.
- And a route none of the four anticipated: **change the buffer policy rather than the
  ladder.** The ratchet exists because `max_sugg` is a hard cap enforced by discarding
  the lowest-scoring candidates. A policy that instead refused to proceed — or that
  reported the shortfall — would let `quadrature` fail honestly rather than silently
  substitute a different cut. That is a library change, not a campaign design, and it is
  out of scope here; it is noted because it is the only remaining way to ask the
  idealised question at the deployed configuration.

*Scope.* Three realisations per cell at one buffer, one arm's control. The effect is not
marginal — 0/192 against 60/192 — so n = 3 is adequate to separate them, but the
per-run ratcheted count (18.7 vs 20.0) is not resolved at this n and nothing is claimed
from that difference beyond "not zero".

---

## 7. What has to change before this is worth running — **superseded by §6.6**

*(References elsewhere of the form "§7.1" mean item 1 of the numbered list below; §7 has
no numbered subsections.)*

> **§6.6 supersedes item 1 below: the ratchet cannot be eliminated for `quadrature`.**
> Item 1 was written assuming a large enough buffer exists. It does not, in this
> configuration. The two-questions paragraph inside it survives and is now the whole of
> the decision; the "start the search at 2^20" instruction was carried out (§6.6) and
> answered in the negative.

1. **Characterise the buffer ratchet, then eliminate it.** This is the whole of the
   problem, and §6.4 is where it stands: the buffer changes outcomes on fixed data, and
   whether it does so asymmetrically between the arms is the open question. Raise `max_sugg` until the 99th percentile of `saturation = ncand /
   max_sugg` is below 0.9 in both arms over a 50-pair pilot, and record it per run as a
   first-class output (`injection_pilot.py` already does). 2^14 is badly wrong and 2^18
   still saturates somewhere in every batch; start the search at 2^20. Any run above 0.9 is not measuring the
   scheme's thresholds — exclude it and report the exclusion count. Re-measure the
   per-pair cost at the final value: it went from 4.6 s to 49.8 s between 2^14 and 2^18
   and will rise further.

   Worth separating two questions the buffer conflates. *Does tiling matter for the
   idealised thresholded search?* needs a non-binding buffer, and is the question this
   branch has been asking. *Does tiling matter for what a user actually runs?* would run
   at the shipped default and treat the ratchet as part of the system — a legitimate
   question, but a different one, and it would have to be reported as such rather than
   as a statement about tiling.
2. **Pre-screen the null stratum** with `nearest_template_cheby.py` (§5.1) rather than
   assuming random positions give one. About **1 in 5** does, after the §5.2 correction
   (it was 1 in 8 before it). Screen against the *corrected* `pi`, not the cached one —
   `rho_check.py --strata` is the version that does.
3. **Stream the realisations** instead of persisting them (§6.1), or budget 160 GB.
4. Keep the ladders frozen and committed. Do not re-derive `DynamicThresholdScheme` per
   run; it is unseeded.

---

## 8. What can be concluded without running it

> ⚠ **Three of the original six bullets are withdrawn.** This is the section most likely
> to be lifted into a summary, so what is gone is kept visible rather than deleted.

**Withdrawn:**

- ~~The equal-`P_d` recalibration works on real data: both arms recovered 3/24 at S/N 12.~~
  **Batch A ran at `max_sugg` = 2^14, where both arms saturate (§6.6), so the realised
  cut in each was the ratchet and not the ladder.** The agreement is real but cannot be
  credited to the recalibration. §10 row 8.
- ~~The tiling geometry is largest where nothing is being decided.~~ **Backwards.** The
  measured per-stage survival profile (§4.3 ⚑) puts the first loss at level 13 and the
  median at 27, with **0%** of losses by stage 10 where the modelled curve put 99%. The
  stages where `pi` is 0.87–0.93 are *inside* the decision window, not outside it. This
  was the structural reason the whole effort expected a small effect, and it does not
  hold.
- ~~A defensible upstream sentence, without any campaign: at equal `P_d` the tiling
  choice moves the decision by ~0.04 against a 0.13 stochastic difference…~~ **Rests on
  both of the above plus a superseded number** — the stochastic term is 0.230, not 0.133
  (§5.2), so it is ~5x the deterministic term rather than 3x, and "the stages where the
  decision is made" are the wrong stages. **Nothing in this document currently supports
  an upstream sentence about tiling.**

**Stands:**

- **`max_sugg` changes search outcomes silently.** On one fixed time series, `aggressive`
  returns 0 candidates at 2^14 and 110 at 2^18 (§6.4); over 50 fixed realisations,
  raising the buffer produced 20 recovery flips, all toward recovery and none away
  (§6.5). The overflow ratchet at `utils/world_tree.py:547` raises the effective cut and
  `PruneStats` records only the nominal threshold (`prune.py:685`). **Patched and
  verified** — `06_upstream_max_sugg_logging.md`.
- **For a high-branching configuration the ratchet is permanent, not occasional.**
  `aggressive` converges at 16 778 candidates and stays there across 2^18–2^21;
  `quadrature`'s count is proportional to the buffer over the same factor of 8 (span
  exponents 1.135 / 1.022 / 0.994) with saturation pinned near 0.77 (§6.6). So a
  high-branching user runs against the ratchet at every buffer they can afford, with the
  log asserting their scheme is applied. **This is the strongest result on the branch
  and it needs no campaign.**
- **Pairing on the noise realisation does not make the comparison deterministic.** The
  template-to-template noise difference is sd **0.230** at ducy 0.10 (§5.2), ~5x the
  amplitude difference it is meant to reveal. Any future design that treats a paired run
  as a clean A/B should be checked against this number first.
- **The discordance rate is 8–30% at plausible operating points**, so a paired design is
  the right one and an unpaired one would be hopeless. Measured, and unaffected by the
  above.
- **The per-stage `pi` profile itself stands** — 0.515–0.558 early, 0.87–0.93 late
  (§4.3). It is geometry plus §5.2's corrected `rho` and contains no survival model. It
  is the *weighting* over those stages that was wrong, not the profile.

---

## 9. The design, pre-committed

**Arms.** `poly_basis="chebyshev"`, `use_moving_grid=True`, `prune_poly_order=4`,
`branch_max=16`; `tiling_strategy` ∈ {`aggressive`, `quadrature`} and nothing else
differs. Config from `phase3_config.make_config`.

**Pairing.** One `TimeSeries` per pair, generated once, persisted as `(ts_e, ts_v, dt)`,
consumed by both arms, deleted after both have run. One arm per process. Verified
bit-deterministic (§6.1).

**Thresholds.** `schemes/cheby_aggressive.npz` and `schemes/cheby_quadrature.npz`,
committed, Viterbi at `P_d = 0.1`, never re-derived.

**Operating point.** Injected S/N chosen so on-grid `P_d ≈ 0.5` in *both* arms, since
discordance is maximised at the steepest part of the recovery curve. At `max_sugg =
2^18` that is **below** S/N 14 — §6.5 measures `P_d` = 0.76 and 0.80 there on 50
realisations — and the curve moves with the buffer (the same S/N gave 7/20 and 4/20 at
2^14), so fix it by a 3-point pilot **at the final `max_sugg`** and then do not change
it.

*Why this criterion is stated in `P_d` and not in S/N, which turns out to matter.* The
`metric-gridding` session initially recommended S/N 15–17 from its own discordance-peak
model, the opposite direction from the measurement above. It withdrew that (D78) on
finding that its model's discordance peaks at model `P_d` = 0.567 — i.e. at this
criterion. The two criteria were the same criterion in different variables, and the
disagreement was entirely in its `S/N ↔ P_d` mapping, which is ~2x too shallow: it
misplaces *which S/N* gives `P_d ≈ 0.5` while getting the `P_d` at the optimum right.
**A rule parameterised in the observable is immune to an error in the mapping from the
input; one parameterised in S/N is not.** That is the general reason to state it this
way, and it is worth more than the specific disagreement it resolved. (Reported by that
session; not re-derived here.)

**`max_sugg`.** Raised until the 99th percentile of `ncand / max_sugg` over a 50-pair
pilot is below 0.9 in both arms. Recorded per run. Runs above 0.9 are excluded and the
exclusion count is reported.

**Strata.** Three, assigned by pre-screening candidate positions with
`nearest_template_cheby.py` over the decision stages (2–28) and binning on the
decision-weighted `pi` of §5.1:
- **null** (`pi < 0.55`, ~1 position in **5** after the §5.2 correction): must show no
  strategy difference.
- **mid** (`pi ≥ 0.55`, the remaining ~4 in 5).
- ~~**effect** (`pi > 0.70`)~~ — **deleted. It is not reachable.** The corrected
  per-position maximum over 24 sampled positions is `pi` = 0.643 (§5.1), so no position
  qualifies. The design is **two strata**, and the secondary "monotonicity of `pi`
  across the three strata" below is correspondingly a two-point comparison.

**Controls.**
- *Built-in null (the important one).* The null stratum must give `n01 ≈ n10`. If it
  does not, the run is measuring threshold-calibration mismatch or the buffer ratchet,
  not coverage, and **the primary result is void regardless of what it says**. This is
  pre-committed: the null is checked and reported before the primary test is looked at.
- *Positive control.* Signal exactly on a base-grid point, `L_A = L_B = 0` at stage 0.
  Both arms must agree; a discordance here is a pipeline defect.
- *Negative control.* Signal displaced beyond any claimed cell. Both arms must fail; if
  either recovers, the recovery criterion (`m_recover = 1.0` in the full-baseline
  metric) is too loose.

**n and stopping rule.** Fixed n, no interim looks: **n = 2 000 pairs**. That needs
`m = 194` discordant pairs for `pi = 0.60`, `85` for 0.65 and `47` for 0.70, so it gives
≥80% power for `pi ≥ 0.65` across the whole measured `p_disc` range (0.08–0.30), and for
`pi = 0.60` provided `p_disc ≥ 0.10` — which it was, at 0.30, at the operating point the
design picks. Covering `pi = 0.60` at the pilot's worst `p_disc = 0.08` would take 2 425.
At the measured 49.8 s/pair, 2 000 pairs is **28 core-hours**; recompute both from §5
once `max_sugg` is settled, since the cost rises with it. **No optional
stopping** — this branch has retracted three results, and a stopping rule that looks at
the data is how a fourth would happen.

**What the corrections do to `n`, computed rather than asserted.** Two of the
corrections above stack onto the primary test: `pi` is now the **mid stratum's**
(median 0.593 over the 19 of 24 positions at `pi` ≥ 0.55), not the all-cell aggregate
0.628, and the null stratum's 21% of screened positions do not enter the primary test,
so the pair count has to be inflated by 1/0.79.

| `p_disc` | mid-stratum pairs | **total pairs** | core-hours at 50 s/pair |
|---|---|---|---|
| **0.30** (measured at the operating point) | 755 | **954** | 13.2 |
| 0.20 | 1 132 | 1 430 | 19.9 |
| 0.15 | 1 509 | 1 906 | 26.5 |
| 0.10 | 2 263 | 2 859 | 39.7 |
| 0.08 (pilot's worst, batch A at S/N 12) | 2 828 | 3 572 | 49.6 |

**And the stage weighting behind all of these is itself an unvalidated model.** `power()`
weights each stage by `reach(s)` = the ladder's own cumulative H1 survival `succ_h1`.
That is a modelled quantity from the same single-leaf Gaussian family that the
`metric-gridding` session has now measured to be badly wrong in the other direction (its
survival model predicts `P_d` = 0.368 at S/N 14 where §6.5 measures **0.760** in the
converged `aggressive` arm, a 4.5-S/N error). Re-running the aggregation under different
stage weightings, changing nothing else:

| stage weighting | `pi` | pairs at `p_disc` = 0.30 | total (÷0.79) |
|---|---|---|---|
| ladder `succ_h1` (what §5.2 and the table above use) | 0.628 | 391 | **954** |
| flat — every stage weighted equally | 0.654 | 268 | 339 |
| `sqrt(succ_h1)` — decay halved in log | 0.642 | 319 | 404 |
| late stages only (> 10) | 0.680 | 194 | 246 |
| **early stages only (≤ 10)** | **0.554** | 2 245 | **2 842** |
| **measured profile** (§4.3 ⚑, inherited) | **0.648** | **291** | **368** |

**So `n` spans roughly 340 to 2 840 — a factor of 8 — on an assumption that has never
been checked, which is larger than the `rho` correction of §5.2.** The direction the
falsification implies is the *favourable* one: if the real search survives far better
than modelled, the weights should decay less steeply than `succ_h1`, which moves `pi`
up and `n` down toward 340. But §4.3's claim that the decision is made early, where
`pi` is only 0.554, is the pessimistic corner and rests on the same modelled profile.
Both cannot be read off a model that is wrong by 4.5 S/N.

**`n = 954` should therefore be quoted as a point in that range, not as the answer** —
and it is a point under the weighting now known to be backwards. Under the *measured*
profile the aggregate gives **368 total pairs**, and the expensive 2 842 corner is
excluded outright, because it assumed the early decision the measurement rules out.

Three cautions against promoting 368 in 954's place, the third measured after the fact.

1. It rests on an inherited measurement from one arm; hardened since (two injections,
   N = 30 each, plus the original N = 40) but still `aggressive`-only.
2. It is the *aggregate* `pi`. The operative figure is the mid-stratum one, which cannot
   be recomputed without also moving `DECISION_STAGES` into the measured window (§4.3 ⚑)
   and re-running the stratification.
3. **It is not stable against the profile's own revision.** Re-deriving the profile under
   an independent criterion (metric mismatch `m ≤ 1.0` rather than a calibrated
   excursion) moved the median loss level from 27 to **15** — inside the range the
   interpolation smooths over. Recomputing `pi` across profiles consistent with that
   revision:

   | reach profile | `pi` | total pairs |
   |---|---|---|
   | excursion criterion, median loss 27 (gives 368) | 0.648 | **368** |
   | metric criterion, gentle, median 15 | 0.633 | 458 |
   | metric criterion, as reported | 0.624 | 533 |
   | metric criterion, steep | 0.609 | **691** |

   **So `n` is 368–691, not 368** — a factor of 1.9 that turns entirely on where the
   survival curve falls, which the two criteria do not agree on.

What survives all three cautions, and is the only thing that should be quoted: **the
early-decision corner (2 842 pairs) is excluded, and `n` is in the high hundreds rather
than the low thousands.** The precise value is not established.

**The campaign remains affordable — 13 to 50 core-hours — but `n = 2 000` no longer
covers the whole measured `p_disc` range.** It covers `p_disc` ≥ **0.145**. Two honest
readings, and they differ:

- *The design's own operating point is chosen to maximise discordance* (on-grid
  `P_d ≈ 0.5` in both arms), and there the pilot measured `p_disc` = 0.30, where 954
  pairs suffice and `n = 2 000` has real headroom. This is the relevant number.
- *The 0.08 corner is batch A at S/N 12*, where both arms recovered 3/24 — far from the
  operating point, and a regime the design explicitly avoids. Quoting it as the
  requirement would be powering against a configuration the campaign will not run in.

So `n = 2 000` stands, on the condition §9 already imposes: the operating point is fixed
by a 3-point pilot at the final `max_sugg` and `p_disc` is re-measured there. If that
pilot returns `p_disc` < 0.145, `n` must be raised to ~3 600 (50 core-hours) or the
design re-scoped — and that check is now a precondition, not a formality.

**Pre-committed analysis.** Primary: two-sided McNemar (exact binomial on the discordant
pairs) of `pi = 1/2`, over the **mid** stratum (the effect stratum is deleted above;
what was "pooled over mid and effect" is now just mid), α = 0.05. Secondary,
descriptive only, not tested: per-stratum `pi`, and its ordering across the **two**
strata. Reported unconditionally: the null-stratum result, the saturation
distribution, the exclusion count, and `n01` and `n10` separately. Effect size of
interest: `pi = 0.60`, i.e. 3 discordances for `quadrature` to every 2 for `aggressive`.

---

## 10. Assumption register, ordered by how much the answer moves

| # | assumption | if wrong, `n` moves by | note |
|---|---|---|---|
| 1 | **`pi ≈ 0.63`** (§5.2; was 0.70), i.e. `mu·(L_A−L_B) / √(2(1−rho_AB))` | **up to 10⁵x** — `n ∝ (pi−½)^-2`, and `pi → ½` makes it unbounded | **still the one to worry about**, and now closer to ½ than it was: row 3 moved it once already. |
| 2 | `L_A` and `L_B` are the losses of the leaf that actually scores | factor of a few, **sign not guaranteed** | Both are *upper* bounds (M3, D52 convention): the search minimises sup-norm, not loss, over 30 retained leaves. If `aggressive`'s true best leaf is better than retained, `pi → ½`; if `quadrature`'s is, `pi` rises. **Partly retired, on trust:** the `metric-gridding` session reports (its `66064c3`) that over stages 1–10 both arms' losses are *exact* — identical at 4M and 16M nodes — so no bound asymmetry can enter the window where ~99% of decisions are made; it stands for the late stages. **Not re-run here**, and §11.1's rule applies: this is an inherited number, not a verified one. It also reports a genuine sign reversal there (`quadrature` worse at stages 2 and 8), which `power()` handles correctly — `_folded_normal_moments` uses the signed mean, so a reversed cell drags `pi` down rather than violating anything. |
| 3 | ~~`rho_AB` from the between-template phase residual is the score correlation~~ **MEASURED, and it was wrong** | **2.4x in `n`, already applied** | §5.2. It is the *profile* overlap; the boxcar bank does not preserve it (`<k²>` ratio 4.10). Median `1−rho` 0.0088 → 0.0265, `pi` 0.697 → 0.628. No longer an assumption at ducy 0.10; still one at 0.05 and 0.20. |
| 4 | M6 — survival is decided by one covering leaf, not a max over many | raises `pi`, lowers `n` | ~~Deliberately conservative for the stated question.~~ **Now measured to be the DOMINANT error, not a margin.** A single-leaf model is pessimistic by ~4.5 S/N against the converged `aggressive` arm (§4.3 box); the `metric-gridding` session recovers most of the gap with `K` partially decorrelated near-covering leaves (`K` = 100, `rho` = 0.90 → 0.703 against a measured 0.760). `aggressive` had 16 778 candidates alive at 2^18, so the supply exists. Still the channel an arm-dependent `max_sugg` bias would ride on (§6.4), and that channel is now known to be wide. |
| 5 | M2 — `mu(s) = snr_final·√((s+1)/nseg)` | ~20% in `pi` | Same model the shipped ladder uses, so an error here is an error in the ladder too. |
| 6 | M5 — single crossing, `p_disc = 0.37` | **linear in `n`, and already measured to be ~2–4x high** | Pilot says 0.08–0.30 over five batches. Does not change the verdict (§5). |
| 7 | duty cycle 0.10 | **more than** factor 32 in `n` across 0.05–0.20 | Dominant uncertainty in the *amplitude* numbers (D55), and `pi` only goes as `√L`. But §5.2 shows the `rho` correction is itself ducy-dependent and **compounds in the same direction**: at larger ducy the losses are smaller *and* the correction is larger (`<k²>` ratio 1.03 → 5.35 across 0.05 → 0.20), both pushing `pi` toward ½. So the factor 32 is now a **lower bound**, not an estimate. Not recomputed at 0.05 or 0.20. |
| 8b | **`reach(s)` = the ladder's `succ_h1` is the right stage weighting** | **factor 8 in `n`** (340 to 2 840) | **NEW, and unchecked.** The weight on each stage is a modelled survival curve of the same single-leaf Gaussian family that is measured wrong by 4.5 S/N elsewhere (§9). It is now the largest unvalidated lever in the calculation — bigger than row 3 was. Testing it needs the real per-stage survival profile, which nothing on this branch measures. |
| 8 | the Viterbi ladders are a fair equal-`P_d` comparison | would void the result | ~~Checked against real data: 3/24 in both arms at S/N 12 (§6.3). The one assumption that has been independently verified.~~ **WITHDRAWN as a verification (§6.6).** Batch A ran at `max_sugg` = 2^14, where **both** arms were saturating, so the realised cut in each was the ratchet and not the ladder. The 3/24 agreement is real but cannot be credited to the ladders, and in the `quadrature` arm the ladder has never been the operative cut at any buffer tried. This assumption is **unverified**, not verified, and no longer the exception in this table. |

---

## Reproducers

- `injection_power.py` — the power calculation. `--geometry` recomputes the per-cell
  losses and `rho_AB` (~5 min), `--strata N` measures the stratification, no flag
  re-derives the tables from `injection_power.json`.
- `injection_pilot.py` — pairing check, per-pair cost, and `p_disc`.
  `--make` / `--arm` / `--report`, with `--max-sugg`.
- `rho_check.py` — §5.2. Measures the score correlation against the profile overlap on
  the same cells and pushes the corrected value back through `power`/`mcnemar_n`.
  Self-checks the sampler against a closed form first. `rho_check.json` holds every cell.
- `paired_max_sugg.py` — §6.5. The four-cell paired buffer experiment, with the
  analysis pre-committed in `plan()` and written to the run directory before the first
  cell runs. `--plan` / `--make` / `--run` / `--analyse` / `--selftest` (which checks
  the analysis detects a planted bias and does *not* fire on a symmetric one).
  `paired_max_sugg_results.json` holds all 200 runs plus the pre-registration.
- `saturation_sweep.py` — §6.6. The buffer sweep, reporting saturation percentiles and
  the growth exponent `d log2(ncand)/d log2(max_sugg)` that decides whether a feasible
  buffer exists. `--source` builds a paired subset from an existing realisation set.
  `saturation_sweep_results.json` holds the rows.
- `rho_strata.json` — §5.1's 24 positions, each with the uncorrected and corrected `pi`.
- `injection_pilot_results.json` — every pilot run behind §6.3, per run, including the
  `saturation` column.
- `schemes/cheby_aggressive.npz`, `schemes/cheby_quadrature.npz` — the recalibrated
  Chebyshev ladders, `P_d = 0.1031` both.

Nothing inherited from the `metric-gridding` branch was modified. `nearest_template.py`,
`nearest_template_cheby.py` and `amplitude_loss.py` are used as-is, and no defect was
found in them. The one library change made on this work — recording the effective
pruning threshold — was built on a **separate branch off `upstream/main`**
(`upstream-max-sugg-logging`) and is [PR #14](https://github.com/pravirkr/pyloki/pull/14);
`src/` in this worktree was never touched while a measurement was running.

---

## 11. Handoff

Written at hand-off, mid-experiment. Read this before trusting anything above.

### 11.1 What I verified myself vs. what I took on trust

**Verified in this session, from the code or by measurement:**

- The score is a unit-variance matched-filter S/N. Checked the boxcar template in
  `scoring.py:136-139` is exactly zero-mean and unit-L2 for several `(size, width)`,
  and read the per-bin whitening at `scoring.py:287`. M1 is not an assumption.
- The threshold comparison is a hard per-leaf cut (`prune.py:294`) and the overflow
  ratchet is real (`world_tree.py:527-551`, read directly).
- The noise is unseeded (`pulse.py:338`, read directly) and **pairing by persisting the
  `TimeSeries` works** — see §11.2.
- The **Chebyshev ladders are mine**, generated this session
  (`schemes/cheby_*.npz`), not inherited. Both land at `P_d = 0.1031`, and batch A of
  the pilot reproduced that on real data (3/24 in both arms at S/N 12).
- The branching patterns, from `generate_branching_pattern(kind="poly_chebyshev_moving")`:
  `aggressive` 1.506e12, `quadrature` 5.405e20.
- Per-pair wall clock, §11.3.
- All the geometry in §4.1 and §5.1 — but see the caveat below on *whose* code produced it.

**Taken on trust from `metric-gridding`, NOT re-derived:**

- `nearest_template.py` / `nearest_template_cheby.py` are correct. I used
  `build_sets_cheby` and `min_excursion` as black boxes. I did not re-run
  `tests/test_nearest_template.py`, did not re-check the Minkowski-sum argument, and did
  not repeat the brute-force agreement or the window-monotonicity sweep that D47 says any
  such measurement must carry. **Everything in §4.1 and §5.1 inherits whatever
  correctness those modules have.** Given D47, a successor should re-run that test suite
  before quoting §4.1.
- `amplitude_loss.snr_ratio` is correct, including the 16-position phase average (D54).
  I used it unmodified and extended it to a use it was not written for — see §11.5.
- D60's ≥1.56x and 87/90, D62's `branch_max` verdicts, and the D55 loss figures. My own
  §4.1 medians (1.90% / 1.10% matched-filter) sit close to the report's boxcar
  2.03% / 1.14%, which is weak corroboration, not a check.
- I did **not** re-derive the Taylor results at all; this arm is Chebyshev only.

### 11.2 Can the pipeline be seeded to reproduce a noise realisation across arms?

**Not by seeding — but yes, by persistence, and this is verified.**

There is no `seed` or `rng` parameter anywhere in the chain: `PulseSignalConfig`,
`PulseSignalConfig.generate`, `DynamicProgramming`, `Pruning`, `prune_dyp_tree` all lack
one, and `pulse.py:338` is a bare `np.random.default_rng()`. `np.random.seed` does not
help — `default_rng` ignores the legacy global seed. The injected *amplitude* also varies
run to run, because `calibrate_scale_on_folds` (`pulse.py:392`) calibrates against the
realised noise.

What works, and what `injection_pilot.py` does: generate the `TimeSeries` once, persist
`(ts_e, ts_v, dt)` with `np.savez`, and rebuild `TimeSeries(ts_e, ts_v, dt)` in each
arm's process. **Verified:** the saved array round-trips to an identical SHA-256, and two
separate processes running the same arm on it produced bit-identical sorted score vectors
(`cc66644060fdc8f2` twice for `aggressive`, `9ccd6c6400eb642a` twice for `quadrature`).
Nothing downstream of the time series is stochastic.

One caveat that must travel with this: `DynamicThresholdScheme.__init__`
(`thresholding.py:740`) is *also* unseeded, so a ladder re-derived per run would differ
between arms. Generate once, commit the array, reuse — which is what `schemes/cheby_*.npz`
are for. Storage is the real constraint: 32 MB per realisation, so a campaign must stream
(generate → run both arms → delete) rather than persist a corpus.

### 11.3 Per-run wall clock: how it was measured, and the number

`injection_pilot.py` times each run from `DynamicProgramming` construction through
periodogram load, and `--report` prints the **median over runs 2..n**, discarding the
first because numba JIT costs ~25 s once per process.

    python docs/metric_gridding/injection_pilot.py --dir D --make --n 20 --snr 14
    python docs/metric_gridding/injection_pilot.py --dir D --arm aggressive --max-sugg 262144
    python docs/metric_gridding/injection_pilot.py --dir D --arm quadrature --max-sugg 262144
    python docs/metric_gridding/injection_pilot.py --dir D --report

| `max_sugg` | `aggressive` | `quadrature` | **pair** |
|---|---|---|---|
| 2^14 | 0.7–1.4 s | 3.1–4.0 s | 3.8–5.1 s |
| 2^18 | 2.3 s | 47.5 s | **49.8 s** |

Plus ~0.15 s to generate each realisation. **Budget from the 2^18 row** (§6.4), and
expect it to rise again at whatever `max_sugg` finally stops binding — `quadrature`'s
cost scaled 13x for a 16x buffer increase, roughly linearly.

### 11.4 The phase-error-to-score mapping: concluded, not blocked

Stated in full as M3 in §3 and in `injection_power.py`'s module docstring:

    score_X(s) = mu(s) · (1 − L_X(s)) + N_X(s),     N_X ~ N(0, 1)

with `L_X` the fractional S/N loss from `amplitude_loss.snr_ratio` — the exact
profile-smearing model, not a sup-norm proxy (D56: a sup-norm overestimates by ~3x).
Multiplicative in amplitude, additive in noise, because folding with a different
ephemeris averages the same samples so the noise level is unchanged.

**The part a successor should attack is not M3, it is M4.** `D(s) = score_B − score_A`
has a stochastic term of sd `√(2(1−rho_AB))` from the two arms scoring *different*
templates, and measured it is **three times the deterministic term** (0.133 vs 0.044 at
ducy 0.10) — **corrected in §5.2 to 0.230 vs 0.044, i.e. ~5x**. That ratio is what sets `pi`, and `n ∝ (pi − ½)^-2`, so it is the single
assumption the whole answer turns on. Two specific weaknesses:

1. ~~`rho_AB` is computed by feeding `delta_A − delta_B` to `snr_ratio`, which gives the
   *profile* overlap.~~ **CHECKED, and it does not hold — see §5.2.** The boxcar bank
   does not preserve the profile overlap: both quantities are the same `S(k)` weighted
   by `|P(k)|²` and `|H(k)|²` respectively, and the selected boxcar's `<k²>` is 4.10x
   the pulse's. Median `1 − rho` at ducy 0.10 is **0.0265, not 0.0088**, and `pi` falls
   from 0.697 to **0.628**, raising the pairs needed at `p_disc = 0.30` from 161 to 391
   on that aggregate — and to **954** once §9's mid-stratum `pi` and the null stratum's
   share are accounted for, which is the operative number. `rho_check.py`.
2. `L_A` and `L_B` are both **upper** bounds (the search minimises sup-norm over 30
   retained leaves, not loss), and the bounds are not symmetric, so their *difference* is
   not sign-proven. This is D52's caveat, inherited.

### 11.5 Defects found in inherited code

**None.** `nearest_template.py`, `nearest_template_cheby.py` and `amplitude_loss.py` were
used unmodified and nothing in them misbehaved. Two things to flag as *use* notes rather
than defects, since another live session owns that code:

- `nearest_template_cheby` + `quadrature` at late stages (S = 54, 62) exhausts its node
  budget and returns incumbents as large as **724 tolerances**. That is correct
  behaviour — a budget-limited result is a valid upper bound and it is flagged
  `exact=False`, exactly as D50 designed — but it means 7 of 90 cells had to be dropped,
  and anyone who ignores the `exact` flag will get nonsense. `injection_power.py` drops
  them.
- `amplitude_loss.snr_ratio` is used here for something it was not written for: the
  overlap *between two templates*, via `snr_ratio(delta_A − delta_B, ...)`. That is
  legitimate because the phase residual is linear in the coefficient vector, but it is my
  extension, not the original author's, and §11.4 (1) is the check it still needs.

The **cached schemes** under `schemes/` deserve a warning: `aggressive.npz`,
`conservative.npz` and `metric.npz` are **Taylor** (`kind="poly_taylor_moving"`) and there
is no `quadrature` entry at all. They are not valid for the Chebyshev arm and using them
would silently compare the wrong ladders. `cheby_aggressive.npz` and
`cheby_quadrature.npz` are the correct ones and are new in this commit.

### 11.6 The single most useful thing to do first — **DONE, see §6.5**

> **Closed.** Run as specified: 50 realisations, all four cells, pre-registered analysis
> committed before the first cell. The answer is split — the buffer *pressure* is
> arm-dependent (p = 0.012) but the *outcome* bias is unresolved with a point estimate
> (+0.080) above the campaign's effect of interest (0.06). **§6.4 is not closed**, the
> campaign must not run at 2^18, and closing it needs n ≈ 211 per cell. The original
> text follows.

**Finish the paired `max_sugg` experiment** — both arms, both buffers, on one fixed set
of ~50 realisations at S/N 14 — and answer the one question §6.4 left open: does the
buffer penalise `quadrature` more than `aggressive`? Everything downstream waits on it.
If it does, no campaign at any `max_sugg` either arm saturates means anything. If it does
not, the design in §9 can be run as written.

The script exists and is committed; the driver that was killed mid-run is the loop in
`injection_pilot.py` called at four `(arm, max_sugg)` combinations on one `--make`
directory, with `arm_*.json` renamed between runs. Budget ~35 min for
`quadrature`-at-2^18 over 20 pairs, ~90 min over 50. Partial results are in
`injection_pilot_results.json` under `paired_max_sugg_test`: `aggressive` is complete at
both buffers (11/20 → 13/20), `quadrature` only at 2^14 (11/20).

Second: the `rho_AB` check in §11.4 (1). It is cheap, it has never been done, and `pi`
— hence `n` — depends on it more than on anything else.

> **Also done — §5.2.** It had never been checked and it was wrong: the code used the
> profile overlap where the model needs the score correlation. `pi` 0.697 → 0.628, and
> through the stratification it deleted one of the three strata (§5.1, §9).

### 11.7 What the second session left open

1. ~~**The buffer.** Raise `max_sugg` past 2^18 until the 99th percentile of saturation
   is below 0.9 in both arms.~~ **Attempted and it fails — §6.6.** There is no such
   buffer for `quadrature`. The de-confound at n ≈ 211 is therefore not reachable
   either, since it would have to run at a buffer where one arm still saturates. The
   stricter ladder, which was the only remaining route keeping both the deployed
   configuration and the idealised question, **was then tested and fails (§6.7)**. What
   is left is §6.6's options 2 and 3 and the buffer-policy change of §6.7.
2. **The operating point.** §6.5 measures `P_d` = 0.76/0.80 at S/N 14 and 2^18, well
   above the 0.5 §9 wants, so the 3-point pilot should bracket *below* S/N 14. ~~This
   contradicts the `metric-gridding` session's D78 (S/N 15–17).~~ **Resolved:** that
   session withdrew D78 on finding its own discordance peak sits at model `P_d` = 0.567
   — the same criterion in different variables, the disagreement being entirely in its
   `S/N ↔ P_d` mapping. §9 now records why the criterion is stated in the observable.
3. **ducy 0.05 and 0.20** are still on the uncorrected `rho_AB` (§5.2). The
   leading-order factor says 0.05 barely moves and 0.20 moves more than 0.10 did.
4. **Row 2 of §10** is retired for stages 1–10 *on the `metric-gridding` session's
   measurement, not mine.* If it has to carry weight, re-run it here.
5. ~~**The upstream logging patch** is written but not applied or verified.~~ **Done and
   shipped** — implemented on `upstream-max-sugg-logging` off `upstream/main`, verified
   both directions with a regression test, and merged upstream as
   [PR #14](https://github.com/pravirkr/pyloki/pull/14).

6. **The survival profile is inherited and unvalidated (§4.3 ⚑).** Its shape carries
   §4.3's reversal, the `n` range and the strata critique; its absolute levels carry
   nothing and should not be used. Replicating it here would settle the `n` range — but
   note it can only be replicated honestly in the `aggressive` arm, since `quadrature`
   cannot be run at a non-binding buffer at all (§6.6), so a "both arms" replication
   would measure the ratchet in one of them.

7. **`DECISION_STAGES` = [2, 6, 10, 14, 20, 28] is misplaced.** It sits almost entirely
   *before* the measured decision window (first loss at level 13, median 27), so §5.1's
   stratification and §9's strata — including the null stratum the design leans on —
   are built over stages where little is decided. Not recomputed here, because that
   would put a new headline on an unverified inherited measurement. First thing to redo
   if the configuration question reopens.

8. ~~**The stricter-ladder test (§6.6, option 4).**~~ **Done — §6.7, and it fails.** A
   tenfold reduction in the detection target buys 4 ratcheted levels out of 192 and
   leaves zero clean runs, because the effective cut is set by `top-K` rather than by
   the scheme: raising the ladder by 0.87 moved the cut that actually ran by −0.06.

### 11.8 How this ended

**The campaign was cancelled** (2026-09-18, assaferan) on the strength of §6.6 and
§6.7: the comparison it was designed to make cannot be made at this configuration,
because `quadrature` cannot run at its own threshold scheme at any affordable buffer
and no recalibration of that scheme changes it. That is a result, not a failure to get
one, and it needed no campaign to establish.

**What went out.** The silent-ratchet defect found along the way is
[PR #14](https://github.com/pravirkr/pyloki/pull/14), merged upstream. The
buffer-policy finding is drafted as `07_upstream_buffer_policy.md` — an *issue* rather
than a patch, since it asks a question about pruning semantics — and is **unposted**,
awaiting review. This branch is pushed to `assaferan/pyloki` so that draft's reproducer
links resolve; nothing has been posted to the upstream tracker.

**What deliberately did not go out.** Anything about tiling. Per §8 nothing in this
document currently supports an upstream sentence about it, and the buffer result is
stronger standing alone. `04_upstream_report.md` on the `metric-gridding` branch is a
separate, still-unposted document owned by that session; if both ever go out they need
to be consistent with each other.

**If someone reopens this**, items 3, 4, 6 and 7 above are the live ones, and item 7 —
`DECISION_STAGES` sitting before the measured decision window — is the first, because
the stratification and the strata in §9 are conditioned on it.
