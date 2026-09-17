# 05_injection_design.md — power calculation for the `tiling_strategy` injection campaign

Status: **design and power calculation only. No production campaign was run**, and this
document recommends that none be run in its current form. Everything below is either
measured in this session or derived from a model whose assumptions are listed, with a
sensitivity, in the assumption register at the end.

Arm analysed: **Chebyshev basis, `aggressive` vs `quadrature`, at the shipped
`branch_max = 16`** — the only pairing with a proven geometric gain (D60: `quadrature`
strictly closer in 87/90 cells, ≥1.56x) that builds at the default (D62), and therefore
the only one that tests behaviour a user gets without changing config. Config throughout
is Phase 3's: 268.4 s / 64 segments / `poly_order = 4` / `N_b = 64` / `eta = 1`.

---

## 1. Verdict

**The statistics are reachable. The systematics are not — and the systematic runs the
other way.**

- **Power:** 80% power at two-sided α = 0.05 needs **≈ 130 to 2 400 paired runs**
  across the plausible range of the effect size, and **≈ 5 × 10⁵ pairs in the worst
  corner I can construct**. At the measured **3.8–4.9 s per pair** that is
  **0.2 to 3 core-hours** in the central case and ~600 core-hours in the worst corner.
  By the standard the brief set — can the effect be resolved at feasible cost — the
  answer is **yes, comfortably**, and it is not close.

- **But the pilot found a larger, arm-dependent artefact that points the wrong way.** At
  the operating point where discordance is greatest (injected S/N 14, `P_d ≈ 0.35`), 20
  paired runs gave `aggressive` 7 recoveries to `quadrature`'s 4, with the discordances
  splitting **4 for `aggressive` against 1 for `quadrature`** — the opposite of the
  hypothesis. The cause is visible in the same run: the candidate buffer saturated at
  **0.978 of `max_sugg`**, and on overflow `world_tree.py:527-551` ratchets the
  effective cut to `max(stage threshold, top-K, median)`. `quadrature` emits ~2^18 more
  branches than `aggressive` (branching-pattern product 5.4e20 against 1.5e12), so it
  saturates first and gets silently thresholded harder than its calibrated ladder says.
  `PruneStats` records only the nominal threshold (`prune.py:685`), so this is invisible
  in the log.

**Recommendation: do not run the campaign yet.** Not because it cannot be resolved —
it can — but because at `max_sugg = 2^14` it would resolve the buffer ratchet, publish
it as a tiling result, and become the fourth retraction on this branch. §7 says what has
to change first. §8 says what can be concluded without running it.

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

### 4.3 Where the decision is actually made — and it is not where the geometry is

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

---

## 5. The power calculation

Aggregating §4 over stages, weighting each by `reach(s) · phi(T(s) − mu(s)(1−L_A))`
(the probability of being in the race times the density of arm A's score at the
threshold), and taking `pi = E[w·D⁺] / E[w·|D|]`:

| ducy | `pi` | `p_disc` (modelled, M5 upper bound) | discordant pairs needed | **pairs needed** | core-hours at 4.5 s/pair |
|---|---|---|---|---|---|
| 0.05 | 0.787 | 0.699 | 22 | **31** | 0.04 |
| 0.10 | **0.697** | 0.370 | 49 | **131** | 0.16 |
| 0.20 | 0.607 | 0.168 | 168 | **997** | 1.2 |

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

1. **`p_disc` barely matters.** Moving it by a factor of 30 (0.30 → 0.01) moves `n` by
   the same factor, i.e. from hours to a day. The verdict does not turn on it — which is
   fortunate, because it is the parameter the model estimates worst (M5).

2. **`pi` is the whole question.** `n ∝ (pi − ½)^-2`, so `pi = 0.52` instead of 0.70
   costs a factor of 105. **`pi` is the single assumption the conclusion is most
   sensitive to**, and §4.1 says what drives it: the ratio
   `mu·(L_A − L_B) / √(2(1−rho_AB))`.

   The reassuring part is that this ratio is **more robust than it looks to the duty
   cycle**, which D55 identified as the dominant uncertainty in the amplitude numbers.
   Both the numerator and the denominator scale with the loss, so the ratio goes like
   `√L` rather than `L`: a factor of 15 in the losses between ducy 0.05 and 0.20 becomes
   a factor of 1.3 in `pi` (0.787 → 0.607) and a factor of 32 in `n` (31 → 997). All
   three remain trivially affordable.

### 5.1 Stratification: measured, and it does not stratify

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

---

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

| injected S/N | `aggressive` s/run | `quadrature` s/run | **s/pair** |
|---|---|---|---|
| 12 | 0.67 | 3.11 | **3.78** |
| 14 | 1.01 | 3.54 | **4.55** |
| 16 | 1.38 | 3.54 | **4.92** |

Plus ~0.15 s to generate each realisation. So **≈ 5 s per pair, single-threaded**, and
10 000 pairs is **14 core-hours**. This is not an expensive experiment.

### 6.3 What the pilot actually found

Three batches of 20–24 pairs at `max_sugg = 2^14` (the value `run_injections.py` uses),
with the per-strategy recalibrated Chebyshev ladders:

| injected S/N | `aggressive` recovered | `quadrature` recovered | discordant | quad-only | agg-only | max buffer saturation |
|---|---|---|---|---|---|---|
| 12 | 3/24 | 3/24 | 2/24 | 1 | 1 | 0.880 |
| 14 | **7/20** | **4/20** | 5/20 | 1 | **4** | **0.978** |
| 16 | 16/20 | 17/20 | 3/20 | 2 | 1 | **0.986** |

Two positives and one showstopper.

*Positive 1 — the ladders are honestly calibrated.* At S/N 12, both arms recovered
exactly 3 of 24 (12.5%), against a nominal `P_d = 0.1031`. The equal-`P_d` recalibration
does what it claims, in both arms, on real data. That is the one thing all three
retracted results never checked.

*Positive 2 — `p_disc` is real and measurable.* 2/24, 5/20, 3/20, i.e. 0.08–0.25 against
the model's 0.37. M5 is high by about 4x, as expected. At `p_disc = 0.10` and
`pi = 0.70`, `n = 467` pairs — 40 minutes on one core.

*Showstopper — the buffer ratchet dominates, and it points the wrong way.* At S/N 14,
the operating point closest to `P_d = 0.5` and therefore the one a campaign would
choose, `aggressive` beat `quadrature` 7 to 4, with discordances splitting 4-to-1 **for
`aggressive`**. That is not the hypothesised effect and it is not noise at the observed
rate; it is the candidate buffer. `world_tree.py:527-551` raises the effective cut to
`max(stage threshold, top-K, median)` on overflow, and `quadrature` — with 2^17.9 times
`aggressive`'s branching — overflows first. Saturation hit 0.978 and 0.986. It is
**silent**: `PruneStats` logs only the nominal threshold (`prune.py:685`), so nothing in
the run output would have flagged it, and a campaign run as designed in the brief would
have reported it as a tiling result with a p-value attached.

---

## 7. What has to change before this is worth running

1. **Eliminate the buffer ratchet, and prove it is gone.** Raise `max_sugg` until
   neither arm saturates in any run, and record `saturation = ncand / max_sugg` per run
   as a first-class output (`injection_pilot.py` already does). Any run with saturation
   above ~0.9 is not measuring the scheme's thresholds and must be excluded or the
   batch rerun. **This is a cost question, not a correctness question**: `quadrature`'s
   branching product is 5.4e20, so a buffer large enough not to bind may make its runs
   far more expensive than the 3.5 s measured here, and the feasibility verdict has to
   be recomputed at whatever `max_sugg` turns out to be sufficient. See §9.
2. **Pre-screen the null stratum** with `nearest_template_cheby.py` (§5.1) rather than
   assuming random positions give one. Only ~1 in 8 does.
3. **Stream the realisations** instead of persisting them (§6.1), or budget 160 GB.
4. Keep the ladders frozen and committed. Do not re-derive `DynamicThresholdScheme` per
   run; it is unseeded.

---

## 8. What can be concluded without running it

- The equal-`P_d` recalibration works on real data: both arms recovered 3/24 at S/N 12
  against a nominal 0.1031. The threshold ladder is not the confound it was in Taylor
  (7.10 vs 7.00 here, against 7.70 vs 9.10 there).
- The discordance rate is 8–25% at plausible operating points, so a paired design is
  the right one and an unpaired one would be hopeless.
- **Pairing on the noise realisation does not make the comparison deterministic.** The
  template-to-template noise difference (sd 0.133 at ducy 0.10) is ~3x the amplitude
  difference it is meant to reveal. Any future design on this branch that treats a
  paired run as a clean A/B should be checked against this number first.
- **The tiling geometry is largest where nothing is being decided.** `pi` per stage runs
  0.52–0.56 over stages 2–10, where 90% of the `P_d` is lost, and 0.87–0.93 over stages
  40–50, where the cumulative survival is already flat. Any mechanism that only bites
  late in the tree cannot move the final `P_d` much, whatever its size. This is a
  structural statement about the moving-grid ladder and it does not need an injection
  campaign to support it.
- A defensible upstream sentence, without any campaign: *at equal `P_d`, the tiling
  choice moves the survival decision by ~0.04 in score at the stages where the decision
  is made, against a 0.13 stochastic difference from the change of template — so its
  effect on end-to-end sensitivity is smaller than the geometric gain (≥1.56x in phase)
  suggests, and is not where `quadrature`'s 2^17.9 cost is going.*

---

## 9. If it is run anyway: the design, pre-committed

**Arms.** `poly_basis="chebyshev"`, `use_moving_grid=True`, `prune_poly_order=4`,
`branch_max=16`; `tiling_strategy` ∈ {`aggressive`, `quadrature`} and nothing else
differs. Config from `phase3_config.make_config`.

**Pairing.** One `TimeSeries` per pair, generated once, persisted as `(ts_e, ts_v, dt)`,
consumed by both arms, deleted after both have run. One arm per process. Verified
bit-deterministic (§6.1).

**Thresholds.** `schemes/cheby_aggressive.npz` and `schemes/cheby_quadrature.npz`,
committed, Viterbi at `P_d = 0.1`, never re-derived.

**Operating point.** Injected S/N chosen so on-grid `P_d ≈ 0.5` in *both* arms, since
discordance is maximised at the steepest part of the recovery curve. From §6.3 that is
between 14 and 16 at this config; fix it by a 3-point pilot **at the final `max_sugg`**
and then do not change it.

**`max_sugg`.** Raised until the 99th percentile of `ncand / max_sugg` over a 50-pair
pilot is below 0.9 in both arms. Recorded per run. Runs above 0.9 are excluded and the
exclusion count is reported.

**Strata.** Three, assigned by pre-screening candidate positions with
`nearest_template_cheby.py` over the decision stages (2–28) and binning on the
decision-weighted `pi` of §5.1:
- **null** (`pi < 0.55`, ~1 position in 8): must show no strategy difference.
- **mid** (0.55 ≤ `pi` ≤ 0.70, the bulk).
- **effect** (`pi > 0.70`, the upper quartile; a `pi > 0.75` stratum is not reachable).

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

**n and stopping rule.** Fixed n, no interim looks: n = 2 000 pairs, which gives ≥80%
power for any `pi ≥ 0.60` at the pilot's `p_disc ≥ 0.08`, and ~3 core-hours at 5 s/pair.
Recompute n from §5 once `max_sugg` is settled and the per-pair cost is re-measured. **No
optional stopping** — this branch has retracted three results and a stopping rule that
looks at the data is how a fourth would happen.

**Pre-committed analysis.** Primary: two-sided McNemar (exact binomial on the discordant
pairs) of `pi = 1/2`, pooled over the mid and effect strata, α = 0.05. Secondary,
descriptive only, not tested: per-stratum `pi`, and the monotonicity of `pi` across the
three strata. Reported unconditionally: the null-stratum result, the saturation
distribution, the exclusion count, and `n01` and `n10` separately. Effect size of
interest: `pi = 0.60`, i.e. 3 discordances for `quadrature` to every 2 for `aggressive`.

---

## 10. Assumption register, ordered by how much the answer moves

| # | assumption | if wrong, `n` moves by | note |
|---|---|---|---|
| 1 | **`pi ≈ 0.70`** (§5), i.e. `mu·(L_A−L_B) / √(2(1−rho_AB)) ≈ 0.4` | **up to 10⁵x** — `n ∝ (pi−½)^-2`, and `pi → ½` makes it unbounded | **the one to worry about.** Everything else is second order. |
| 2 | `L_A` and `L_B` are the losses of the leaf that actually scores | factor of a few, **sign not guaranteed** | Both are *upper* bounds (M3, D52 convention): the search minimises sup-norm, not loss, over 30 retained leaves. If `aggressive`'s true best leaf is better than retained, `pi → ½`; if `quadrature`'s is, `pi` rises. The two are not symmetric and the difference is not proven. |
| 3 | `rho_AB` from the between-template phase residual is the score correlation | factor ~3 in `(pi−½)`, so ~10x in `n` | The smearing model is exact for the *profile* overlap (D53); treating it as the *score* correlation assumes the boxcar filter preserves it. Not checked. |
| 4 | M6 — survival is decided by one covering leaf, not a max over many | raises `pi`, lowers `n` | Deliberately conservative for the stated question, but it is also the channel the pilot's buffer artefact rides on. |
| 5 | M2 — `mu(s) = snr_final·√((s+1)/nseg)` | ~20% in `pi` | Same model the shipped ladder uses, so an error here is an error in the ladder too. |
| 6 | M5 — single crossing, `p_disc = 0.37` | **linear in `n`, and already measured to be ~4x high** | Pilot says 0.08–0.25. Does not change the verdict (§5). |
| 7 | duty cycle 0.10 | factor 32 in `n` across 0.05–0.20 | Dominant uncertainty in the *amplitude* numbers (D55), but `pi` only goes as `√L`, so it does not threaten feasibility. |
| 8 | the Viterbi ladders are a fair equal-`P_d` comparison | would void the result | Checked against real data: 3/24 in both arms at S/N 12 (§6.3). The one assumption that has been independently verified. |

---

## Reproducers

- `injection_power.py` — the power calculation. `--geometry` recomputes the per-cell
  losses and `rho_AB` (~5 min), `--strata N` measures the stratification, no flag
  re-derives the tables from `injection_power.json`.
- `injection_pilot.py` — pairing check, per-pair cost, and `p_disc`.
  `--make` / `--arm` / `--report`.
- `schemes/cheby_aggressive.npz`, `schemes/cheby_quadrature.npz` — the recalibrated
  Chebyshev ladders, `P_d = 0.1031` both.

Nothing inherited from the `metric-gridding` branch was modified. `nearest_template.py`,
`nearest_template_cheby.py` and `amplitude_loss.py` are used as-is, and no defect was
found in them.
