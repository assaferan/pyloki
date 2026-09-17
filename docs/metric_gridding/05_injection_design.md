# 05_injection_design.md — power calculation for the `tiling_strategy` injection campaign

> ## ⚠ WORK IN PROGRESS — handed off mid-experiment
>
> Written in one session and stopped before the last experiment finished. The power
> calculation (§2–§5) is complete and self-consistent. The empirical part (§6) is not:
> one confound was found, partly characterised, and **not** pinned down, and an early
> draft of this file over-claimed it from a comparison that could not support the claim
> (§6.4 says exactly what and why). **§11 is the handoff** — read it before trusting any
> number here, because it separates what was verified in this session from what was
> taken on trust from the `metric-gridding` branch.
>
> No production campaign was run. No inherited code was modified.

Status: **design and power calculation only.** Everything below is either measured in
this session or derived from a model whose assumptions are listed, with a sensitivity,
in the assumption register (§10).

Arm analysed: **Chebyshev basis, `aggressive` vs `quadrature`, at the shipped
`branch_max = 16`** — the only pairing with a proven geometric gain (D60: `quadrature`
strictly closer in 87/90 cells, ≥1.56x) that builds at the default (D62), and therefore
the only one that tests behaviour a user gets without changing config. Config throughout
is Phase 3's: 268.4 s / 64 segments / `poly_order = 4` / `N_b = 64` / `eta = 1`.

---

## 1. Verdict

**Reachable. ~30 core-hours for the design in §9 — but the measurement is not yet
trustworthy, because `max_sugg` changes outcomes on identical data and it has not been
shown to do so symmetrically between the arms.**

- **Power.** At the discordance rate the pilot measured (`p_disc` = 0.08–0.30), 80%
  power at two-sided α = 0.05 needs **156 pairs if `pi` = 0.70, 647 if 0.60, 2 609 if
  0.55**, and 1.6 × 10⁴ in the pessimistic corner `pi` = 0.52 (all at `p_disc` = 0.30;
  scale inversely for smaller). The pre-committed design is **n = 2 000 pairs**.

- **`pi` is 0.63, not 0.70 — `rho_AB` was the wrong quantity (§5.2, new).** The code set
  `rho_AB` to the *profile* overlap `snr_ratio(delta_A − delta_B)`; M4 needs the *score*
  correlation. They are the same smearing factor `S(k)` weighted by the pulse's spectrum
  and the filter's respectively, and the bank's selected boxcar has 4.10x the pulse's
  `<k²>`. Median `1 − rho` at ducy 0.10 is **0.0265, not 0.0088**; `pi` falls 0.697 →
  **0.628** and the pairs needed at `p_disc` = 0.30 rise **161 → 391** (1 464 at the
  pilot's worst `p_disc` = 0.08). **`n = 2 000` still covers it**, so this makes the
  design adequate rather than comfortable. It was §10's row 3 and is now measured; row 1
  is still the one to worry about.

- **Cost, measured.** At `max_sugg = 2^18` (the library default) a pair costs **49.8 s**
  — `aggressive` 2.3 s, `quadrature` 47.5 s — so 2 000 pairs is **28 core-hours** and
  the pessimistic corner ~230. At `max_sugg = 2^14` it is 3.8–5.1 s/pair, but §6.3 says
  that setting is not sound. There is no statistical or computational obstacle.

- **What is *established* about the confound.** On a **single fixed noise realisation**,
  changing only `max_sugg` from 2^14 to 2^18 changed `aggressive` from **0 surviving
  candidates to 110** and `quadrature` from 8 553 to 134 418. Across a **fixed set of 20
  realisations**, `aggressive` recovered **11/20 at 2^14 and 13/20 at 2^18**, flipping 2
  runs, both toward recovery. So the candidate buffer is not a neutral knob: it changes
  the search outcome on data the search has already seen, at both the shipped default
  and the value the existing injection driver uses.

  The mechanism is `world_tree.py:527-551`: on overflow the effective cut ratchets to
  `max(stage threshold, top-K, median)`, and `PruneStats` logs only the nominal
  threshold (`prune.py:685`), so it is silent. A direct probe showed the *final*-stage
  cut was nominal in all four runs tested (minimum surviving score below
  `thresholds[-1]`), so the ratchet is firing at **intermediate** stages, which is where
  it matters and where it is hardest to see.

- **What is NOT established, and what an early draft of this file wrongly claimed.** That
  the confound is *arm-dependent* — i.e. that it penalises `quadrature` more than
  `aggressive` because `quadrature`'s branching product is 2^28.4 times larger. It is
  plausible from the code and from `quadrature`'s branching, but the experiment that
  would show it (both arms, both buffers, identical realisations) **was killed before
  the `quadrature`-at-2^18 arm finished**. See §6.4.

**Recommendation: do not run the campaign until §6.4 is closed.** Not because the effect
is out of reach — it is not — but because the one systematic large enough to reverse the
answer is currently uncharacterised, and 20-pair batches are far too noisy to
characterise it (§6.4 shows two batches at *identical settings* giving `quadrature`
4/20 and 11/20).

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

| ducy | `pi` | `p_disc` (modelled, M5 upper bound) | discordant pairs needed | **pairs needed** | core-hours at 50 s/pair |
|---|---|---|---|---|---|
| 0.05 | 0.787 | 0.699 | 22 | **31** | 0.4 |
| 0.10 | **0.697** | 0.370 | 49 | **131** | 1.8 |
| 0.20 | 0.607 | 0.168 | 168 | **997** | 13.8 |

At the pilot's *measured* `p_disc = 0.30` (§6.3) rather than the modelled one, the same
three rows are n = 74 / 163 / 560, i.e. 1.0 / 2.3 / 7.8 core-hours.

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
needed at the measured `p_disc` rise by ~2.4x. The pre-committed `n = 2 000` (§9) still
covers it — 391 pairs at `p_disc = 0.30` and 1 464 at the pilot's worst 0.08 — so the
verdict moves from comfortable to **adequate**, not from reachable to unreachable. Only
ducy 0.10 was recomputed; the 0.05 and 0.20 rows of §5's table are uncorrected and
should be assumed to move by a similar factor.

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

## 7. What has to change before this is worth running

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
   assuming random positions give one. Only ~1 in 8 does.
3. **Stream the realisations** instead of persisting them (§6.1), or budget 160 GB.
4. Keep the ladders frozen and committed. Do not re-derive `DynamicThresholdScheme` per
   run; it is unseeded.

---

## 8. What can be concluded without running it

- The equal-`P_d` recalibration works on real data: both arms recovered 3/24 at S/N 12
  against a nominal 0.1031. The threshold ladder is not the confound it was in Taylor
  (7.10 vs 7.00 here, against 7.70 vs 9.10 there).
- The discordance rate is 8–30% at plausible operating points, so a paired design is
  the right one and an unpaired one would be hopeless.
- **`max_sugg` changes search outcomes silently, and that is worth reporting upstream on
  its own.** On one fixed time series, `aggressive` returns 0 candidates at
  `max_sugg = 2^14` and 110 at 2^18 (§6.4). The overflow ratchet at
  `world_tree.py:527-551` raises the effective cut, `PruneStats` records only the nominal
  threshold (`prune.py:685`), and the final-stage cut stays nominal so nothing in the
  output betrays it. Logging the effective threshold when `prune_on_overload` fires would
  be a three-line fix. Whether it biases one `tiling_strategy` more than another is a
  separate question and is open (§6.4).
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
2^18` that is near S/N 14 (13/20 and 15/20, §6.3), but the curve moves with the buffer —
the same S/N gave 7/20 and 4/20 at 2^14 — so fix it by a 3-point pilot **at the final
`max_sugg`** and then do not change it.

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

**n and stopping rule.** Fixed n, no interim looks: **n = 2 000 pairs**. That needs
`m = 194` discordant pairs for `pi = 0.60`, `85` for 0.65 and `47` for 0.70, so it gives
≥80% power for `pi ≥ 0.65` across the whole measured `p_disc` range (0.08–0.30), and for
`pi = 0.60` provided `p_disc ≥ 0.10` — which it was, at 0.30, at the operating point the
design picks. Covering `pi = 0.60` at the pilot's worst `p_disc = 0.08` would take 2 425.
At the measured 49.8 s/pair, 2 000 pairs is **28 core-hours**; recompute both from §5
once `max_sugg` is settled, since the cost rises with it. **No optional
stopping** — this branch has retracted three results, and a stopping rule that looks at
the data is how a fourth would happen.

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
| 1 | **`pi ≈ 0.63`** (§5.2; was 0.70), i.e. `mu·(L_A−L_B) / √(2(1−rho_AB))` | **up to 10⁵x** — `n ∝ (pi−½)^-2`, and `pi → ½` makes it unbounded | **still the one to worry about**, and now closer to ½ than it was: row 3 moved it once already. |
| 2 | `L_A` and `L_B` are the losses of the leaf that actually scores | factor of a few, **sign not guaranteed** | Both are *upper* bounds (M3, D52 convention): the search minimises sup-norm, not loss, over 30 retained leaves. If `aggressive`'s true best leaf is better than retained, `pi → ½`; if `quadrature`'s is, `pi` rises. The two are not symmetric and the difference is not proven. |
| 3 | ~~`rho_AB` from the between-template phase residual is the score correlation~~ **MEASURED, and it was wrong** | **2.4x in `n`, already applied** | §5.2. It is the *profile* overlap; the boxcar bank does not preserve it (`<k²>` ratio 4.10). Median `1−rho` 0.0088 → 0.0265, `pi` 0.697 → 0.628. No longer an assumption at ducy 0.10; still one at 0.05 and 0.20. |
| 4 | M6 — survival is decided by one covering leaf, not a max over many | raises `pi`, lowers `n` | Deliberately conservative for the stated question. It is also the channel any arm-dependent `max_sugg` bias would ride on (§6.4). |
| 5 | M2 — `mu(s) = snr_final·√((s+1)/nseg)` | ~20% in `pi` | Same model the shipped ladder uses, so an error here is an error in the ladder too. |
| 6 | M5 — single crossing, `p_disc = 0.37` | **linear in `n`, and already measured to be ~2–4x high** | Pilot says 0.08–0.30 over five batches. Does not change the verdict (§5). |
| 7 | duty cycle 0.10 | factor 32 in `n` across 0.05–0.20 | Dominant uncertainty in the *amplitude* numbers (D55), but `pi` only goes as `√L`, so it does not threaten feasibility. |
| 8 | the Viterbi ladders are a fair equal-`P_d` comparison | would void the result | Checked against real data: 3/24 in both arms at S/N 12 (§6.3). The one assumption that has been independently verified. |

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
- `paired_max_sugg.py` — §11.6. The four-cell paired buffer experiment, with the
  analysis pre-committed in `plan()` and written to the run directory before the first
  cell runs. `--plan` / `--make` / `--run` / `--analyse`.
- `injection_pilot_results.json` — every pilot run behind §6.3, per run, including the
  `saturation` column.
- `schemes/cheby_aggressive.npz`, `schemes/cheby_quadrature.npz` — the recalibrated
  Chebyshev ladders, `P_d = 0.1031` both.

Nothing inherited from the `metric-gridding` branch was modified. `nearest_template.py`,
`nearest_template_cheby.py` and `amplitude_loss.py` are used as-is, and no defect was
found in them.

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
ducy 0.10). That ratio is what sets `pi`, and `n ∝ (pi − ½)^-2`, so it is the single
assumption the whole answer turns on. Two specific weaknesses:

1. ~~`rho_AB` is computed by feeding `delta_A − delta_B` to `snr_ratio`, which gives the
   *profile* overlap.~~ **CHECKED, and it does not hold — see §5.2.** The boxcar bank
   does not preserve the profile overlap: both quantities are the same `S(k)` weighted
   by `|P(k)|²` and `|H(k)|²` respectively, and the selected boxcar's `<k²>` is 4.10x
   the pulse's. Median `1 − rho` at ducy 0.10 is **0.0265, not 0.0088**, and `pi` falls
   from 0.697 to **0.628**, raising the pairs needed at `p_disc = 0.30` from 161 to 391.
   `n = 2 000` still covers it. `rho_check.py`.
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

### 11.6 The single most useful thing to do first

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
