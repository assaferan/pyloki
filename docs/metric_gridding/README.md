# injection-design — what this branch established, and why the campaign was not run

**Status: CLOSED, 2026-09-18.** The `tiling_strategy` injection campaign this branch was
created to design was **cancelled, not run** — the comparison it makes cannot be made at
this configuration. Reopened twice since only to correct the record (2026-09-22).

Entry point. `05_injection_design.md` is 1 500+ lines and **deliberately keeps retracted
claims visible**, because this branch's failure mode is a correct measurement carrying a
framing that does not hold; struck-through text and ⚠/⚑ boxes mark them. Read §1, then
§6.6, then §11. `DECISIONS.md` here stops at 2026-09-17 and does **not** cover the
injection work at all — it is the shared Phase 0–3 log, not this branch's record.

Relation to `metric-gridding`. The two branches share no commits: each rebased its own
copy of the Phase 0–3 work onto `upstream/main` at `18d04b3`. `survival_profile.py` and
that branch's own `README.md` live only there; everything listed under **Files** below
lives only here. Where this document relies on its measurements, it says so and marks
them inherited.

Arm analysed throughout: **Chebyshev basis, `aggressive` vs `quadrature`, at the shipped
`branch_max = 16`** — the only pairing with a proven geometric gain that builds at the
default. Config: 268.4 s / 64 segments / `poly_order = 4` / `N_b = 64` / `eta = 1`.

## How it ended

The campaign was designed, pre-registered, powered, costed — and then blocked, by
something that is not cost:

> **`quadrature` cannot run at its own threshold scheme at any affordable buffer, and no
> recalibration of that scheme changes it.** The comparison requires both arms to run
> against the scheme rather than against the candidate buffer. `aggressive` can;
> `quadrature` cannot, because its candidate count rises with the buffer it is given.

That is a result, and it needed no campaign to establish. It is also the branch's one
transferable finding, and the reason the branch ends here rather than running something
cheaper and calling it an answer.

## Live results

**The candidate set does not converge** (`saturation_sweep.py`, §6.6). Median `ncand`
over **one fixed set of 10 realisations** (the first 10 of §6.5's, so every column is the
same data) swept across four buffers:

| arm | 2^18 | 2^19 | 2^20 | 2^21 | saturation at 2^21 |
|---|---|---|---|---|---|
| `aggressive` | 16 778.5 | **16 778.5** | **16 778.5** | **16 778.5** | 0.008 |
| `quadrature` | 152 007 | 390 575 | 808 116 | **1 610 052** | 0.768 |

`aggressive` is converged by 2^19 — its **per-realisation** counts are bit-identical
across 2^19–2^21, so the buffer is provably irrelevant there, not merely
indistinguishable. `quadrature` grows by 10.6x over 8x of buffer (span exponents 1.135,
1.022, 0.994) with saturation pinned near 0.77. The overflow ratchet relaxes the cut to
keep the buffer full at whatever size it is given, so there is no escape by spending more.

**Confirmed per level, not inferred from counts** (`ratchet_probe.py`, §6.7). Using the
`threshold_eff` field: `aggressive` ratchets on **0 of 192** levels, `quadrature` on
**60 of 192**, with the effective cut roughly double the nominal. This is the direct
measurement; candidate counts cannot show it (see the memory fact
`saturation-is-not-evidence-of-ratcheting`).

**The cut is not set by the scheme at all** (§6.7). A tenfold reduction in the detection
target buys **4** ratcheted levels out of 192 and leaves zero clean runs: raising the
ladder by 0.87 moved the cut that actually ran by **−0.06**. The operative cut is
`top-K`. No `P_d` target fixes a cut that is not set by `P_d`.

**Buffer pressure is arm-dependent; outcome bias is unresolved** (`paired_max_sugg.py`,
§6.5 — 50 realisations, all four cells, analysis pre-registered before the first cell
ran). Two answers that must not be collapsed into one:

- *Pressure*: at 2^18, **10/50** `quadrature` runs exceed 0.9 saturation against **1/50**
  for `aggressive`, paired McNemar **p = 0.012**.
- *Outcomes*: the pre-registered sign test gives **p = 0.34** — a failure to reject and
  **not** a null — with mean Δ = **+0.080 [−0.043, +0.203]** against a campaign effect of
  interest of **0.06**. The point estimate of the confound exceeds the effect.
- Raising the buffer produced **20 recovery flips, all toward recovery, none away**
  (`aggressive` gained 8 lost 0, p = 0.008; `quadrature` gained 12 lost 0, p < 0.001).
  The ratchet is purely destructive to recovery.

Closing the outcome question needs **n ≈ 211 per cell** — and is unreachable, because it
would have to run at a buffer where one arm still saturates.

**Pairing works, and seeding does not exist** (§11.2, and the memory fact
`wt-injection-design-injections-cannot-be-seeded`). No `seed`/`rng` parameter anywhere in
the chain; `pulse.py:338` is a bare `np.random.default_rng()`, which ignores
`np.random.seed`. Pair instead by persisting `(ts_e, ts_v, dt)` and rebuilding the
`TimeSeries` per arm — **verified**: identical SHA-256 on round-trip, and two processes
running the same arm produced bit-identical sorted score vectors. `DynamicThresholdScheme`
is also unseeded, so ladders must be generated once and committed. A realisation is
**32 MB**, so a campaign must stream rather than persist a corpus.

**`rho_AB` was the wrong quantity** (`rho_check.py`, §5.2). The code used the *profile*
overlap where the model needs the *score* correlation — the same smearing factor weighted
by the pulse's spectrum rather than the filter's, and the selected boxcar has **4.10x**
the pulse's `<k²>`. At ducy 0.10, median `1 − rho` **0.0088 → 0.0265**; `pi`
**0.697 → 0.628**. Still uncorrected at ducy 0.05 and 0.20.

**The stratification inherited that error** (§5.1). Corrected, the **`effect` stratum does
not exist** — no sampled position reaches `pi` > 0.70 — so §9's three strata collapse to
two, and the `null` stratum, the design's main control, roughly doubles to about **1
position in 5**.

**Cost, measured** (§11.3). Per-pair wall clock at 2^18 is **49.8 s** (`aggressive` 2.3 s,
`quadrature` 47.5 s), median over runs 2..n to exclude ~25 s of one-off numba JIT.
`quadrature`'s cost scaled ~13x for a 16x buffer increase, roughly linearly.

## What shipped, and what is held

- **Shipped.** The silent-ratchet defect found along the way is
  [PR #14](https://github.com/pravirkr/pyloki/pull/14), **merged upstream** — it adds the
  `threshold_eff` field that made every measurement above possible. Implemented on a
  separate branch off `upstream/main`, not here; verified in both directions by a
  regression test. Write-up: `06_upstream_max_sugg_logging.md`.
- **Held, unposted.** `07_upstream_buffer_policy.md` — the buffer-policy question, framed
  as an *issue* rather than a patch because the fix is a design call. The branch is pushed
  to the `assaferan` fork **only** so that draft's reproducer links resolve. Nothing has
  been posted to the upstream tracker, and nothing should be without assaferan reading it
  first.
- **Deliberately not shipped.** Anything about tiling. Per §8 nothing in this document
  supports an upstream sentence about it, and the buffer result is stronger standing
  alone. If `metric-gridding`'s `04_upstream_report.md` and this ever both go out, they
  must be consistent with each other.

**No inherited library code was modified by the design work on this branch.**

## Withdrawn — do not cite these

| withdrawn | replaced by |
|---|---|
| §4.3, "the decision is made where the geometry isn't" | **Reversed.** The measured survival profile puts first loss at level 13 and none before it, against a modelled 99% of losses by stage 10. The per-stage `pi` row survives — it is geometry — but the weighting over it does not |
| §10 row 8, the equal-`P_d` ladder check | Void. Batch A's 3/24 agreement ran at 2^14 where **both** arms saturate, so it cannot be credited to the ladders. The register now has **no** independently verified row |
| The "defensible upstream sentence" about tiling (§8) | Nothing in this document currently supports one |
| `n` = 954, and `n` = 368 | **368–691** (§9 caution 3). Both are marked in place |
| `44/70 = 0.629` and its `p` = 0.164 | **29/40 = 0.725 vs 38/50 = 0.760, p = 0.81.** The pooled figure added batch 1's *excursion* numerator to batch 2's *metric* numerator and silently dropped a third group; all three would be 61/100. Do not re-test it |
| "`m <= 1.0` is stricter than end-to-end recovery, so a lower rate is expected" | No ordering exists. Different norms of the same `delta`, independently normalised; and recovery is evaluated after the final ascend and resolve, so neither is nested in the other |

## Not established, and why

- **Whether the geometry converts into detections.** The question the branch existed to
  answer. Blocked as above, at this configuration, permanently.
- **`n`, precisely.** It is **368–691**; the **368** end is the excursion criterion
  (median loss level 27), the **691** end the metric `m <= 1.0` criterion at its steepest
  (median 15). The early-decision corner (2 842 pairs) is excluded. The spread is the
  survival profile's own, and the two criteria do not agree on where the curve falls.
- **The survival profile itself.** Inherited from `metric-gridding`, **uncontradicted but
  not validated**. Its *shape* — no losses before level 13 — is robust across two
  criteria, two parameter sets and three groups spanning two batches (N = 100). Its
  *absolute levels* carry nothing. Quote the shape; name the criterion before quoting any
  fraction or median from it.
- **ducy 0.05 and 0.20**, still on the uncorrected `rho_AB`.

## If this is reopened, in order

1. **`DECISION_STAGES = [2, 6, 10, 14, 20, 28]` is misplaced.** It sits almost entirely
   *before* the measured decision window (first loss at level 13, median 27), so §5.1's
   stratification and §9's strata — including the null stratum the design leans on — are
   built over stages where little is being decided. Not recomputed here, because that
   would put a new headline on an inherited, unverified measurement. **First thing to do.**
2. **Replicate the survival profile here** to settle the `n` range. Note it can only be
   replicated honestly in the `aggressive` arm; a "both arms" replication would measure
   the ratchet in one of them.
3. **Correct `rho_AB` at ducy 0.05 and 0.20.** Leading order says 0.05 barely moves and
   0.20 moves more than 0.10 did.
4. **§10 row 2** is retired for stages 1–10 on the *other* branch's measurement. If it has
   to carry weight, re-run it here.

The three routes that remain open at all are in §1: ask the user-facing question instead
(answerable now, but it is a claim about pyloki-as-shipped, not about tiling), test a
smaller configuration (a regime nobody deploys), or change the buffer policy upstream
(`07`, the only remaining route to the idealised question at the deployed configuration).

## Running anything here

**`import pyloki` in a worktree resolves to the MAIN checkout's `src` unless you set
`PYTHONPATH`.** The shared venv has `pyloki` installed editable against `main`. There is
no bare `python` on PATH, and no `uv` or `timeout`.

    PYTHONPATH=$PWD/src /Users/assaferan/Documents/GitHub/pyloki/.venv/bin/python <script>

For anything that spawns subprocesses, **assert** the resolution rather than trusting the
variable — `paired_max_sugg.py:_env()` runs `python -c "import pyloki; print(pyloki.__file__)"`
and refuses to start if the path is not under the worktree's `src`. Copy that guard. A
module present on both branches with different contents is picked up **silently**, and a
campaign run under the wrong library looks exactly like a campaign.

**The cached schemes are not interchangeable.** `schemes/aggressive.npz`,
`conservative.npz` and `metric.npz` are **Taylor** (`kind="poly_taylor_moving"`) and there
is no `quadrature` entry among them. Using them for the Chebyshev arm silently compares
the wrong ladders. The correct ones are `schemes/cheby_*.npz`.

## Methodological conclusions worth outliving this branch

1. **Measure the mechanism, not a proxy for it.** `ncand/max_sugg` is read at the end of a
   level; the ratchet fires *during* one and never relaxes. Two sections were written
   leaning on saturation before the direct measurement existed, and both had to be
   corrected. The direct measure separated the arms cleanly where saturation had not.
2. **A number can be correct in its own frame and wrong in the frame it is quoted in.**
   Every defect this branch produced is that one move: a fraction whose denominator
   counted something else, a one-sample test against another session's *estimate*, a
   pooled rate mixing two predicates, a cross-reference into a section that no longer
   exists.
3. **Quote no fraction without its denominator and its predicate.** Not for tidiness — it
   is what makes the next cross-check possible. A fact reading "two batches" beside three
   figures was visibly inconsistent but not *diagnosable*; a reader could tell something
   was wrong and not what.
4. **A predicate evaluated at one point in a pipeline is not nested in one evaluated after
   two more stages of it**, however intuitive the ordering feels. Neither direction bounds
   the other unless you have shown it.
5. **A correction is itself an unreviewed edit**, written with more confidence than the
   text around it and given less scrutiny. Both corrections on 2026-09-22 introduced a new
   defect while fixing an old one. Read the table, then the sentence above it — and the
   sentence below it.
6. **Void at the granularity of the datum, not the result.** Over-withdrawing has its own
   cost: the conclusion behind the withdrawn `44/70` survived intact on the comparable arm
   alone, and retracting it wholesale would have destroyed a valid constraint.
7. **Failure to reject is not a null.** The one that cancelled the campaign: p = 0.34 with
   a point estimate above the effect of interest is not permission to proceed.

## Files

    05_injection_design.md          the design, the verdict, the handoff — START HERE (§1, §6.6, §11)
    06_upstream_max_sugg_logging.md the logging patch — SHIPPED as PR #14
    07_upstream_buffer_policy.md    the buffer-policy issue — DRAFT, UNPOSTED, held for review

    saturation_sweep.py             buffer sweep; the non-convergence result (§6.6)
    ratchet_probe.py                per-level threshold_eff; the direct ratchet measure (§6.7)
    paired_max_sugg.py              the pre-registered paired experiment (§6.5); _env() guard
    injection_pilot.py              realisation generation, persistence, per-run timing
    injection_power.py              the power calculation; DECISION_STAGES lives here (see item 1)
    rho_check.py                    score correlation vs profile overlap (§5.2)
    schemes/cheby_*.npz             the Chebyshev ladders — the only valid ones for this arm
    schemes/{aggressive,conservative,metric}.npz   TAYLOR — invalid here, see above

    DECISIONS.md                    Phase 0-3 log, ends 2026-09-17; does NOT cover this work
    03_results.md, 04_upstream_report.md           inherited; 04 is owned by metric-gridding
