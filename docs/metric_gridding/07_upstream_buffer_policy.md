# A high-branching search cannot run at its configured threshold scheme

**Status: draft for review. Not posted, not pushed.** Proposed as an *issue* — a
question about pruning semantics — rather than a patch, because the fix is a design
call. Successor to [PR #14](https://github.com/pravirkr/pyloki/pull/14), which made the
behaviour below visible; this asks whether it should be avoidable.

*All line numbers refer to `upstream/main` at `18d04b3`. All measurements use the
`threshold_eff` field added by PR #14.*

## Summary

When the candidate buffer overflows, pruning ratchets the effective cut to
`max(nominal, top-K, median)` and never lowers it again for that level
(`utils/world_tree.py:547`). For a configuration whose branching pattern is large
relative to `max_sugg`, **this is not an occasional event that a bigger buffer avoids —
it is the steady state**, and the `threshold_scheme` the user configured is never the
operative cut.

Recalibrating the scheme does not help, because the scheme is not what is cutting.

This is a **reporting-and-policy** issue, not a correctness bug. Ratcheting to stay
inside a fixed buffer is a reasonable way to bound memory. The question is whether a
user should be able to choose to fail instead.

## Measurements

One search configuration — 64 segments, `poly_order = 4`, `branch_max = 16`, Chebyshev
basis, `nbins = 64` — differing only in `tiling_strategy`, which changes the branching
pattern product:

| strategy | product of branching pattern |
|---|---|
| low-branching | 1.5 × 10¹² |
| high-branching | 5.4 × 10²⁰ |

### 1. The buffer is the steady state for the high-branching case

Candidate count and saturation on **one fixed set of 10 noise realisations**, swept
across four buffers so every row is the same data:

| strategy | 2^18 | 2^19 | 2^20 | 2^21 |
|---|---|---|---|---|
| low-branching, median `ncand` | 16 778.5 | **16 778.5** | **16 778.5** | **16 778.5** |
| low-branching, median saturation | 0.064 | 0.032 | 0.016 | 0.008 |
| high-branching, median `ncand` | 152 007 | 390 575 | 808 116 | **1 610 052** |
| high-branching, median saturation | 0.580 | 0.745 | 0.771 | 0.768 |

(The 2^18 column comes from a 50-realisation run restricted to the same 10; 2^19–2^21
are from the sweep itself. Same realisations throughout.)

The low-branching search converges: its **per-realisation** counts are bit-identical
across 2^19–2^21, so the buffer is provably irrelevant, not merely indistinguishable.
The high-branching one never converges — over a factor of 8 in buffer its count grows by
a factor of 10.6, span exponents `d log2(ncand) / d log2(max_sugg)` = 1.135 / 1.022 /
0.994, with saturation pinned near 0.77. Doubling the buffer doubles the count.

### 2. The scheme is not the operative cut

At `max_sugg` = 2^20, three fixed realisations, 64 prune levels each, counting levels
where `threshold_eff > threshold`:

| strategy | ladder `P_d` | thresholds | levels ratcheted | runs that ran the configured scheme |
|---|---|---|---|---|
| low-branching | 0.1031 | 2.10–7.70 | **0 / 192** | **3 / 3** |
| high-branching | 0.1031 | 2.00–7.90 | 60 / 192 (20.0 per run) | **0 / 3** |
| high-branching | 0.1083 | 1.40–8.10 | 60 / 192 (20.0 per run) | **0 / 3** |
| high-branching | **0.0100** | 2.40–8.50 | 56 / 192 (18.7 per run) | **0 / 3** |

The low-branching row is the control, and it is perfectly clean. The high-branching
search ratchets on roughly 30% of levels, with the realised cut **about double** what
was configured: max excess over nominal 4.19 / 4.98 / 5.32 on thresholds of order 2–8.

### 3. Recalibrating the scheme is not the lever

The last two high-branching rows come from the **same** `DynamicThresholdScheme` run,
backtracked at different `P_d`, so they differ only in the target. A tenfold reduction
in the detection target moves the ratcheted count from 20.0 to 18.7 per run and leaves
**zero** clean runs. The candidate count does not fall either — median 730 462 against
949 921, i.e. it rises, because the count is pinned to the buffer.

Per level, on one realisation under both ladders:

| level | nominal (base) | nominal (strict) | Δ nominal | effective (base) | effective (strict) | **Δ effective** |
|---|---|---|---|---|---|---|
| 4 | 1.80 | 3.30 | +1.50 | 2.93 | 3.30 | +0.37 |
| 5 | 2.30 | 3.50 | +1.20 | 4.35 | 4.40 | +0.05 |
| 7 | 2.90 | 4.20 | +1.30 | 5.39 | 5.38 | −0.01 |
| 10 | 3.90 | 5.00 | +1.10 | 5.92 | 5.92 | **0.00** |
| 13 | 4.60 | 5.80 | +1.20 | 6.36 | 6.37 | +0.01 |
| 19 | 6.00 | 6.70 | +0.70 | 7.29 | 7.31 | +0.02 |

Over all 20 ratcheted levels:

    mean change in NOMINAL threshold (strict − base):  +0.870   sd 0.344
    mean change in EFFECTIVE cut                    :  −0.063   sd 0.175
    |Δ effective| / |Δ nominal| = 0.07

**Raising the ladder by 0.87 moved the cut that actually ran by −0.06.** At level 10 the
effective cut is identical to two decimal places while the nominal rose by 1.10.

## Mechanism — `max(nominal, top-K, median)` behaving exactly as written

The effective cut is `max(nominal, top-K, median)` where `top-K` is the score of the
`max_sugg`-th best candidate. Once the buffer fills, `top-K` exceeds the nominal
threshold and sets the cut; it is a function of the buffer size and the score
distribution, and the configured scheme drops out of the maximum entirely.

The one level the strict ladder recovered confirms this rather than complicating it.
At level 4 the strict nominal (3.30) rises **above** the baseline's effective cut
(2.93), so the nominal becomes binding again and `eff == thresh` exactly. Recovering all
20 levels would require nominal > top-K at *every* level — thresholds of 6–8 from the
first stage, which is not a stricter operating point but no search at all.

## The question

`max_sugg` is currently a hard cap enforced by discarding the lowest-scoring candidates,
and the discarding is silent in its effect on the threshold semantics (PR #14 makes it
visible in the log, but it still happens). Would it be reasonable to let a user choose a
different policy — for example:

- **fail fast**: raise when the realised cut departs from the scheme by more than a
  configured tolerance, so a configuration that cannot run at its own thresholds says so
  rather than silently running at different ones;
- **or report the shortfall** as a first-class output — "this run applied a cut of 7.3
  where the scheme asked for 2.5 at level 10" — so a user comparing configurations knows
  which comparison they actually made.

Either would let a high-branching configuration fail honestly instead of quietly
becoming a different search. We have no view on which is right; it touches pruning
semantics and is a design call for the maintainer.

## Scope — what this does NOT claim

- It does **not** claim the ratchet is a correctness bug. Bounding memory by discarding
  the weakest candidates is reasonable, and the search still recovers signals: every run
  in the table above recovered the injected signal.
- It does **not** claim the ratchet favours one `tiling_strategy` over another in a way
  that biases a comparison between them. That was measured separately and is
  **unresolved** — a paired test over 50 realisations gave a two-sided p of 0.34, with
  the point estimate large enough that we would not rely on it either way. The case here
  stands without it.
- It does **not** depend on any result about tiling, grid geometry, or template
  placement. The two strategies appear only as a convenient way to vary the branching
  pattern within one configuration.

## Reproducing

All of this is single-threaded and needs no special data.

- The convergence sweep and its growth exponents: `saturation_sweep.py`, results in
  `saturation_sweep_results.json`.
- The ratchet counts, the stricter-ladder comparison and the per-level table:
  `ratchet_probe.py`, results in `ratchet_probe_results.json`. `--ladders` regenerates
  the matched `P_d` ladders from one scheme run (`DynamicThresholdScheme.__init__` is
  unseeded, so ladders from separate runs are not comparable).

Every number in this document has been checked against those two files
programmatically rather than by eye (32 claims, all passing).

**Limitations, stated plainly.** The ratchet counts are three realisations per cell at
one buffer, and the per-level table is one realisation. The effect is not marginal —
0/192 against 60/192, with 20/20/20 consistency run to run — so the separation is not in
doubt, but the per-run ratcheted count (18.7 against 20.0) is not resolved at this n and
nothing is claimed from that difference beyond "not zero". One injected signal, at one
S/N, in one configuration family.
