# metric-gridding — what this branch established, and what it withdrew

Entry point. `DECISIONS.md` is 2 600+ lines over 41 sessions with 116 numbered
decisions and roughly a dozen retractions; read this first or you will cite something
that was withdrawn. Every claim below points at the decision that carries it.

## How it ended

The branch set out to replace the axis-aligned box grid with a metric-based covering,
following the closing remark of Kumar & Zackay (2026) §5.2.4. That premise did not
survive: the shipped box grid is already metric-derived (D38), the box strategies leave
no coverage gap in the Taylor basis (D36), and a metric covering costs ~1070x at equal
`P_d` for a guarantee the box grid already provides (D35).

What replaced it is a sharper question and a negative answer to it:

> **`quadrature`'s geometric advantage over `aggressive` is exactly established — 89/90
> cells, `>= 2.30x` in Taylor and `>= 1.56x` in Chebyshev, 100% of cells in the measured
> decision window — and is worth about 0.05 in score, against a configuration that
> cannot run at its own threshold scheme at any affordable buffer, and no recalibration
> of that scheme changes it.** (D52, D60, D106, D114-D116)

The advantage is real, exactly quantified, and swamped. Nothing failed in the geometry;
what failed is that the configuration the geometry favours is precisely the one that
cannot run at its own calibration.

## Live results

**Exact nearest-template search** (`nearest_template.py`, `nearest_template_cheby.py`,
`tests/test_nearest_template.py`). Child offsets from `branch_param_padded` are
parent-independent, so the leaf set at a stage is *exactly* a Minkowski sum of small
per-stage offset sets — 1e34 leaves in a few hundred points per stage. Branch and bound
with an admissible bound; agrees with brute-force enumeration to 1e-12 on every
enumerable case in both bases. Returns `(value, exact)` and never passes a truncated
search off as a minimum. (D50, D58)

- `aggressive`'s nearest template: **1.069** `eta/N_b` (Taylor), **0.904** (Chebyshev),
  exact in all 90 cells. (D51, D60)
- `quadrature` strictly closer in **89/90** cells (Taylor), **87/90** (Chebyshev). (D52, D60)
- In the *measured* decision window: **2.1-2.5x** (Taylor), **~1.5x** (Chebyshev),
  closer in **100%** of cells. (D106)
- Threshold-independent by construction — no threshold, score or buffer quantity enters
  it, and it enumerates the full tree rather than the surviving set. (D114)

**Grid closed forms.** The corner of the shipped optimal grid costs exactly
`(2^(k_max-1) - 1/2) * eta/N_b`, verified for orders 2-8, independent of `t_s` and
`f_max` when evaluated at the span that enters the step formula. 7.5x at `poly_order=4`.
The Chebyshev analogue is `k_max/2` — linear rather than geometric. (D46, D59)

**But the corner is not a covering radius** (D48): a search enumerates every grid point,
and minimising over neighbours gives ~0.5. `eta` is *not* optimistic by 7.5x, and the
`2^(k-1)` coarsening is economization working rather than a defect. The nominal corner
also predicts the amplitude ordering **wrongly** — Chebyshev is 3.75x better on it and
slightly worse in S/N. (D61)

**Amplitude conversion** (`amplitude_loss.py`). No metric: the folded profile is
`P(k) * S(k)` with `S(k)` the characteristic function of the phase residual, scored with
the shipped boxcar bank and cross-checked against a matched filter. At 10% duty,
`aggressive` loses **1.73%** and `quadrature` **0.75%** (Taylor). Sup-norm phase error
*overestimates* the loss by ~3x, because the residual attains its peak only briefly.
(D55, D56)

**Measured per-stage survival** (`survival_profile.py`). Drives the real pruning loop and
records the nearest survivor to the truth at every level. **No losses before level 13**,
across two criteria, two parameter sets and two batches — where a single-leaf model put
99% of losses in stages 1-10. The decision is made in the middle and late stages. (D101,
D110)

**`branch_max` at the shipped default of 16**: `aggressive` needs 7 (Taylor) and 9
(Chebyshev); `quadrature` needs 28 and 16; `conservative` 27 and 20. The guard is strict,
so **Chebyshev + `quadrature` is the only non-`aggressive` configuration that builds**,
and it does so with no headroom. (D62, amended)

**Report figures are generated, not typed** (`report_numbers.py`,
`tests/test_report_numbers.py`): the report is bound to computed values, including the
`branch_max` verdicts as booleans. It caught two defects on its first run. (D63-D65)

## Withdrawn — do not cite these

| withdrawn | replaced by |
|---|---|
| D39's conclusion, "transport can buy cost and never sensitivity" | D52. Its *mechanism* (the cell is criterion-set) is correct and survives as **D41** |
| D43's magnitude, "3.1x" | D52's `>= 2.30x`. The method was unsound (**D47**: a keep-best-N cap, minima not monotone in search width) |
| D44, "0.34% in S/N" | D55 |
| D45's sensitivity half, "`conservative` is dominated" | withdrawn; its **cost** side (2^15) stands |
| D46's *interpretation*, "`eta` is optimistic by 7.5x" | D48. The arithmetic stands |
| D70's framing, "the gain sits where detection is not decided" | **reversed** by D101/D102. The early/late split itself stands as geometry |
| D68/D71's pair counts (~290, ~770) | the survival model behind them is **retired** (D95), falsified against a converged measurement (D94) |
| D77, the pilot "validation" | void — it measured the ratchet, not the ladder (D92) |
| D78, "inject at S/N 15-17" | D86. Fitted an offset to one point, which cannot separate offset from shape (D85) |
| D109's "discrepant" label | **D111**: the correct test is two-sample, `p = 0.16`, consistent. The profile is *unvalidated*, not contradicted |
| `pruning_multiplicity.py` | `nearest_template.py`. The module is marked superseded in its own docstring |

## Not established, and why

- **Whether the geometry converts into detections.** `quadrature`'s candidate count is
  proportional to `max_sugg` over a factor of 8 with no convergence, so `prune_on_overload`
  is always active for it and the scheme's thresholds are never the operative cut. There
  is no affordable buffer at which both arms run at their own calibration. (D91)
- **The profile's absolute level.** Consistent with an external recovery rate but not
  validated by it; one batch agreeing was never evidence. The *shape* carries the
  conclusions. (D109, D111)
- **A `quadrature` survival profile**, blocked by the same non-convergence (302 s/run at
  2^21).
- **The circular basis.** Untested throughout.

## Methodological conclusions worth outliving this branch

1. **Fix an operating point by the quantity you can measure, not the one you have to
   model.** Two sessions appeared to disagree about where to inject; the criteria turned
   out to be identical in different variables, and the one parameterised in *measured*
   `P_d` was immune to a 2x error in the other's `S/N -> P_d` mapping. (D89)
2. **A number can be correct in its own frame and wrong in the frame it is presented in.**
   Four instances here: a table with two denominators, a one-sample test against someone
   else's estimate, line-number citations against the wrong branch, and a model right
   about the `P_d` at the optimum and wrong about which S/N delivers it. (D107, D111, D89)
3. **Bound, don't cap.** A one-sided bound is harmless when comparing two numbers and
   fatal inside a model that takes an extremum over many of them — a loose per-stage bound
   inverted a whole power calculation. Replace truncation with an admissible bound, and
   report `(value, exact)`. (D47, D66)
4. **Void at the granularity of the datum, not the result.** Over-withdrawing has a cost:
   retracting a comparison wholesale left a valid constraint unused, and finding it later
   falsified a model outright. (D96)
5. **One batch agreeing is not a validation.** Run the second batch before calling it one.
   And when comparing against someone else's measured rate, test two-sample — their
   number has an interval too. (D109, D111)
6. **Check what the comparison datum measures before interpreting the residual.** Three
   successive characterisations of one model's error were all computed against data
   incapable of measuring it; one valid datum settled it immediately. (D94)
7. **What smears a profile is the residual's *distribution*, not its peak.** This flipped
   a conclusion three times. (D56, D61)
8. **When asked whether a result depends on X, grep for X** rather than reasoning about
   whether it should. (D114)

## Files

    nearest_template.py       exact search, Taylor            tests/test_nearest_template.py
    nearest_template_cheby.py exact search, Chebyshev         (same)
    amplitude_loss.py         phase residual -> S/N loss      (same)
    survival_profile.py       measured per-stage survival
    report_numbers.py         every figure in the report      tests/test_report_numbers.py
    sensitivity_loss.py       own-cell excursions; see D43/D47 before quoting
    pruning_multiplicity.py   SUPERSEDED AND UNSOUND (D47)
    04_upstream_report.md     draft, unposted
    05_injection_design.md    pre-registered design; campaign not run
    03_results.md             Phase 3; carries a banner, partially superseded
    DECISIONS.md              full log, 41 sessions
