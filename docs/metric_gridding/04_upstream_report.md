# 04_upstream_report.md — draft report for pravirkr/pyloki

Status: **draft, not yet posted.** Upstream has Discussions disabled, so this would be an
issue. Rewritten 2026-09-16 around D43-D46; the previous headline ("`tiling_strategy`
buys no sensitivity") was refuted by D43 and must not be posted. Numbers verified in
session (q).

---

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
