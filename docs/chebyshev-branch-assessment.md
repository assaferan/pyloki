# Chebyshev branch: grid closed forms, and what is publishable

Working note. Not part of the package. Written 2026-09-16 after rebasing the
branch onto `upstream/main` (9f8101f).

## 1. Corner-error closed form in both bases

Conventions, stated once because the branch history is ambiguous about them:

- `k_max == poly_order == ` the number of **branched** coefficients. In both bases
  the constant term is carried but never branched (`poly_chebyshev_seed`: "we never
  branch on d0"), so the branched set is orders `j = 1..k_max`.
- `j` below is the **polynomial order**, 1-based. `psr_utils.poly_taylor_step_f`
  uses `k = np.arange(nparams)`, a 0-based *array index*, with `fact(k+1)` and a
  `(k+1)` power; its `2**k` coarsening is therefore `2**(j-1)` in order terms.
  Both spellings name the same factor.
- "Corner error" is the sup-norm phase distance from a cell centre to the cell
  corner, every branched coefficient displaced by half a grid step. This is the
  quantity `scratch/econ_experiment/grid_guarantee.py` tabulates, in units of the
  nominal tolerance `eta/N_b`.

Per branched order `j`, the corner contribution is

| basis | per-order | total over `j = 1..k_max` |
|---|---|---|
| Taylor, `use_cheby=True` (shipped) | `2**(j-2) * eta/N_b` | `(2**(k_max-1) - 1/2) * eta/N_b` |
| Taylor, `use_cheby=False` | `1/2 * eta/N_b` | `(k_max/2) * eta/N_b` |
| Chebyshev | `1/2 * eta/N_b` | `(k_max/2) * eta/N_b` |

Verified against `poly_cheb_step_vec` and `poly_taylor_step_d_vec` for
`poly_order` 2..8 at machine precision, and for Chebyshev shown independent of
`nbins`, `eta`, `f_max` and `t_s` (worst relative deviation 1.3e-14 over a grid of
all four). Reproducer: `scratch/econ_experiment/grid_corner_closed_form.py`.
The Taylor closed form is not this branch's result; it was established
independently on `metric-gridding` and is reproduced here only as the control.

`k_max/2` reproduces the 2.0 and 2.5 that `grid_guarantee.py` already reported for
Chebyshev at `poly_order` 4 and 5.

### What that means for the "linear vs geometric" conjecture

The conjecture — that `|T_k| <= 1` makes the Chebyshev corner sum linear in
`k_max` rather than geometric — is **confirmed in form and wrong in constant and
in mechanism**.

- Constant: it is `k_max/2`, not `k_max`. The factor 2 is the half-step, which the
  Taylor form `(2**(k_max-1) - 1/2)` already carries (`2**(j-2)` is half of
  `2**(j-1)`) and which the conjecture dropped on the Chebyshev side.
- Mechanism: `|T_k| <= 1` is not what produces the linear behaviour. The Taylor
  grid with the coarsening removed is *also* exactly `1/2` per order, because
  `poly_taylor_step_f` sets `dparams_f[k] = dphi * (k+1)!/T**(k+1)`, which is
  calibrated so that a half-step in `d_j` alone contributes exactly `dphi/2` of
  phase at the far end of the span. Both bases are "one axis = half the
  tolerance". The **entire** Taylor/Chebyshev corner gap is the `2**k` factor.

So the corner closed form is not an argument for the Chebyshev basis. It is a
quantification of the `use_cheby` coarsening factor, and it says nothing about the
basis that `use_cheby=False` would not also say. Commit 36e460e ("Taylor with the
factor removed lands exactly on the Chebyshev values") already observed this
numerically at two orders; the closed form makes it an identity at all orders.

## 2. The corner error is not a coverage deficit — and this contradicts 36e460e

36e460e reads the 7.5x/15.5x corner numbers as the Taylor grid failing to honour
its advertised tolerance, and concludes "the Chebyshev basis explores more leaves
because it keeps the stated tolerance, not because it is less efficient". That
inference does not survive checking, in either direction.

The corner error measures the distance from a point to **the centre of the cell
containing it**. A grid search does not need that point: it needs *some* evaluated
grid point within tolerance, and the branching tiles the space contiguously
(`branch_param`: "zero overlap and zero gaps"), so every neighbour exists. The
operative quantity is the **covering radius**: the max over points in a cell of
the distance to the *nearest* grid point. Measured (random points in a cell,
nearest lattice point over integer step offsets; section B of the reproducer):

| `k_max` | Taylor `2**k` | Taylor raw | Chebyshev |
|---|---|---|---|
| 2 | 0.50 | 0.50 | 0.99 |
| 3 | 0.50 | 0.50 | 1.00 |
| 4 | 0.50 | 0.50 | 1.19 |
| 5 | 0.50 | 0.50 | 1.31 |

(plain sup-norm, multiples of `eta/N_b`. Absorbing the free constant term — which
is legitimate, since `d_0`/`alpha_0` are never branched and a constant phase
offset only rotates the folded profile — gives 0.27 / 0.25 / 0.77-1.19 with the
same ordering.)

Three consequences, all of which cut against 36e460e:

1. The shipped Taylor grid covers at 0.5 `eta/N_b` at every order tested, not 7.5
   or 15.5. It does not under-cover. The `2**k` factor is doing exactly the job
   its name implies: it is the Chebyshev economization factor, coarsening the
   high orders by the amount that the low orders can absorb.
2. The `2**k` factor is close to free. Dropping it takes the covering radius from
   0.50 to 0.50 while multiplying the cell count by `2**(k_max(k_max-1)/2)` — 64x
   at `poly_order=4`, 1024x at 5.
3. Chebyshev is the **looser** of the two grids, and is the only one that exceeds
   its own nominal tolerance (1.19 at `poly_order=4`, 1.31 at 5).

Point 3 is consistent with a volume computation and inconsistent with the branch's
leaf-count story. Mapping both cells into a common coordinate system with the
codebase's own `cheby_to_taylor`, the Chebyshev cell is `2**(k_max(k_max+1)/2)`
times **larger** than the shipped Taylor cell — 1024x at `poly_order=4`,
independent of `T`. A coarser grid cannot be producing more leaves because it is
finer. The extra Chebyshev leaves must come from the error-propagation path
(`shift_cheby_full` mixing orders inflates the carried `dparam`, so each leaf
branches more), which is precisely the mechanism 00ff1b4 identified for the
`tiling_strategy` effect. 00ff1b4 already supersedes 36e460e's conclusion; this
supplies the arithmetic for why.

**Nothing here touches 00ff1b4.** That result is a measured paired experiment on
EP score and does not depend on the grid-guarantee analysis. It stands.

Caveats on the covering-radius numbers: they are an empirical max over random
points with the lattice search truncated to a few steps per axis, so each is an
upper bound on the true covering radius — which only strengthens the Taylor side.
And the model is idealised: the pruning tree can remove the covering neighbour at
an earlier stage, and per-leaf propagated `dparam` makes the real grid only
approximately a lattice. Settling those is the obvious next step, and until they
are settled the 7.5x/15.5x under-coverage claim should not go upstream.

## 3. Rebase

Rebased onto `upstream/main` (9f8101f) with no conflicts. The merge commit
7aae975 (`fix-prune-segfault`) was dropped by the rebase because its content is
already in `upstream/main` via 2b4b80c; 11 non-merge commits replayed cleanly, and
the branch now has a single merge base. The only tree change relative to the
pre-rebase tip is upstream's one-line notebook fix (c181aea). A safety ref
`backup-pre-rebase-Chebyshev` points at the old tip.

`pytest`: 35 passed, 1 warning (a pre-existing `NumbaPendingDeprecationWarning`
about reflected lists in `generate_bp_poly_taylor`). Everything the branch added
still passes.

## 4. What is publishable

`git diff --stat upstream/main...HEAD` is 2974 insertions across 37 files. It
decomposes as:

| part | lines | verdict |
|---|---|---|
| `src/pyloki/utils/transforms.py` | +67 | the whole of the shippable change |
| `tests/test_transforms.py` | +90 | ship with it |
| `.gitignore` | +1/-1 | trailing newline only; drop or keep, immaterial |
| `HANDOFF.md` | +624 | session note, do not ship |
| `scratch/econ_experiment/` | +2192, 33 files | do not ship |

The `scratch/` tree carries `.npy` and `.jsonl` data files and scripts that hard-code
this fork's experiment layout. Its own README says "nothing here is part of the
package". None of it belongs in a PR.

So the minimal defensible contribution is **one new function plus its tests**:
`transforms.economize_taylor_params`, the near-minimax degree-reduction of a
Taylor expansion via Chebyshev economization, with six tests covering the no-op
case, agreement with naive truncation when there is no high-order content, the
error reduction, the round trip, and the parity limitation.

Two things a reviewer will raise, and they should be answered before the PR, not
during it:

1. **It has no caller.** e7d51b8 introduced it in the `resolve` step and e561317
   took it back out again, so on the current tip `economize_taylor_params` is
   reachable only from tests and from `scratch/`. Upstream is being asked to carry
   a public utility that the package does not use. Either land it together with a
   caller, or propose it explicitly as a utility and say so in the PR body.
2. **Its documented value is narrow.** The docstring's own Notes say that at
   `poly_order=3` — the regime where the truncation error is largest — the parity
   structure means economization changes nothing that matters, and the ~0.244x
   sup-norm improvement only appears at `poly_order>=4`. That is honest and it is
   the right thing to have written down, but it is also the reason the change was
   reverted from `resolve`, and it should be in the PR description rather than
   left for a reviewer to discover in the Notes.

The grid analysis in sections 1 and 2 is a separate matter from this PR. It is a
result about `poly_taylor_step_f` and `poly_cheb_step_vec`, it changes no code, and
section 2 contradicts a conclusion this branch itself recorded. It should be
settled — specifically, the pruning caveats above — before any of it is written up.

### Branch shape

`git diff main...Chebyshev` warned about multiple merge bases before the rebase:
the branch forked at 6b11aba, then merged 2634e81 (`fix-prune-segfault`) at
7aae975, and 2634e81 later reached `main` independently. After the rebase there is
a single merge base and the warning is gone. Any diff taken against the *old* tip
should be re-taken.
