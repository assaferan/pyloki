# Questions about the paper (arXiv:2607.07700)

Running list of questions to revisit later, cross-referenced against the pyloki codebase.

## Section 3.3 — Phase-Coherent Folding

1. The time-domain folded profile is defined by binning on `Phi(t_n)` - the instantaneous
   rotational phase. It says it is derived from the ephemeris. Do we derive it from the
   time series `T_n` in some direct manner? Or is it an unknown quantity depending on our
   estimate of `Phi(t)`?

   > **Answered (2026-09-08):** `Phi(t_n)` is computed purely from trial/hypothesized
   > phase-model parameters (a grid point in the search space), applied to sample times
   > `t_n = n*t_s`. It is not derived from or adapted to `T_n` — see
   > `src/pyloki/utils/psr_utils.py:13-56` (`get_phase_idx`/`get_phase_idx_int`) and
   > `src/pyloki/core/fold.py` (`brutefold_single`, `brutefold`, `brutefold_start`), where
   > `T_n` only enters as the amplitude accumulated into each bin.

## Section 3.2 — Polynomial phase model (and appendices on Chebyshev/orthogonal bases)

1. The way we use to approximate a function in the search space is via its derivatives at
   points (the Taylor approximation). Is this choice optimal? What are we trying to
   minimize? Is there a better approximation (e.g. Fourier/Chebyshev) that would achieve
   better results?

   > **Answered (2026-09-08, revised again — treating `Phi(t)` as a genuinely unknown
   > smooth function, not an assumed Taylor series):**
   >
   > **Taylor and Chebyshev, truncated at the same order `k_max`, are not competing
   > families.** A degree-≤`k_max` polynomial in the monomial basis `(t-t_c)^k` and the
   > same polynomial in the Chebyshev basis `T_k(x)` are the identical `(k_max+1)`-dim
   > vector space, just different coordinates for it. So switching basis doesn't change
   > which `Φ(t)` the grid can reach — only how efficiently (how few points) it can
   > guarantee reaching them, since monomial coordinates are badly conditioned/correlated
   > while Chebyshev coordinates are nearly diagonal. That's what Appendix D buys: a
   > grid-*efficiency* gain, not a modeling gain. The real question is why a degree-`k_max`
   > polynomial family at all (in any basis), for approximating an unknown smooth `Φ(t)`
   > on `[0, Tobs]`.
   >
   > **That does have a real (non-assumed) answer from approximation theory:** for a
   > smooth/analytic function on a finite, *non-periodic* interval, polynomial
   > approximation is close to the best possible finite-parameter scheme (Jackson-type
   > theorems: best degree-`n` polynomial error decays geometrically in `n` for analytic
   > functions). Truncated Fourier/trigonometric approximation is optimal only for
   > genuinely *periodic* functions; on a finite non-periodic window it suffers
   > Gibbs-type boundary artifacts and only converges algebraically — i.e. it's worse for
   > a generic smooth function on a bounded interval. (Chebyshev approximation is in fact
   > equivalent to a cosine-Fourier series under `t = cos θ`, which handles the
   > non-periodic boundary correctly by reflection — "Fourier done right for a finite
   > window.") So polynomial-degree-`k_max` is close to optimal independent of any
   > physical assumption about `Φ(t)`.
   >
   > **Remaining gap:** even granting the polynomial family, each grid point is still
   > built by literal Taylor-point-matching (matching `k_max` derivatives at `t_ref`),
   > not by fitting the true minimax/Chebyshev-projected polynomial of the same degree.
   > These generically differ, and the minimax one has strictly smaller worst-case
   > sup-norm error for the same degree (the "Chebyshev economization" result — Taylor
   > truncation is *not* the best degree-`n` approximation, Chebyshev truncation nearly
   > is). Appendix D fixes the conditioning/correlation problem (grid spacing), but does
   > not replace derivative-matching with a genuine minimax fit at the level of what each
   > candidate `Φ(t)` actually is — so there is, in principle, additional unexploited
   > headroom (a coarser grid, or lower `k_max`, for the same tolerance).

## Section 3.1/3.2 — Coherent/polynomial phase model

1. In the formula for `Phi(t)`, we need `Phi_ref`. Do we get it from the data? Do we
   enumerate it? Do we assume it is 0?

## Section 3.3.1 — Fourier-Domain Folding

1. The formula in (13) [eq:fourier_folding] does not quantize `Phi(t_n)` (not even to be
   in `N_b^{-1} Z`). With such a quantization, I can see how an IDFT will yield `P(b)`
   directly, but without it, cancellation won't be exact, and requires some quantifying.
2. The IDFT itself is not standard, as it uses only `floor(N_b/2)` harmonics instead of
   `N_b`, and then uses the real part and the value at 0. This seems to use some symmetry.
   It might be standard, but I would like to know the derivation. This might be related to
   the sentence - "the Nyquist-frequency term requires separate treatment when `N_b` is
   even."
3. It says the primary drawback of this method is the steep computational cost, measured
   as `O(N_s N_h)`. Shouldn't we be able to do it faster using standard FFT/DFT tricks?

## Section 3.3.1 (cont.) — Optimal Detection Statistics / template T(b)

1. Where is the template shape `T(b)` coming from? If the `T_n` were actually i.i.d., and
   the model were correct, there shouldn't be a need for such a template. Does it
   compensate for a caveat in the model? Maybe the right thing is to model the `T_n` with
   some dependency (e.g. a covariance matrix)?
2. More on the previous question - maybe the model used is actually
   `T_n = mu_n (1 + A T(b_n)) + eps_n`, assuming a template `T(b_n)`, and then it makes
   sense to try to get an MLE estimate for `A`, which seems to be what's happening here.
   Why is this a good model?
3. The paper states it defers more sophisticated template families and Fourier-domain
   detection statistics to future work. Are these already implemented in the current
   code?
4. The paper states it implements the `SNR_alpha` statistic, but there is a tradeoff for
   using the `SNR_beta` statistic. Has this been tried?
5. What is a boxcar template? Is it a step function (characteristic function)?

## Section 3.4 — Search Grid Design

1. Is `Delta Phi(t) = Phi(t) - Phi(t_ref)` or is it
   `Delta Phi(t_n) = Phi(t_n) - Phi(t_{n-1})`?
2. Why don't we determine the grid spacing for increasing `k` adaptively? I.e. after
   determining `delta f_k` grid spacing, for any given range of `f_k`, we can determine
   the maximal value of `Delta f_{k+1}` given the range of `f_0, ..., f_k`. This way we
   get an adaptive grid of points that zooms in or out according to our needs.
3. Why does the fact that the monomial basis functions are not orthogonal over the
   observation span lead to strong correlation between the model parameters?
4. It could be that the Chebyshev approach makes the adaptive idea (question 2 above)
   redundant. But then why not try to predict the Chebyshev coefficients from the
   beginning? Also, appendix D3 seems to be relevant regardless of the polynomial basis
   being used. Last point about this - Pravir had told me that the Chebyshev basis for the
   grid seems to not work as well - is that true? Why?

## Section 4.1 — Algorithm Description (P-FFA)

1. Algorithm 1 seems highly susceptible to parallelization - are any FFT-style speed-ups
   implemented?

## Section 4.3 — Sensitivity Loss and Error Bounds

1. Can't we solve the grid discretization issue by forcing our chosen grids such that
   `G_{i+1}` always refines `G_i`?

## Section 4.3.2 — (Fourier-domain folding trade-off)

1. Fourier-domain folding prevents phase-shift quantization errors, but costs more work.
   Can we find a more continuous trade-off of work <-> error between time-domain and
   Fourier-domain folding (guessing: by keeping a varying number of harmonic sums)?

## Section 5.1 — Algorithm Description (Extreme Pruning)

1. Why do we add the same size of segment in every stage and not increase their sizes
   dyadically?
