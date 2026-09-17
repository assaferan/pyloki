"""Power calculation for the proposed `tiling_strategy` injection campaign.

Arm under analysis: **Chebyshev basis, `aggressive` vs `quadrature`, `branch_max = 16`**
-- the only pairing with a proven geometric gain (D60) that builds at the shipped
default (D62), and therefore the only one that tests behaviour a user gets for free.

The quantity powered against is NOT a mean shift in S/N. The hypothesised mechanism is
bimodal: a signal whose covering leaf is thresholded away at some stage is lost
entirely, while one that stays clear of every threshold is unaffected. So the design is
a paired **McNemar** test and the effect size is the *discordance* structure

    p_disc = P(recovered_A != recovered_B),
    pi     = P(quadrature wins | discordant),

with the null pi = 1/2. This module estimates both from measured geometry plus one
explicitly stated score model, and converts them into a required n.

--------------------------------------------------------------------------------
THE SCORE MODEL, STATED IN FULL (this is where a hidden assumption would do the most
damage, so it is written out rather than buried in code)
--------------------------------------------------------------------------------

M1. **The score is a matched-filter S/N in units of the noise sigma.**
    `scoring.snr_score_batch_func` whitens per bin (`ts_e/sqrt(ts_v)`) and correlates
    with a zero-mean, unit-L2-norm boxcar, so under H0 the per-leaf score is N(0, 1)
    and the `threshold_scheme` values are directly comparable to it. Not an assumption
    about the physics -- read off the implementation.

M2. **Accumulation.** `common.shift_add_batch` adds `ts_e` and `ts_v` linearly across
    segments, so after s+1 of `nseg` segments the on-grid mean score is

        mu(s) = snr_final * sqrt((s + 1) / nseg).

    This is exactly the model the shipped ladder is built on (`schemes.bound_scheme`
    sets S/N^2 linear in the stage index), so using it here keeps the power calculation
    on the same footing as the thresholds it is powered against.

M3. **Phase error enters multiplicatively in amplitude, not additively in score.**
    A phase residual dPhi(t) smears the folded profile; `amplitude_loss.snr_ratio`
    computes the *exact* smeared profile via the characteristic function of dPhi and
    scores it with the pipeline's own filter bank (D53-D54). Writing L_X(s) for the
    resulting fractional S/N loss of arm X's best covering leaf at stage s,

        score_X(s) = mu(s) * (1 - L_X(s)) + N_X(s),   N_X ~ N(0, 1).

    THIS IS THE PHASE-ERROR-TO-SCORE MAPPING. It is multiplicative because the loss is
    a loss of *recovered amplitude*, and the noise level is unchanged by folding with a
    different ephemeris (the same samples are averaged either way). A sup-norm phase
    error must NOT be used in its place: D56 measured that doing so overestimates the
    loss by ~3x.

M4. **The two arms' noise is CORRELATED, and the correlation is the effect's
    denominator.** This is the point an unpaired power calculation gets wrong. Pairing
    fixes the data, but the two arms score *different templates*, so

        D(s) = score_B(s) - score_A(s) = mu(s) * (L_A - L_B)  +  eps,
        eps ~ N(0, 2 * (1 - rho_AB(s))),

    where rho_AB is the overlap of the two arms' templates. rho_AB is computed with the
    same smearing machinery applied to the phase residual *between the two templates*,
    delta_A - delta_B, which is legitimate because the residual is linear in the
    coefficient vector. Pairing removes the data realisation from the comparison; it
    does not remove the template-to-template noise difference, and that term turns out
    to dominate the deterministic one.

M5. **Single-crossing.** A pair is counted as discordant when, at some stage, the
    threshold separates the two arms' scores. Survival is a conjunction over stages;
    this sums the per-stage marginal probability weighted by the probability of having
    survived that far, which is an upper bound on p_disc (it double-counts pairs that
    are marginal at more than one stage) and is therefore ANTI-conservative for n.
    Flagged as such in the report.

M6. **One covering leaf per arm per stage.** Both arms have many leaves near the
    signal; survival is really `max` over them. Ignoring the max understates
    `quadrature`, whose redundancy gives it more draws. That is a separate mechanism
    from template proximity, it is *not* what this campaign is designed to measure, and
    the equal-P_d recalibration is what is supposed to null it. See the report's
    "confounds" section -- it is the reason the equal-P_d recalibration is not optional.

--------------------------------------------------------------------------------
Usage
--------------------------------------------------------------------------------

    python injection_power.py --geometry     # recompute the per-cell geometry (~20 min)
    python injection_power.py                # power table from the cached geometry

Writes/reads `injection_power.json` next to this file.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy import stats

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import phase3_config as P  # noqa: E402
from amplitude_loss import snr_ratio  # noqa: E402
from nearest_template_cheby import build_sets_cheby, cheby_basis  # noqa: E402
from nearest_template_cheby import min_excursion as min_excursion_cheby  # noqa: E402

from pyloki.utils.misc import C_VAL  # noqa: E402

JSON_PATH = HERE / "injection_power.json"
SCHEME_PATH = HERE / "schemes" / "cheby_{}.npz"

F0 = 1.0 / P.PERIOD
PO = P.POLY_ORDER
NBINS = P.NBINS
NSEG = 64
SEED = 20260917          # the same six signal positions the report uses
N_SIGNALS = 6
# The row 04_upstream_report.md quotes; 0.05 is marginal at N_b = 64.
DUCY = 0.10
DUCIES = (0.05, 0.10, 0.20)   # D55: duty cycle is the dominant uncertainty, ~10x across
SNR_FINAL = 10.0         # the scheme's calibration point
ARMS = ("aggressive", "quadrature")

# Denser at small s on purpose: the Viterbi ladder loses essentially all of its P_d
# before stage ~30 (cumulative H1 success is flat afterwards), so that is where the
# survival decision -- and therefore the whole effect -- actually lives.
STAGES = [2, 4, 6, 8, 10, 12, 14, 17, 20, 24, 28, 32, 40, 50, 62]


def truths() -> np.ndarray:
    cfg = P.make_config("aggressive")
    dp = np.array(cfg.get_dparams_actual(cfg.niters_ffa), dtype=float)
    dp[-1] *= C_VAL / F0
    rng = np.random.default_rng(SEED)
    return np.array([0.001, 0.05, 1.0, 0.0]) + rng.uniform(
        -0.5, 0.5, size=(N_SIGNALS, PO)) * dp


# --------------------------------------------------------------------------- geometry

def _best_leaf(strategy: str, truth: np.ndarray, stage: int, keep: int = 30):
    cfg = P.make_config(strategy)
    tol = cfg.eta / cfg.nbins
    sets, x, sc = build_sets_cheby(cfg, strategy, PO, F0, truth, stage, n_seed=1)
    sup, exact, _, top = min_excursion_cheby(
        sets, x, sc, PO, tol, max_nodes=300_000, keep=keep)
    return sup, bool(exact), x, top


def _pick_by_loss(top, x, ducy=DUCY):
    """The leaf among the candidates that actually loses least S/N.

    The search minimises the sup-norm, not the loss, so the sup-norm-nearest leaf need
    not be the best-scoring one (D56). Picking by loss over the 30 retained candidates
    tightens the upper bound the same way `amplitude_loss.loss_for_leaves` does.
    """
    ratios = [snr_ratio(d, x, F0, PO, NBINS, ducy, filt="matched",
                        basis_fn=cheby_basis) for _, d in top]
    i = int(np.argmax(ratios))
    return top[i][1], 1.0 - ratios[i]


def geometry() -> list[dict]:
    """Per (signal, stage): L_aggressive, L_quadrature, and the CROSS overlap rho_AB."""
    rows = []
    for ti, tr in enumerate(truths()):
        for s in STAGES:
            sup_a, ex_a, x, top_a = _best_leaf("aggressive", tr, s)
            sup_b, ex_b, _, top_b = _best_leaf("quadrature", tr, s)
            row = {"signal": ti, "stage": s,
                   "sup_aggressive": float(sup_a), "sup_quadrature": float(sup_b),
                   "exact_aggressive": ex_a, "exact_quadrature": ex_b}
            for du in DUCIES:
                d_a, l_a = _pick_by_loss(top_a, x, du)
                d_b, l_b = _pick_by_loss(top_b, x, du)
                # rho_AB: the overlap of the two arms' templates, from the phase
                # residual BETWEEN them. The residual is linear in the coefficient
                # vector, so delta_A - delta_B is the right argument to the same
                # smearing model that produced the losses.
                rho = snr_ratio(d_a - d_b, x, F0, PO, NBINS, du,
                                filt="matched", basis_fn=cheby_basis)
                row[f"ducy_{du:.2f}"] = {
                    "L_aggressive": float(l_a), "L_quadrature": float(l_b),
                    "rho_AB": float(rho), "one_minus_rho_AB": float(1.0 - rho)}
            rows.append(row)
            g = row[f"ducy_{DUCY:.2f}"]
            print(f"  sig{ti} S={s:2d} supA={sup_a:7.3f} supB={sup_b:7.3f} "
                  f"L_A={g['L_aggressive']:.5f} L_B={g['L_quadrature']:.5f} "
                  f"1-rho={g['one_minus_rho_AB']:.5f} "
                  f"{'' if (ex_a and ex_b) else '[NOT EXACT]'}", flush=True)
    return rows


# ------------------------------------------------------------------------ power model

def _folded_normal_moments(mean: np.ndarray, sd: np.ndarray):
    """E[D+] and E[|D|] for D ~ N(mean, sd^2), elementwise.

    E[D+]  = sd*phi(d) + mean*Phi(d)
    E[|D|] = 2*sd*phi(d) + mean*(2*Phi(d) - 1)        with d = mean/sd
    """
    sd = np.maximum(sd, 1e-9)
    d = mean / sd
    e_pos = sd * stats.norm.pdf(d) + mean * stats.norm.cdf(d)
    e_abs = 2.0 * sd * stats.norm.pdf(d) + mean * (2.0 * stats.norm.cdf(d) - 1.0)
    return e_pos, e_abs


def stage_weights(
    strategy: str = "aggressive",
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(thresholds, mu, survival-to-stage) for the recalibrated Chebyshev ladder."""
    z = np.load(str(SCHEME_PATH).format(strategy))
    thr = np.asarray(z["thresholds"], dtype=float)
    surv = np.asarray(z["succ_h1"], dtype=float)          # cumulative, per stage
    n = len(thr)
    mu = SNR_FINAL * np.sqrt((np.arange(n) + 1.0) / NSEG)
    return thr, mu, surv


def power(rows: list[dict], strategy: str = "aggressive", ducy: float = DUCY) -> dict:
    key = f"ducy_{ducy:.2f}"
    thr, mu_all, surv_all = stage_weights(strategy)
    ok = [r for r in rows if r["exact_aggressive"] and r["exact_quadrature"]]
    by_stage: dict[int, list[dict]] = {}
    for r in ok:
        by_stage.setdefault(r["stage"], []).append(r)

    num_pos = num_abs = p_disc = 0.0
    per_stage = []
    for s in sorted(by_stage):
        cells = by_stage[s]
        mu = mu_all[s]
        t = thr[s]
        l_a = np.array([c[key]["L_aggressive"] for c in cells])
        l_b = np.array([c[key]["L_quadrature"] for c in cells])
        omr = np.array([c[key]["one_minus_rho_AB"] for c in cells])
        delta = mu * (l_a - l_b)                       # deterministic part of D
        sd = np.sqrt(2.0 * np.maximum(omr, 0.0))       # template-noise part of D
        e_pos, e_abs = _folded_normal_moments(delta, sd)
        # Density of arm A's score at the stage threshold. The score is
        # N(mu*(1-L_A), 1), so the density at t is phi(t - mu*(1-L_A)).
        dens = stats.norm.pdf(t - mu * (1.0 - l_a))
        # Probability of still being in the race at stage s.
        reach = surv_all[s - 1] if s >= 1 else 1.0
        w = reach * dens
        s_pos, s_abs = float(np.mean(w * e_pos)), float(np.mean(w * e_abs))
        num_pos += s_pos
        num_abs += s_abs
        p_disc += s_abs
        per_stage.append({
            "stage": s, "cells": len(cells), "threshold": round(float(t), 3),
            "mu": round(float(mu), 3), "reach": round(float(reach), 4),
            "density_at_threshold": round(float(np.mean(dens)), 5),
            "median_delta": round(float(np.median(delta)), 5),
            "median_sd": round(float(np.median(sd)), 5),
            "median_delta_over_sd": round(
                float(np.median(delta / np.maximum(sd, 1e-9))), 4),
            "mean_E_absD": round(float(np.mean(e_abs)), 5),
            "stage_discordance": round(s_abs, 6),
            "stage_pi": round(s_pos / s_abs, 4) if s_abs > 0 else None,
        })

    pi = num_pos / num_abs if num_abs > 0 else 0.5
    # STAGES is a subsample of the 63 stages; each sampled stage stands for the
    # 63/len(STAGES) stages around it. Scaling p_disc by that factor is the honest
    # (and, under M5, generous) way to aggregate.
    scale = 63.0 / len(STAGES)
    return {
        "ducy": ducy, "cells_used": len(ok), "cells_total": len(rows),
        "pi": round(float(pi), 5),
        "p_disc_sampled": round(float(p_disc), 6),
        "p_disc": round(float(min(p_disc * scale, 1.0)), 5),
        "per_stage": per_stage,
    }


def mcnemar_n(pi: float, p_disc: float, alpha: float = 0.05,
              beta: float = 0.20) -> dict:
    """Pairs needed for a two-sided McNemar test at (alpha, 1-beta).

    Conditional on the m discordant pairs, McNemar is a binomial test of pi = 1/2, so

        m = (z_{alpha/2}*sqrt(1/4) + z_beta*sqrt(pi*(1-pi)))^2 / (pi - 1/2)^2
        n = m / p_disc.
    """
    za = stats.norm.isf(alpha / 2.0)
    zb = stats.norm.isf(beta)
    eff = abs(pi - 0.5)
    if eff < 1e-9:
        return {"pi": pi, "p_disc": p_disc, "m_discordant": float("inf"),
                "n_pairs": float("inf")}
    m = (za * 0.5 + zb * np.sqrt(pi * (1.0 - pi))) ** 2 / eff**2
    return {"pi": round(pi, 5), "p_disc": round(p_disc, 5),
            "m_discordant": float(np.ceil(m)),
            "n_pairs": float(np.ceil(m / p_disc))}


# ---------------------------------------------------------------------------- strata

# Stages 1..~30 carry essentially all of the ladder's P_d loss (the cumulative H1
# success is flat after ~30), so a stratum has to be defined on these and not on the
# late stages where the tiling geometry is largest but nothing is being decided.
DECISION_STAGES = [2, 6, 10, 14, 20, 28]


def strata(n_positions: int = 24, ducy: float = DUCY) -> dict:
    """Is the 'null' stratum constructible, or only assumed?

    The brief asks for stratification by within-cell position, *measured*. A signal is
    placed once in 4-D parameter space; its per-stage geometry then follows and is not
    separately dial-able. So the honest construction is: sample positions, measure the
    decision-weighted effect size for each, and bin. This reports whether the bins are
    populated -- in particular whether positions exist with Delta ~ 0 at every decision
    stage, which is what the built-in null requires.
    """
    thr, mu_all, surv_all = stage_weights("aggressive")
    cfg = P.make_config("aggressive")
    dp = np.array(cfg.get_dparams_actual(cfg.niters_ffa), dtype=float)
    dp[-1] *= C_VAL / F0
    rng = np.random.default_rng(SEED + 1)
    base = np.array([0.001, 0.05, 1.0, 0.0])
    positions = base + rng.uniform(-0.5, 0.5, size=(n_positions, PO)) * dp

    out = []
    for pi_, tr in enumerate(positions):
        num_pos = num_abs = 0.0
        detail = []
        for st in DECISION_STAGES:
            sup_a, ex_a, x, top_a = _best_leaf("aggressive", tr, st)
            sup_b, ex_b, _, top_b = _best_leaf("quadrature", tr, st)
            if not (ex_a and ex_b):
                continue
            d_a, l_a = _pick_by_loss(top_a, x, ducy)
            d_b, l_b = _pick_by_loss(top_b, x, ducy)
            rho = snr_ratio(d_a - d_b, x, F0, PO, NBINS, ducy,
                            filt="matched", basis_fn=cheby_basis)
            mu, t = mu_all[st], thr[st]
            delta = mu * (l_a - l_b)
            sd = float(np.sqrt(max(2.0 * (1.0 - rho), 0.0)))
            e_pos, e_abs = _folded_normal_moments(
                np.array([delta]), np.array([sd]))
            w = surv_all[st - 1] * stats.norm.pdf(t - mu * (1.0 - l_a))
            num_pos += float(w * e_pos[0])
            num_abs += float(w * e_abs[0])
            detail.append({"stage": st, "sup_aggressive": round(float(sup_a), 4),
                           "sup_quadrature": round(float(sup_b), 4),
                           "delta": round(float(delta), 5), "sd": round(sd, 5)})
        pi_val = num_pos / num_abs if num_abs > 0 else 0.5
        out.append({"position": pi_, "pi_decision_weighted": round(float(pi_val), 4),
                    "weighted_E_absD": round(float(num_abs), 6),
                    "max_abs_delta": round(
                        float(max((abs(d["delta"]) for d in detail), default=0.0)), 5),
                    "stages": detail})
        print(f"  pos{pi_:2d} pi={pi_val:6.3f}  w*E|D|={num_abs:.5f}  "
              f"max|Delta|={out[-1]['max_abs_delta']:.4f}",
              flush=True)
    return {"ducy": ducy, "decision_stages": DECISION_STAGES, "positions": out}


# -------------------------------------------------------------------------------- CLI

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--geometry", action="store_true",
                    help="recompute the per-cell geometry (~5 min)")
    ap.add_argument("--strata", type=int, default=0, metavar="N",
                    help="measure the decision-weighted effect for N sampled positions")
    ap.add_argument("--sec-per-pair", type=float, default=6.0,
                    help="measured wall clock per paired run, warm (see the report)")
    args = ap.parse_args()

    if args.geometry:
        rows = geometry()
        JSON_PATH.write_text(json.dumps({"config": {
            "basis": "chebyshev", "arms": list(ARMS), "branch_max": 16,
            "poly_order": PO, "nbins": NBINS, "nsegments": NSEG, "eta": P.ETA,
            "ducies": list(DUCIES), "snr_final": SNR_FINAL, "signals": N_SIGNALS,
            "seed": SEED, "stages": STAGES},
            "geometry": rows}, indent=1) + "\n")
    data = json.loads(JSON_PATH.read_text())
    if args.strata:
        data["strata"] = strata(args.strata)
        q = np.array([p["pi_decision_weighted"] for p in data["strata"]["positions"]])
        m = np.array([p["max_abs_delta"] for p in data["strata"]["positions"]])
        print(f"\n--- strata over {len(q)} positions (ducy = {DUCY:.2f}) ---")
        print(f"  decision-weighted pi: min {q.min():.3f} "
              f"q25 {np.percentile(q, 25):.3f} median {np.median(q):.3f} "
              f"q75 {np.percentile(q, 75):.3f} max {q.max():.3f}")
        print(f"  positions with pi < 0.55 (NULL):   {(q < 0.55).sum()}/{len(q)}")
        print(f"  positions with pi > 0.75 (EFFECT): {(q > 0.75).sum()}/{len(q)}")
        print(f"  max|Delta| over decision stages: median {np.median(m):.4f} "
              f"min {m.min():.4f} max {m.max():.4f}")
        JSON_PATH.write_text(json.dumps(data, indent=1) + "\n")
    rows = data["geometry"]
    ok = [r for r in rows if r["exact_aggressive"] and r["exact_quadrature"]]

    summaries, powers = {}, {}
    for du in DUCIES:
        k = f"ducy_{du:.2f}"
        l_a = np.array([r[k]["L_aggressive"] for r in ok])
        l_b = np.array([r[k]["L_quadrature"] for r in ok])
        omr = np.array([r[k]["one_minus_rho_AB"] for r in ok])
        summaries[k] = {
            "median_L_aggressive": round(float(np.median(l_a)), 6),
            "median_L_quadrature": round(float(np.median(l_b)), 6),
            "paired_advantage": round(float(np.median((1 - l_b) / (1 - l_a)) - 1.0), 6),
            "median_one_minus_rho_AB": round(float(np.median(omr)), 6),
            "median_sigma_template_noise": round(float(np.sqrt(2 * np.median(omr))), 6),
        }
        powers[k] = power(rows, ducy=du)

    print(f"\ncells used = {len(ok)}/{len(rows)} (a cell is dropped when either "
          "branch-and-bound hit its node budget, because a budget-limited value "
          "is only an upper bound)")
    print("\n--- geometry, by duty cycle ---")
    hdr = ("ducy", "L_agg", "L_quad", "paired adv", "1-rho_AB", "sigma_noise")
    print("{:>6} {:>9} {:>9} {:>11} {:>10} {:>12}".format(*hdr))
    for du in DUCIES:
        g = summaries[f"ducy_{du:.2f}"]
        print(f"{du:6.2f} {g['median_L_aggressive']:9.5f} "
              f"{g['median_L_quadrature']:9.5f} {g['paired_advantage']:11.5f} "
              f"{g['median_one_minus_rho_AB']:10.5f} "
              f"{g['median_sigma_template_noise']:12.4f}")

    print(f"\n--- per stage (ducy = {DUCY:.2f}) ---")
    print(f"{'S':>3} {'T':>5} {'mu':>6} {'reach':>6} {'dens':>7} {'Delta':>8} "
          f"{'sd':>7} {'D/sd':>7} {'E|D|':>7} {'p_s':>8} {'pi_s':>6}")
    for r in powers[f"ducy_{DUCY:.2f}"]["per_stage"]:
        print(f"{r['stage']:3d} {r['threshold']:5.2f} {r['mu']:6.2f} {r['reach']:6.3f} "
              f"{r['density_at_threshold']:7.4f} {r['median_delta']:8.4f} "
              f"{r['median_sd']:7.4f} {r['median_delta_over_sd']:7.3f} "
              f"{r['mean_E_absD']:7.4f} {r['stage_discordance']:8.5f} "
              f"{(r['stage_pi'] or float('nan')):6.3f}")

    print("\n--- effect size and required n, 80% power, two-sided alpha = 0.05 ---")
    print(f"{'ducy':>6} {'pi':>7} {'p_disc':>8} {'m_disc':>9} {'n_pairs':>10} "
          f"{'core-h':>9}")
    for du in DUCIES:
        pw = powers[f"ducy_{du:.2f}"]
        r = mcnemar_n(pw["pi"], pw["p_disc"])
        print(f"{du:6.2f} {pw['pi']:7.4f} {pw['p_disc']:8.4f} "
              f"{r['m_discordant']:9.0f} {r['n_pairs']:10.0f} "
              f"{r['n_pairs'] * args.sec_per_pair / 3600.0:9.1f}")

    print("\n--- n across a range of pi and p_disc "
          "(the point estimates are one cell) ---")
    print(f"{'pi':>6} | " + " ".join(f"{p:>9.3f}" for p in (0.01, 0.03, 0.10, 0.30)))
    print("       | " + " ".join(f"{'n pairs':>9}" for _ in range(4)))
    for pi in (0.52, 0.55, 0.60, 0.65, 0.70, 0.80, 0.90):
        cells = [f"{mcnemar_n(pi, pd)['n_pairs']:9.0f}"
                 for pd in (0.01, 0.03, 0.10, 0.30)]
        print(f"{pi:6.2f} | " + " ".join(cells))
    print("\nwall clock = n_pairs * sec_per_pair; "
          f"at the measured {args.sec_per_pair:.1f} s/pair, 10 000 pairs = "
          f"{10000 * args.sec_per_pair / 3600.0:.1f} core-hours")

    data["summary"] = summaries
    data["power"] = powers
    data["sec_per_pair"] = args.sec_per_pair
    JSON_PATH.write_text(json.dumps(data, indent=1) + "\n")


if __name__ == "__main__":
    main()
