"""Is `rho_AB` the SCORE correlation, or only the PROFILE overlap? — §11.4 (1).

M4 of the power model says the paired difference of the two arms' scores carries a
stochastic term of sd `sqrt(2(1 - rho_AB))`, and that term is measured to be three times
the deterministic one. `n ∝ (pi - 1/2)^-2`, so `rho_AB` is the single number the whole
campaign size turns on (§10, row 3). It has never been checked.

What `injection_power.py` actually computes is

    rho_prof = snr_ratio(delta_A - delta_B, ..., filt="matched")
             = <p, p*s> / <p, p>

with `p` the true pulse and `s` the smearing kernel of the phase residual between the
two arms' templates. That is the overlap of two *signals*. What M4 needs is the
correlation of two *scores* under noise. Writing the score as a linear filter `h` on the
folded profile, and using that folding on two ephemerides assigns the same samples to
bins differing by dPhi(t),

    rho_score = sum_k |H(k)|^2 S(k) / sum_k |H(k)|^2
    rho_prof  = sum_k |P(k)|^2 S(k) / sum_k |P(k)|^2

— the SAME smearing factor S(k), weighted by the filter's power spectrum instead of the
pulse's. The two agree only if the filter IS the pulse. The pipeline's filter is a
boxcar, and at ducy 0.10 with 64 bins the bank selects width 3, whose spectrum is far
broader than a Gaussian pulse's. A broader filter weights high k harder, S(k) is smaller
there, so the substitution should OVERSTATE rho — understate the noise on the paired
difference, overstate pi, and understate n. The direction was predicted before
measuring; the magnitude is what this script measures.

Four quantities, on the SAME cells and through the SAME code path as
`injection_power.py`:

1. `rho_prof`    — what is used today.
2. `rho_fixed`   — closed form, for the single boxcar the bank selects for the pulse.
3. `rho_maxbank` — simulated under H0 with the pipeline's own scorer, max over widths
                   and phase. With no signal the two arms' maxima are free to land on
                   different filters, so this is a LOWER bound on the correlation.
4. `rho_signal`  — simulated with a signal present at the stage's own mu, which is the
                   campaign's actual situation: the signal anchors both maxima to the
                   same filter. **This is the one the power model should use.**

The simulation draws the two folded noise profiles from their exact joint law. Under H0
with white data, bin k of arm A and bin j of arm B have covariance proportional to the
number of samples they share, whose Fourier transform is exactly S(k). So the pair is
generated in Fourier as

    P_B(k) = S(k) P_A(k) + sqrt(1 - |S(k)|^2) Z(k),

which is exact rather than a Monte-Carlo of the time series. It assumes only what the
smearing model already assumes — that the fast phase is uniform and independent of the
slowly varying residual — so it isolates the filter question and nothing else. The
sampler is checked against the closed form on a one-width bank before anything else
runs (`selftest`).

Finally the corrected per-cell `1 - rho` is substituted into the cached geometry and
pushed back through `injection_power.power` and `mcnemar_n`, so the output is not a
correlation coefficient but the campaign size it implies.

    python rho_check.py                      # all 90 cells, ~10 min
    python rho_check.py --cells 6 --draws 5000    # quick
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import injection_power as IP  # noqa: E402
from amplitude_loss import smearing_factor, snr_ratio  # noqa: E402
from nearest_template_cheby import cheby_basis  # noqa: E402

from pyloki.detection.scoring import (  # noqa: E402
    boxcar_snr_1d,
    boxcar_snr_2d,
    generate_box_width_trials,
)
from pyloki.simulation.pulse import generate_folded_profile  # noqa: E402

DUCY_MAX = 0.2


# --------------------------------------------------------------------- spectra


def _full_spectrum(s_half: np.ndarray, nbins: int) -> np.ndarray:
    """Hermitian extension of S(k), k = 0..nbins//2, to the full fft grid."""
    full = np.empty(nbins, dtype=complex)
    full[: nbins // 2 + 1] = s_half
    full[nbins // 2 + 1:] = np.conj(s_half[1:(nbins + 1) // 2][::-1])
    return full


def _boxcar(width: int, nbins: int) -> np.ndarray:
    """The zero-mean unit-L2 boxcar the scorer uses (scoring.py:136-139)."""
    size_w = nbins - width
    height = np.sqrt(size_w / (nbins * width))
    h = np.full(nbins, -width * height / size_w)
    h[:width] = height
    return h


def _best_width(ducy: float, nbins: int, widths: np.ndarray) -> int:
    prof = generate_folded_profile(nbins=nbins, ducy=ducy, center=0.5)
    return int(widths[int(np.argmax(boxcar_snr_1d(prof.astype(np.float32),
                                                  widths, 1.0)))])


def _weighted_rho(power: np.ndarray, s_full: np.ndarray) -> float:
    """sum_k |W(k)|^2 S(k) / sum_k |W(k)|^2, real by Hermitian symmetry."""
    return float(np.real((power * s_full).sum() / power.sum()))


# ------------------------------------------------------------------ simulation


def _draw_pair(s_full: np.ndarray, nbins: int, draws: int,
               rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    """Two folded unit-variance noise profiles with cross-spectrum S(k)."""
    amp = np.sqrt(np.maximum(1.0 - np.abs(s_full) ** 2, 0.0))
    a = rng.standard_normal((draws, nbins))
    z = rng.standard_normal((draws, nbins))
    fb = s_full[None, :] * np.fft.fft(a, axis=1) + amp[None, :] * np.fft.fft(z, axis=1)
    return a, np.real(np.fft.ifft(fb, axis=1))


def _maxbank(x: np.ndarray, widths: np.ndarray) -> np.ndarray:
    return boxcar_snr_2d(np.ascontiguousarray(x, dtype=np.float32),
                         widths, 1.0).max(1)


def simulate(s_full: np.ndarray, widths: np.ndarray, nbins: int, draws: int,
             signal: np.ndarray | None, rng: np.random.Generator) -> float:
    """Correlation of the pipeline's score (max over widths and phase).

    `signal` is added identically to both arms. The two arms' signals differ only by
    L_A - L_B, which is ~1% and cannot matter for the correlation of the noise term;
    what the signal is here for is to anchor both maxima to the same filter, which is
    what makes this the campaign's situation rather than a pure-noise one.
    """
    a, b = _draw_pair(s_full, nbins, draws, rng)
    if signal is not None:
        a = a + signal[None, :]
        b = b + signal[None, :]
    return float(np.corrcoef(_maxbank(a, widths), _maxbank(b, widths))[0, 1])


def selftest(nbins: int, draws: int, rng: np.random.Generator) -> None:
    """The sampler against the closed form, on a single fixed filter.

    With one boxcar and no phase max the score is a fixed linear functional, so
    `_weighted_rho` is exact and the sampler must reproduce it. If this fails, nothing
    below it means anything.
    """
    w = 8
    s_half = np.exp(-0.5 * (0.35 * np.arange(nbins // 2 + 1)) ** 2).astype(complex)
    s_full = _full_spectrum(s_half, nbins)
    h = _boxcar(w, nbins)
    a, b = _draw_pair(s_full, nbins, draws, rng)
    sim = float(np.corrcoef(a @ h, b @ h)[0, 1])
    exact = _weighted_rho(np.abs(np.fft.fft(h)) ** 2, s_full)
    se = (1 - sim**2) / np.sqrt(draws)
    print(f"selftest  simulated rho={sim:.5f}  closed form={exact:.5f}  "
          f"diff={abs(sim - exact):.5f}  (MC se {se:.5f})")
    if abs(sim - exact) > max(5 * se, 0.005):
        raise SystemExit("sampler disagrees with the closed form")


# -------------------------------------------------------------------- measure


def measure(cells: int, draws: int, ducy: float) -> tuple[list[dict], dict]:
    nbins = IP.NBINS
    widths = generate_box_width_trials(nbins, ducy_max=DUCY_MAX)
    w_best = _best_width(ducy, nbins, widths)
    pw_box = np.abs(np.fft.fft(_boxcar(w_best, nbins))) ** 2
    prof0 = generate_folded_profile(nbins=nbins, ducy=ducy, center=0.5)
    pw_prof = np.abs(np.fft.fft(prof0 - prof0.mean())) ** 2
    rng = np.random.default_rng(20260917)
    selftest(nbins, draws, rng)

    k = np.minimum(np.arange(nbins), nbins - np.arange(nbins))
    k2_box = float((pw_box * k**2).sum() / pw_box.sum())
    k2_prof = float((pw_prof * k**2).sum() / pw_prof.sum())
    print(f"nbins={nbins} ducy={ducy} bank={list(map(int, widths))} "
          f"selected width={w_best}")
    print(f"<k^2>  boxcar={k2_box:.2f}  pulse={k2_prof:.2f}  "
          f"ratio={k2_box / k2_prof:.2f}  "
          f"<- predicted (1-rho) inflation, small-smearing limit")

    # the signal used to anchor the max, scaled to the stage's own mu
    s0 = float(boxcar_snr_1d(prof0.astype(np.float32), widths, 1.0).max())
    _, mu_all, _ = IP.stage_weights("aggressive")

    rows = []
    for ti, tr in enumerate(IP.truths()):
        for s in IP.STAGES:
            if len(rows) >= cells:
                break
            _, ex_a, x, top_a = IP._best_leaf("aggressive", tr, s)
            _, ex_b, _, top_b = IP._best_leaf("quadrature", tr, s)
            if not (ex_a and ex_b):
                print(f"  sig{ti} S={s:2d}  skipped (node budget, not exact)",
                      flush=True)
                continue
            d_a, _ = IP._pick_by_loss(top_a, x, ducy)
            d_b, _ = IP._pick_by_loss(top_b, x, ducy)
            delta = d_a - d_b
            rho_prof = float(snr_ratio(delta, x, IP.F0, IP.PO, nbins, ducy,
                                       filt="matched", basis_fn=cheby_basis))
            s_full = _full_spectrum(
                smearing_factor(delta, x, IP.F0, IP.PO, nbins, basis_fn=cheby_basis),
                nbins)
            rho_fixed = _weighted_rho(pw_box, s_full)
            rho_h0 = simulate(s_full, widths, nbins, draws, None, rng)
            sig = prof0 * (float(mu_all[s]) / s0)
            rho_sig = simulate(s_full, widths, nbins, draws, sig, rng)
            rows.append({
                "signal": ti, "stage": s, "mu": float(mu_all[s]),
                "rho_prof": rho_prof, "rho_fixed_boxcar": rho_fixed,
                "rho_maxbank_h0": rho_h0, "rho_signal": rho_sig,
                "one_minus": {"prof": 1 - rho_prof, "fixed_boxcar": 1 - rho_fixed,
                              "maxbank_h0": 1 - rho_h0, "signal": 1 - rho_sig},
            })
            print(f"  sig{ti} S={s:2d}  1-rho  prof={1 - rho_prof:.5f} "
                  f"boxcar={1 - rho_fixed:.5f} signal={1 - rho_sig:.5f} "
                  f"H0={1 - rho_h0:.5f}", flush=True)
    return rows, {"boxcar_width": w_best, "k2_boxcar": k2_box,
                  "k2_pulse": k2_prof, "k2_ratio": k2_box / k2_prof,
                  "ducy": ducy, "draws": draws}


# ------------------------------------------------------------------ propagate


def propagate(rows: list[dict], ducy: float) -> dict:
    """Substitute the corrected 1-rho into the cached geometry and re-run the model.

    Only cells measured here are replaced; cells not measured keep their cached value,
    and the count of each is reported, because a partially corrected table is not a
    corrected table and must not be quoted as one.
    """
    cached = json.loads((HERE / "injection_power.json").read_text())
    geo = cached["geometry"]
    key = f"ducy_{ducy:.2f}"
    out = {}
    for variant in ("prof", "fixed_boxcar", "signal", "maxbank_h0"):
        by_cell = {(r["signal"], r["stage"]): r["one_minus"][variant] for r in rows}
        patched, n_sub = [], 0
        for g in geo:
            g2 = {k: (dict(v) if isinstance(v, dict) else v) for k, v in g.items()}
            omr = by_cell.get((g["signal"], g["stage"]))
            if omr is not None:
                g2[key]["one_minus_rho_AB"] = float(max(omr, 0.0))
                g2[key]["rho_AB"] = 1.0 - float(max(omr, 0.0))
                n_sub += 1
            patched.append(g2)
        pw = IP.power(patched, "aggressive", ducy)
        out[variant] = {
            "cells_substituted": n_sub, "cells_total": len(geo),
            "pi": pw["pi"], "p_disc": pw["p_disc"],
            "n_at_model_p_disc": IP.mcnemar_n(pw["pi"], pw["p_disc"])["n_pairs"],
            "n_at_measured_p_disc_0.30": IP.mcnemar_n(pw["pi"], 0.30)["n_pairs"],
            "n_at_measured_p_disc_0.08": IP.mcnemar_n(pw["pi"], 0.08)["n_pairs"],
        }
    return out


def strata(n_positions: int, draws: int, ducy: float) -> dict:
    """§5.1's stratification, recomputed with the corrected score correlation.

    `injection_power.strata` builds its per-position pi from the same
    `snr_ratio(delta_A - delta_B)` this script just showed is the wrong quantity
    (`injection_power.py:337`), so §5.1's distribution -- and the stratum boundaries
    §9 pre-commits to -- inherit the error. This recomputes it both ways on the SAME
    sampled positions, so the shift is visible and the old numbers are reproduced as a
    control rather than asserted.
    """
    nbins = IP.NBINS
    widths = generate_box_width_trials(nbins, ducy_max=DUCY_MAX)
    prof0 = generate_folded_profile(nbins=nbins, ducy=ducy, center=0.5)
    s0 = float(boxcar_snr_1d(prof0.astype(np.float32), widths, 1.0).max())
    thr, mu_all, surv_all = IP.stage_weights("aggressive")
    cfg = IP.P.make_config("aggressive")
    dp = np.array(cfg.get_dparams_actual(cfg.niters_ffa), dtype=float)
    dp[-1] *= IP.C_VAL / IP.F0
    rng = np.random.default_rng(IP.SEED + 1)
    positions = np.array([0.001, 0.05, 1.0, 0.0]) + rng.uniform(
        -0.5, 0.5, size=(n_positions, IP.PO)) * dp

    out = []
    for idx, tr in enumerate(positions):
        acc = {"prof": [0.0, 0.0], "signal": [0.0, 0.0]}
        for st in IP.DECISION_STAGES:
            _, ex_a, x, top_a = IP._best_leaf("aggressive", tr, st)
            _, ex_b, _, top_b = IP._best_leaf("quadrature", tr, st)
            if not (ex_a and ex_b):
                continue
            d_a, l_a = IP._pick_by_loss(top_a, x, ducy)
            d_b, l_b = IP._pick_by_loss(top_b, x, ducy)
            delta_c = d_a - d_b
            rho_p = float(snr_ratio(delta_c, x, IP.F0, IP.PO, nbins, ducy,
                                    filt="matched", basis_fn=cheby_basis))
            s_full = _full_spectrum(
                smearing_factor(delta_c, x, IP.F0, IP.PO, nbins,
                                basis_fn=cheby_basis), nbins)
            rho_s = simulate(s_full, widths, nbins, draws,
                             prof0 * (float(mu_all[st]) / s0), rng)
            mu, t = mu_all[st], thr[st]
            delta = mu * (l_a - l_b)
            w = surv_all[st - 1] * IP.stats.norm.pdf(t - mu * (1.0 - l_a))
            for tag, rho in (("prof", rho_p), ("signal", rho_s)):
                sd = float(np.sqrt(max(2.0 * (1.0 - rho), 0.0)))
                e_pos, e_abs = IP._folded_normal_moments(
                    np.array([delta]), np.array([sd]))
                acc[tag][0] += float(w * e_pos[0])
                acc[tag][1] += float(w * e_abs[0])
        row = {"position": idx}
        for tag in ("prof", "signal"):
            num_pos, num_abs = acc[tag]
            row[f"pi_{tag}"] = float(num_pos / num_abs) if num_abs > 0 else 0.5
        out.append(row)
        print(f"  pos{idx:2d}  pi  profile={row['pi_prof']:.4f}  "
              f"corrected={row['pi_signal']:.4f}", flush=True)

    res = {"n_positions": n_positions, "ducy": ducy, "positions": out}
    for tag in ("prof", "signal"):
        v = np.array([r[f"pi_{tag}"] for r in out])
        res[tag] = {
            "min": float(v.min()), "q25": float(np.percentile(v, 25)),
            "median": float(np.median(v)), "q75": float(np.percentile(v, 75)),
            "max": float(v.max()),
            "frac_null_below_0.55": float((v < 0.55).mean()),
            "frac_effect_above_0.70": float((v > 0.70).mean()),
        }
    print("\n" + f"{'':10}{'min':>8}{'q25':>8}{'median':>8}{'q75':>8}{'max':>8}"
          f"{'<0.55':>8}{'>0.70':>8}")
    for tag, label in (("prof", "profile"), ("signal", "corrected")):
        d = res[tag]
        print(f"{label:10}{d['min']:>8.3f}{d['q25']:>8.3f}{d['median']:>8.3f}"
              f"{d['q75']:>8.3f}{d['max']:>8.3f}"
              f"{d['frac_null_below_0.55']:>8.2f}{d['frac_effect_above_0.70']:>8.2f}")
    return res


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--strata", type=int, default=0, metavar="N",
                    help="recompute §5.1's stratification with the corrected rho")
    ap.add_argument("--cells", type=int, default=10**6)
    ap.add_argument("--draws", type=int, default=20000)
    ap.add_argument("--ducy", type=float, default=IP.DUCY)
    ap.add_argument("--out", type=Path, default=HERE / "rho_check.json")
    args = ap.parse_args()

    if args.strata:
        res = strata(args.strata, args.draws, args.ducy)
        args.out.write_text(json.dumps(res, indent=1))
        print(f"\nwritten to {args.out}")
        return

    rows, meta = measure(args.cells, args.draws, args.ducy)
    if not rows:
        raise SystemExit("no exact cells")
    prop = propagate(rows, args.ducy)

    def med(v):
        return float(np.median([r["one_minus"][v] for r in rows]))

    ok = [r for r in rows if r["one_minus"]["prof"] > 1e-4]
    infl = {v: (float(np.median([r["one_minus"][v] / r["one_minus"]["prof"]
                                 for r in ok])) if ok else float("nan"))
            for v in ("fixed_boxcar", "signal", "maxbank_h0")}
    summary = {**meta, "cells": len(rows),
               "cells_with_measurable_smearing": len(ok),
               "median_one_minus_rho": {v: med(v) for v in
                                        ("prof", "fixed_boxcar", "signal",
                                         "maxbank_h0")},
               "median_inflation_vs_prof": infl,
               "power": prop}
    args.out.write_text(json.dumps({"summary": summary, "cells": rows}, indent=1))

    print(f"\ncells {len(rows)}  (with measurable smearing: {len(ok)})")
    print(f"median 1-rho   profile overlap   {med('prof'):.5f}  "
          "<- what injection_power.py uses")
    for v, label in (("fixed_boxcar", "fixed boxcar    "),
                     ("signal", "with signal     "),
                     ("maxbank_h0", "max over bank H0")):
        print(f"median 1-rho   {label}  {med(v):.5f}  "
              f"({infl[v]:.2f}x)")
    print("\nwhat it does to the campaign (ducy "
          f"{args.ducy}, substituted into the cached geometry):")
    print(f"{'variant':<14}{'pi':>9}{'p_disc':>9}{'n @ p_disc=0.30':>18}"
          f"{'n @ 0.08':>12}")
    for v, d in prop.items():
        print(f"{v:<14}{d['pi']:>9.4f}{d['p_disc']:>9.4f}"
              f"{d['n_at_measured_p_disc_0.30']:>18.0f}"
              f"{d['n_at_measured_p_disc_0.08']:>12.0f}")
    print(f"\nwritten to {args.out}")


if __name__ == "__main__":
    main()
