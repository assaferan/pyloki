"""The paired `max_sugg` experiment: is the candidate-buffer bias arm-dependent?

This is §11.6 of `05_injection_design.md`, the one thing that has to be settled before
the campaign in §9 means anything. `injection_pilot.py` already established that the
buffer changes the search outcome on data the search has already seen (§6.4). What it
could not establish is whether it does so *symmetrically* between the arms. If it does
not -- if `quadrature`, whose branching product is 2^28.4 times `aggressive`'s, is
penalised harder by the overflow ratchet at `world_tree.py:527-551` -- then a campaign
run at any buffer either arm saturates is measuring the buffer, not the tiling.

Design: **one** fixed set of n realisations at a single S/N, run at all four
(arm, max_sugg) combinations. 2 arms x 2 buffers x n runs. Each combination is a
separate process, because building several configs in one interpreter is unstable (the
note in `run_injections.py`); this module only drives `injection_pilot.py --arm` and
renames its output, so the measurement code is the same code the pilot used.

THE ANALYSIS IS PRE-COMMITTED AND IS PRINTED BY `--plan` BEFORE ANY RESULT EXISTS.
Three results on this branch died because a conclusion was drawn from a measurement too
small to support it, and §6.4 records a fourth that was caught in draft. So the decision
rule, and the power the design actually has, are fixed in `plan()` below and are written
to the output directory before the first run starts.

    python paired_max_sugg.py --dir D --plan          # pre-registration, no compute
    python paired_max_sugg.py --dir D --make --n 50 --snr 14
    python paired_max_sugg.py --dir D --run           # 4 combinations, resumable
    python paired_max_sugg.py --dir D --analyse
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from itertools import product
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PILOT = HERE / "injection_pilot.py"
SRC = HERE.parent.parent / "src"
ARMS = ("aggressive", "quadrature")
BUFFERS = (2**14, 2**18)

# The campaign's effect of interest, from §9: pi = 0.60 at p_disc = 0.30, i.e. a
# per-pair excess of quadrature-only over aggressive-only discordances of
#   p_disc * (2*pi - 1) = 0.30 * 0.20 = 0.06.
# An arm-dependent buffer bias rides on exactly the same channel, so this is the scale
# at which the interaction has to be bounded for the campaign to be safe.
EFFECT_OF_INTEREST = 0.06


# --------------------------------------------------------------------------- plan


def plan(n: int) -> dict:
    """The pre-committed analysis, and an honest statement of what n buys.

    Contrast per realisation i and buffer m:

        d_i(m) = 1{quadrature recovered} - 1{aggressive recovered}   in {-1, 0, +1}
        Delta_i = d_i(2^18) - d_i(2^14)                              in {-2, ..., +2}

    If the buffer acts on both arms alike, raising it moves both arms' recoveries the
    same way and the *contrast* is unchanged in distribution, so Delta is symmetric
    about 0. An arm-dependent bias breaks that symmetry in a definite direction.

    PRIMARY (the screen): two-sided exact binomial (sign test) on sign(Delta_i) over the
    realisations with Delta_i != 0, H0: p = 1/2, alpha = 0.05.

    SECONDARY, descriptive, not tested: per-arm buffer flip counts (gained/lost on
    raising the buffer, with exact CIs), the four recovery counts, and the saturation
    distribution per cell -- the last being the mechanism, not the effect.
    """
    # What a sign test on n realisations can resolve. The fraction of realisations with
    # Delta != 0 is not known in advance; the pilot's only complete strand (aggressive,
    # 20 realisations, 2 flips) suggests it is small, so this is scanned.
    rows = []
    for frac in (0.10, 0.20, 0.40):
        m = frac * n
        # sd of the mean of Delta over n realisations, taking |Delta| = 1 when nonzero
        sd_mean = np.sqrt(frac / n)
        rows.append({
            "frac_nonzero": frac,
            "expected_m": m,
            "sd_of_mean_Delta": sd_mean,
            "half_width_95": 1.96 * sd_mean,
        })
    return {
        "n_realisations": n,
        "primary": "two-sided exact binomial on sign(Delta_i), Delta_i != 0, H0 p=1/2",
        "alpha": 0.05,
        "effect_of_interest_per_pair": EFFECT_OF_INTEREST,
        "resolution": rows,
        "interpretation": {
            "reject": (
                "DECISIVE AGAINST THE CAMPAIGN as specified. The buffer bias is "
                "arm-dependent, so no comparison at a buffer either arm saturates "
                "measures tiling. §7.1 (raise max_sugg until the 99th percentile of "
                "saturation is below 0.9 in both arms) becomes mandatory, not advisory."
            ),
            "fail_to_reject": (
                "NOT a green light, and must not be reported as one. n=50 bounds the "
                "per-pair interaction to roughly +-0.13 at 95% (see `resolution`), "
                "which is about twice the 0.06 effect the campaign is powered for. So "
                "a null here excludes a gross arm-dependence and nothing more; the "
                "remedy in §7.1 is unchanged by it."
            ),
        },
    }


# ---------------------------------------------------------------------------- run


def _env() -> dict:
    """Force `pyloki` to resolve to THIS worktree's `src`, not the installed one.

    The venv's editable install points at the main checkout, which is a different
    branch. Here that fails loudly (`pyloki.core.metric` does not exist on `main`), but
    a module that exists on both branches with different contents would be picked up
    silently, and a campaign comparing two tilings under the wrong search code would
    look exactly like a campaign comparing two tilings. So it is pinned, and checked.
    """
    import os
    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join(
        [str(SRC), *(p for p in (env.get("PYTHONPATH", "").split(os.pathsep)) if p)])
    resolved = subprocess.run(
        [sys.executable, "-c", "import pyloki; print(pyloki.__file__)"],
        env=env, capture_output=True, text=True, check=True).stdout.strip()
    if not resolved.startswith(str(SRC)):
        raise SystemExit(f"pyloki resolves to {resolved}, expected under {SRC}")
    return env


def out_path(outdir: Path, arm: str, buf: int) -> Path:
    return outdir / f"arm_{arm}_{buf}.json"


def run(outdir: Path) -> None:
    """Drive the four combinations, one process each. Resumable: completed cells skip."""
    env = _env()
    for arm, buf in product(ARMS, BUFFERS):
        dest = out_path(outdir, arm, buf)
        if dest.exists():
            print(f"skip {arm} {buf} (exists)", flush=True)
            continue
        print(f"=== {arm} max_sugg={buf} ===", flush=True)
        subprocess.run(
            [sys.executable, str(PILOT), "--dir", str(outdir), "--arm", arm,
             "--max-sugg", str(buf)],
            env=env, check=True)
        # injection_pilot.py always writes arm_{arm}.json; rename so the four cells
        # coexist. This is the renaming the killed driver did by hand (§11.6).
        (outdir / f"arm_{arm}.json").rename(dest)
        print(f"wrote {dest.name}", flush=True)


# ------------------------------------------------------------------------ analyse


def _binom_two_sided(k: int, n: int) -> float:
    """Exact two-sided binomial p-value against p = 1/2, by the symmetry of the null."""
    if n == 0:
        return float("nan")
    from math import comb
    tail = sum(comb(n, i) for i in range(n + 1)
               if abs(i - n / 2) >= abs(k - n / 2) - 1e-12)
    return min(1.0, tail / 2**n)


def _wilson(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    if n == 0:
        return 0.0, 1.0
    p = k / n
    d = 1 + z**2 / n
    c = (p + z**2 / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / d
    return max(c - h, 0.0), min(c + h, 1.0)


def analyse(outdir: Path) -> dict:
    cells, files = {}, None
    for arm, buf in product(ARMS, BUFFERS):
        p = out_path(outdir, arm, buf)
        if not p.exists():
            raise SystemExit(f"missing cell {p.name}; run --run first")
        rows = json.loads(p.read_text())
        names = [r["file"] for r in rows]
        if files is None:
            files = names
        elif names != files:
            raise SystemExit(f"pairing broken in {p.name}")
        cells[arm, buf] = rows

    n = len(files)
    rec = {k: np.array([r["recovered"] for r in v]) for k, v in cells.items()}
    lo, hi = BUFFERS

    # primary: the arm contrast, and how the buffer moves it
    d_lo = rec["quadrature", lo].astype(int) - rec["aggressive", lo].astype(int)
    d_hi = rec["quadrature", hi].astype(int) - rec["aggressive", hi].astype(int)
    delta = d_hi - d_lo
    nz = delta[delta != 0]
    k_pos = int((nz > 0).sum())
    m = int(nz.size)
    p_value = _binom_two_sided(k_pos, m)
    mean_delta = float(delta.mean())
    se_delta = float(delta.std(ddof=1) / np.sqrt(n)) if n > 1 else float("nan")

    # secondary: per-arm buffer flips
    flips = {}
    for arm in ARMS:
        a, b = rec[arm, lo], rec[arm, hi]
        gained, lost = int((~a & b).sum()), int((a & ~b).sum())
        flips[arm] = {
            "recovered_lo": int(a.sum()), "recovered_hi": int(b.sum()),
            "gained_raising_buffer": gained, "lost_raising_buffer": lost,
            "net_per_realisation": (gained - lost) / n,
            "mcnemar_p": _binom_two_sided(gained, gained + lost),
        }

    sat = {f"{arm}_{buf}": {
        "median": float(np.median([r["saturation"] for r in cells[arm, buf]])),
        "p99": float(np.percentile([r["saturation"] for r in cells[arm, buf]], 99)),
        "max": float(max(r["saturation"] for r in cells[arm, buf])),
        "frac_above_0.9": float(np.mean(
            [r["saturation"] > 0.9 for r in cells[arm, buf]])),
    } for arm, buf in product(ARMS, BUFFERS)}

    res = {
        "n_realisations": n,
        "recovered": {f"{arm}_{buf}": int(rec[arm, buf].sum())
                      for arm, buf in product(ARMS, BUFFERS)},
        "primary": {
            "statistic": "sign test on Delta_i = d_i(hi) - d_i(lo)",
            "m_nonzero": m, "k_positive": k_pos, "p_value": p_value,
            "ci95_positive_fraction": _wilson(k_pos, m),
            "mean_Delta": mean_delta, "se_mean_Delta": se_delta,
            "ci95_mean_Delta": [mean_delta - 1.96 * se_delta,
                                mean_delta + 1.96 * se_delta],
            "reject_at_0.05": bool(p_value < 0.05) if m else False,
        },
        # What n WOULD bound the interaction below the campaign's effect of interest.
        # The half-width goes as 1/sqrt(n), so n_needed = n * (halfwidth / effect)^2.
        # The formula is fixed here before any cell was analysed; only the measured
        # se that feeds it comes from the data. This exists because the cheap move on
        # seeing a null is to call it good enough, and the price of an experiment that
        # could actually clear the campaign belongs on the table next to it.
        "n_to_bound_interaction": (
            int(np.ceil(n * (1.96 * se_delta / EFFECT_OF_INTEREST) ** 2))
            if np.isfinite(se_delta) and se_delta > 0 else None),
        "per_arm_buffer_flips": flips,
        "saturation": sat,
        "discordance": {
            f"buffer_{b}": {
                "n01_quad_only": int((~rec["aggressive", b]
                                      & rec["quadrature", b]).sum()),
                "n10_agg_only": int((rec["aggressive", b]
                                     & ~rec["quadrature", b]).sum()),
            } for b in BUFFERS},
    }
    return res


def _print(res: dict) -> None:
    n = res["n_realisations"]
    print(f"\nrealisations           {n}")
    print("recovered              " + "  ".join(
        f"{k}={v}" for k, v in res["recovered"].items()))
    p = res["primary"]
    print(f"\nPRIMARY  sign test on Delta_i (does the buffer move the arm contrast?)")
    print(f"  nonzero Delta        {p['m_nonzero']}/{n}   "
          f"(positive = toward quadrature at the high buffer)")
    print(f"  k positive           {p['k_positive']}")
    print(f"  exact two-sided p    {p['p_value']:.4f}   "
          f"{'REJECT' if p['reject_at_0.05'] else 'no rejection'} at alpha=0.05")
    lo_, hi_ = p["ci95_mean_Delta"]
    print(f"  mean Delta           {p['mean_Delta']:+.3f}  "
          f"[95% {lo_:+.3f}, {hi_:+.3f}]   effect of interest {EFFECT_OF_INTEREST}")
    nb = res.get("n_to_bound_interaction")
    print(f"  n to bound it < {EFFECT_OF_INTEREST}   "
          f"{nb if nb else 'n/a'}   pairs per cell, i.e. this experiment "
          f"{(nb / n):.1f}x over" if nb else "")
    print("\nper-arm buffer flips (raising 2^14 -> 2^18)")
    for arm, f in res["per_arm_buffer_flips"].items():
        print(f"  {arm:11s} {f['recovered_lo']:>3d} -> {f['recovered_hi']:<3d} "
              f"gained {f['gained_raising_buffer']}, lost {f['lost_raising_buffer']}, "
              f"net/real {f['net_per_realisation']:+.3f}, "
              f"McNemar p={f['mcnemar_p']:.3f}")
    print("\nsaturation (ncand / max_sugg)")
    for k, s in res["saturation"].items():
        print(f"  {k:22s} median {s['median']:.3f}  p99 {s['p99']:.3f}  "
              f"max {s['max']:.3f}  frac>0.9 {s['frac_above_0.9']:.2f}")
    print("\ndiscordance by buffer")
    for k, d in res["discordance"].items():
        print(f"  {k:14s} n01(quad-only)={d['n01_quad_only']}  "
              f"n10(agg-only)={d['n10_agg_only']}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", type=Path, required=True)
    ap.add_argument("--plan", action="store_true")
    ap.add_argument("--make", action="store_true")
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--analyse", action="store_true")
    ap.add_argument("--n", type=int, default=50)
    ap.add_argument("--snr", type=float, default=14.0)
    args = ap.parse_args()

    if args.plan:
        pl = plan(args.n)
        args.dir.mkdir(parents=True, exist_ok=True)
        (args.dir / "preregistration.json").write_text(json.dumps(pl, indent=1))
        print(json.dumps(pl, indent=1))
    elif args.make:
        subprocess.run(
            [sys.executable, str(PILOT), "--dir", str(args.dir), "--make",
             "--n", str(args.n), "--snr", str(args.snr)], env=_env(), check=True)
    elif args.run:
        run(args.dir)
    elif args.analyse:
        res = analyse(args.dir)
        (args.dir / "paired_result.json").write_text(json.dumps(res, indent=1))
        _print(res)


if __name__ == "__main__":
    main()
