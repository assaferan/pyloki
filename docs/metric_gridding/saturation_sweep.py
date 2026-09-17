"""Stage 1 of the approved path: find a `max_sugg` at which the buffer stops binding.

§7.1 requires the 99th percentile of `saturation = ncand / max_sugg` to be below 0.9
**in both arms**. §6.5 measured, on 50 realisations at 2^18:

    aggressive  p99 0.915   median 0.165   median ncand  43 314
    quadrature  p99 0.975   median 0.724   median ncand 189 828

so BOTH arms fail, `aggressive` only just and `quadrature` badly. Until they clear it,
the realised cut is `max(scheme threshold, top-K, median)` and the campaign is not
measuring the ladder it claims to.

**The thing this sweep is really testing, and it may not have a happy answer.** Between
2^14 and 2^18 `aggressive`'s candidate count went 8 937 -> 43 314, a factor 4.8 for a
factor 16 of buffer: it is converging on a true count near 43k and its saturation fell
0.545 -> 0.165. `quadrature`'s went 10 703 -> 189 828, a factor **17.7** for the same
factor 16 — its saturation did not fall at all (0.653 -> 0.724). That is the signature of
a candidate count still growing as fast as the buffer, i.e. a search whose true output
size is nowhere near 2^18. If that continues, there is no feasible buffer at which
`quadrature` stops saturating, and the campaign as designed is impossible rather than
merely expensive. **Establishing that cheaply is a perfectly good outcome for this
sweep** — it is the same logic as §11.6: a design that cannot work is worth discovering
before 28 core-hours, not after.

So the sweep reports, per arm and buffer, the saturation percentiles AND the growth
exponent `d log2(ncand) / d log2(max_sugg)` between successive buffers. An exponent near
0 means converged; near 1 means the buffer is still the only thing setting the answer.

Cost scales roughly linearly with `max_sugg` (measured: `quadrature` 13x for 16x between
2^14 and 2^18), so this runs on a SUBSET of the realisations and escalates one buffer at
a time rather than committing to the whole ladder up front.

    python saturation_sweep.py --dir D --source S --n 10          # build the subset
    python saturation_sweep.py --dir D --buffer 524288            # one buffer, both arms
    python saturation_sweep.py --dir D --report
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PILOT = HERE / "injection_pilot.py"
SRC = HERE.parent.parent / "src"
ARMS = ("aggressive", "quadrature")
CRITERION = 0.9


def _env() -> dict:
    """Pin `pyloki` to this worktree's src; see paired_max_sugg._env."""
    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join(
        [str(SRC), *(p for p in env.get("PYTHONPATH", "").split(os.pathsep) if p)])
    got = subprocess.run(
        [sys.executable, "-c", "import pyloki; print(pyloki.__file__)"],
        env=env, capture_output=True, text=True, check=True).stdout.strip()
    if not got.startswith(str(SRC)):
        raise SystemExit(f"pyloki resolves to {got}, expected under {SRC}")
    return env


def build(outdir: Path, source: Path, n: int) -> None:
    """Symlink the first n realisations of an existing set, so the sweep is paired."""
    outdir.mkdir(parents=True, exist_ok=True)
    files = sorted(source.glob("tim_*.npz"))[:n]
    if not files:
        raise SystemExit(f"no realisations in {source}")
    for f in files:
        dest = outdir / f.name
        if not dest.exists():
            dest.symlink_to(f.resolve())
    (outdir / "meta.json").write_text(
        json.dumps({"n": len(files), "source": str(source)}))
    print(f"linked {len(files)} realisations from {source}")


def run(outdir: Path, buffer: int) -> None:
    env = _env()
    for arm in ARMS:
        dest = outdir / f"arm_{arm}_{buffer}.json"
        if dest.exists():
            print(f"skip {arm} {buffer} (exists)", flush=True)
            continue
        print(f"=== {arm} max_sugg={buffer} ===", flush=True)
        subprocess.run(
            [sys.executable, str(PILOT), "--dir", str(outdir), "--arm", arm,
             "--max-sugg", str(buffer)], env=env, check=True)
        (outdir / f"arm_{arm}.json").rename(dest)


def report(outdir: Path) -> dict:
    cells = {}
    for p in sorted(outdir.glob("arm_*_*.json")):
        stem = p.stem[len("arm_"):]
        arm, _, buf = stem.rpartition("_")
        if arm in ARMS:
            cells[arm, int(buf)] = json.loads(p.read_text())
    if not cells:
        raise SystemExit("no cells")

    rows = []
    for arm in ARMS:
        bufs = sorted(b for (a, b) in cells if a == arm)
        prev = None
        for b in bufs:
            r = cells[arm, b]
            sat = np.array([x["saturation"] for x in r])
            nc = np.array([x["ncand"] for x in r], dtype=float)
            secs = [x["seconds"] for x in r]
            row = {
                "arm": arm, "max_sugg": b, "log2_max_sugg": int(np.log2(b)),
                "n": len(r),
                "median_saturation": float(np.median(sat)),
                "p99_saturation": float(np.percentile(sat, 99)),
                "max_saturation": float(sat.max()),
                "median_ncand": float(np.median(nc)),
                "recovered": int(sum(x["recovered"] for x in r)),
                "median_seconds": float(np.median(secs[1:] or secs)),
                "clears_criterion": bool(np.percentile(sat, 99) < CRITERION),
            }
            if prev is not None:
                d_nc = np.log2(row["median_ncand"] / prev["median_ncand"])
                d_b = row["log2_max_sugg"] - prev["log2_max_sugg"]
                row["growth_exponent"] = float(d_nc / d_b) if d_b else None
            rows.append(row)
            prev = row

    print(f"\n{'arm':<11}{'max_sugg':>10}{'med sat':>9}{'p99 sat':>9}"
          f"{'med ncand':>12}{'growth':>8}{'s/run':>8}{'clears<0.9':>12}")
    for r in rows:
        g = r.get("growth_exponent")
        print(f"{r['arm']:<11}2^{r['log2_max_sugg']:<8d}"
              f"{r['median_saturation']:>9.3f}{r['p99_saturation']:>9.3f}"
              f"{r['median_ncand']:>12.0f}"
              f"{(f'{g:.2f}' if g is not None else '-'):>8}"
              f"{r['median_seconds']:>8.1f}"
              f"{('YES' if r['clears_criterion'] else 'no'):>12}")
    print("\ngrowth = d log2(median ncand) / d log2(max_sugg) vs the previous buffer.")
    print("  ~0 => the candidate count has converged and the buffer no longer binds.")
    print("  ~1 => the buffer is still the only thing setting the answer.")
    out = {"criterion": CRITERION, "rows": rows}
    (outdir / "sweep_report.json").write_text(json.dumps(out, indent=1))
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", type=Path, required=True)
    ap.add_argument("--source", type=Path)
    ap.add_argument("--n", type=int, default=10)
    ap.add_argument("--buffer", type=int)
    ap.add_argument("--report", action="store_true")
    a = ap.parse_args()
    if a.source:
        build(a.dir, a.source, a.n)
    elif a.buffer:
        run(a.dir, a.buffer)
    elif a.report:
        report(a.dir)


if __name__ == "__main__":
    main()
