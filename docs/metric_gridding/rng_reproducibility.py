"""Reproducer: pyloki cannot reproduce a run, and `np.random.seed` does not help.

Self-contained — no branch code, no data files, no `PYTHONPATH` games beyond importing
the library under test. Runs in a few seconds.

Three checks:

1. `DynamicThresholdScheme` built with identical arguments produces different threshold
   ladders. This is the site that matters most, because a ladder is an *input* to a
   search: two searches configured identically are not running the same calibration.
2. Seeding the legacy global RNG does not fix it. `np.random.default_rng()` ignores
   `np.random.seed`, which is what makes this defect easy to miss.
3. The unseeded constructions are enumerated from the installed source, so the list in
   the write-up cannot drift from the code.

Usage:

    PYTHONPATH=<worktree>/src <repo>/.venv/bin/python docs/metric_gridding/rng_reproducibility.py
"""

from __future__ import annotations

import logging
import pathlib
import re
import tempfile
import warnings

logging.disable(logging.CRITICAL)
warnings.filterwarnings("ignore")

import numpy as np  # noqa: E402

import pyloki  # noqa: E402
from pyloki.detection import schemes  # noqa: E402
from pyloki.detection.thresholding import DynamicThresholdScheme  # noqa: E402

# A deliberately small scheme: 8 stages, B(s) = 2. Nothing here depends on the
# configuration -- it is chosen only so the whole script runs in seconds.
SCHEME_KW = {
    "branching_pattern": np.array([2] * 8),
    "ref_ducy": 0.1,
    "nbins": 32,
    "ntrials": 1024,
    "nprobs": 10,
    "nthresholds": 100,
    "snr_final": 8,
}
N_REPEATS = 20
P_D = 0.1


def one_ladder() -> np.ndarray:
    """Build one scheme and backtrack a ladder at P_D, exactly as a user would."""
    dyn = DynamicThresholdScheme(**SCHEME_KW)
    dyn.run(thres_neigh=11)
    with tempfile.TemporaryDirectory() as tmp:
        analyser = schemes.DynamicThresholdSchemeAnalyser.from_file(dyn.save(outdir=tmp))
        best = analyser.backtrack_best(min_probs=[P_D])[0]
    return np.array([entry.threshold for entry in best.entries])


def report(label: str, *, legacy_seed: bool) -> None:
    rows = []
    for _ in range(N_REPEATS):
        if legacy_seed:
            np.random.seed(42)  # noqa: NPY002 -- the point is that this does nothing
        rows.append(one_ladder())
    arr = np.asarray(rows)
    distinct = len({r.tobytes() for r in arr})
    spread = arr.max(axis=0) - arr.min(axis=0)
    print(f"  {label}")
    print(f"      distinct ladders out of {N_REPEATS} : {distinct}")
    print(f"      max per-stage spread             : {spread.max():.3f}")
    print(f"      mean per-stage spread            : {spread.mean():.3f}")
    print(f"      worst stage                      : {int(spread.argmax())}")


def unseeded_sites() -> list[str]:
    """Enumerate `np.random.default_rng()` with no argument, from the installed source."""
    root = pathlib.Path(pyloki.__file__).parent
    pattern = re.compile(r"default_rng\(\s*\)")
    hits = []
    for path in sorted(root.rglob("*.py")):
        for lineno, line in enumerate(path.read_text().splitlines(), 1):
            if pattern.search(line):
                hits.append(f"{path.relative_to(root)}:{lineno}  {line.strip()}")
    return hits


def main() -> None:
    print(f"library under test: {pyloki.__file__}\n")

    print("1/2. DynamicThresholdScheme ladders, identical arguments every time")
    report("unseeded (as shipped)", legacy_seed=False)
    report("after np.random.seed(42) before each construction", legacy_seed=True)

    print("\n3. Unseeded default_rng() constructions in the installed source")
    sites = unseeded_sites()
    for site in sites:
        print(f"      {site}")
    print(f"\n      {len(sites)} sites")


if __name__ == "__main__":
    main()
