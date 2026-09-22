"""Reproducer for the unseeded-RNG defect, and the check that the fix closes it.

Companion to `FINDINGS_rng_seeding.md`. Originally written on the `injection-design`
branch to demonstrate the defect; extended here to also demonstrate the fix, so one
script covers the whole story.

Self-contained -- no data files, no fixtures. Runs in well under a minute.

Five checks:

1. `DynamicThresholdScheme` built with identical arguments and **no seed** produces
   different threshold ladders. That is the shipped behaviour and it must survive the
   fix: seeding was added as an option, not as a new default.
2. Seeding the legacy global RNG does not change that. `np.random.default_rng()`
   ignores `np.random.seed`, which is what makes this defect easy to miss -- the usual
   reflex appears to work and changes nothing.
3. The same construction **with a seed**, on the default thread count: the ladders
   **still differ**. This is the limitation, not the fix. `run_stage_legacy` is
   `@njit(parallel=True)` and draws from the shared generator inside a `prange`, so a
   seed fixes the stream but not the order threads consume it in.
4. The same, **single-threaded**: now identical. This is the part the seed does fix,
   and it locates the residual non-determinism precisely -- in the thread scheduling,
   not in the seeding.
5. Every `default_rng()` construction is enumerated from the installed source and
   reported as seeded or unseeded, so the site list in the write-up cannot drift from
   the code.

A threshold ladder is an *input* to a search, not an output of one, which is why this
site is the one worth demonstrating: two searches a user believes are identically
calibrated were not.

Note what checks 3 and 4 mean together: seeding `DynamicThresholdScheme.__init__` is
necessary but **not sufficient** to reproduce a ladder under parallel execution. Every
other seeded path in the library -- `determine_scheme`, `evaluate_scheme` and all four
`PulseSignalConfig.generate*` methods -- reproduces on 14 threads as well as on one;
this one site does not. Closing it needs per-iteration generators inside the kernel.

Usage:

    PYTHONPATH=<worktree>/src <repo>/.venv/bin/python rng_reproducibility.py

Against a tree without the fix, checks 3 and 4 fail with a TypeError on the `seed`
keyword, and check 5 reports eight unseeded sites.
"""

from __future__ import annotations

import ast
import logging
import pathlib
import tempfile
import warnings

logging.disable(logging.CRITICAL)
warnings.filterwarnings("ignore")

import numba  # noqa: E402
import numpy as np  # noqa: E402

import pyloki  # noqa: E402
from pyloki.detection import schemes  # noqa: E402
from pyloki.detection.thresholding import DynamicThresholdScheme  # noqa: E402

# A deliberately small scheme: 8 stages, B(s) = 2. Nothing here depends on the
# configuration -- it is chosen only so the whole script runs quickly.
SCHEME_KW = {
    "branching_pattern": np.array([2] * 8),
    "ref_ducy": 0.1,
    "nbins": 32,
    "ntrials": 1024,
    "nprobs": 10,
    "nthresholds": 100,
    "snr_final": 8,
}
N_UNSEEDED = 10
N_SEEDED = 5  # identity needs far fewer repeats than variation does
N_SERIAL = 3  # single-threaded runs are slower; 3 is enough to show identity
P_D = 0.1


def one_ladder(seed: int | None = None) -> np.ndarray:
    """Build one scheme and backtrack a ladder at P_D, exactly as a user would."""
    kwargs = dict(SCHEME_KW)
    if seed is not None:
        kwargs["seed"] = seed
    dyn = DynamicThresholdScheme(**kwargs)
    dyn.run(thres_neigh=11)
    with tempfile.TemporaryDirectory() as tmp:
        saved = dyn.save(outdir=tmp)
        analyser = schemes.DynamicThresholdSchemeAnalyser.from_file(saved)
        best = analyser.backtrack_best(min_probs=[P_D])[0]
    return np.array([entry.threshold for entry in best.entries])


def report(label: str, ladders: list[np.ndarray], *, want_identical: bool) -> bool:
    arr = np.asarray(ladders)
    distinct = len({r.tobytes() for r in arr})
    spread = arr.max(axis=0) - arr.min(axis=0)
    ok = (distinct == 1) if want_identical else (distinct > 1)
    print(f"  {label}")
    print(f"      distinct ladders out of {len(ladders):<2}          : {distinct}")
    print(f"      max per-stage spread               : {spread.max():.3f}")
    print(f"      mean per-stage spread              : {spread.mean():.3f}")
    want = "identical" if want_identical else "varying"
    verdict = "OK" if ok else "UNEXPECTED"
    print(f"      expected {want:<9}            : {verdict}")
    return ok


def default_rng_sites() -> tuple[list[str], list[str]]:
    """Enumerate every `default_rng(...)` in the installed source, via the AST.

    Returns (unseeded, seeded). Parsed rather than grepped so that a call split over
    several lines is classified correctly.
    """
    root = pathlib.Path(pyloki.__file__).parent
    unseeded, seeded = [], []
    for path in sorted(root.rglob("*.py")):
        tree = ast.parse(path.read_text(), filename=str(path))
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            func = node.func
            if not (isinstance(func, ast.Attribute) and func.attr == "default_rng"):
                continue
            where = f"{path.relative_to(root)}:{node.lineno}"
            if node.args or node.keywords:
                seeded.append(f"{where}  {ast.unparse(node)}")
            else:
                unseeded.append(f"{where}  {ast.unparse(node)}")
    return unseeded, seeded


def main() -> None:
    print(f"library under test: {pyloki.__file__}\n")
    ok = True

    print(f"1. DynamicThresholdScheme, identical arguments, no seed ({N_UNSEEDED}x)")
    ok &= report(
        "unseeded -- must still vary, the default is unchanged",
        [one_ladder() for _ in range(N_UNSEEDED)],
        want_identical=False,
    )

    print("\n2. The same, after np.random.seed(42) before each construction")
    ladders = []
    for _ in range(N_UNSEEDED):
        np.random.seed(42)  # noqa: NPY002 -- the point is that this does nothing
        ladders.append(one_ladder())
    ok &= report(
        "legacy global seed -- does NOT help, default_rng ignores it",
        ladders,
        want_identical=False,
    )

    print(f"\n3. The same, with seed=42, on {numba.get_num_threads()} threads "
          f"({N_SEEDED}x)")
    ok &= report(
        "seeded but parallel -- STILL VARIES, this is the limitation",
        [one_ladder(seed=42) for _ in range(N_SEEDED)],
        want_identical=False,
    )

    print(f"\n4. The same, with seed=42, forced single-threaded ({N_SERIAL}x)")
    original = numba.get_num_threads()
    numba.set_num_threads(1)
    try:
        serial = [one_ladder(seed=42) for _ in range(N_SERIAL)]
    finally:
        numba.set_num_threads(original)
    ok &= report(
        "seeded and serial -- identical, so the residue is thread order",
        serial,
        want_identical=True,
    )

    print("\n5. default_rng() constructions in the installed source")
    unseeded, seeded = default_rng_sites()
    for site in unseeded:
        print(f"      UNSEEDED  {site}")
    for site in seeded:
        print(f"      seeded    {site}")
    print(f"\n      {len(unseeded)} unseeded, {len(seeded)} seeded")
    ok &= not unseeded

    print(f"\n{'ALL CHECKS AS EXPECTED' if ok else 'SOMETHING IS NOT AS EXPECTED'}")


if __name__ == "__main__":
    main()
