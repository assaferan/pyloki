"""Every RNG in the library must be steerable from its caller.

The library used to call `np.random.default_rng()` with no argument in eight
places and expose no way to seed any of them, so a run could not be repeated.
`np.random.seed` does not help: `default_rng` ignores the legacy global seed,
which is what made the defect easy to miss -- the usual reflex appears to work
and changes nothing.

Each public entry point now takes `seed`, accepting an int, a
`np.random.Generator`, or `None`. **`None` is still the default and still draws
fresh entropy**, so these tests also pin that the fix did not quietly make the
library deterministic for callers who never asked for it.

`test_no_bare_default_rng_in_library` enumerates the call sites from the source
rather than from a list typed here, so a new unseeded generator fails this file
instead of going unnoticed.

**Known limitation.** A seed does not make `DynamicThresholdScheme.run()`
reproducible under parallel execution -- see
`TestThresholding.test_run_reproduces_single_threaded`. Every other seeded path
here reproduces on 14 threads as well as on one.
"""

from __future__ import annotations

import ast
import pathlib
from typing import TYPE_CHECKING

import numba
import numpy as np
import pytest

import pyloki
from pyloki.detection import thresholding

# aliased: pytest would otherwise try to collect the library's `Test*` class
from pyloki.sensitivity.sim_ffa import TestFFASensitivity as FFASensitivitySim
from pyloki.simulation.pulse import PulseSignalConfig

if TYPE_CHECKING:
    from pyloki.detection.schemes import StatesInfo

# Small enough that the whole file runs in a few seconds; nothing here depends
# on the configuration, only on whether two runs of it agree.
NSAMPS = 2**16
BRANCHING_PATTERN = np.full(4, 2.0)
SCHEME_KW = {
    "branching_pattern": BRANCHING_PATTERN,
    "ref_ducy": 0.1,
    "nbins": 64,
    "ntrials": 256,
    "snr_final": 8.0,
}
GENERATE_METHODS = ("generate", "generate_simple", "generate_noise", "generate_old")


def _config(seed: int | np.random.Generator | None) -> PulseSignalConfig:
    return PulseSignalConfig(
        period=0.01,
        dt=64e-6,
        nsamps=NSAMPS,
        snr=20,
        ducy=0.1,
        mod_kwargs={"acc": 100.0},
        seed=seed,
    )


def _survival(info: StatesInfo) -> np.ndarray:
    return info.get_info("success_h1_cumul")


class TestPulseSignalConfig:
    """`simulation/pulse.py` -- the four noise draws."""

    @pytest.mark.parametrize("method", GENERATE_METHODS)
    def test_same_seed_reproduces(self, method: str) -> None:
        first = getattr(_config(42), method)().ts_e
        second = getattr(_config(42), method)().ts_e
        np.testing.assert_array_equal(
            first, second, err_msg=f"{method}() is not reproducible under seed=42",
        )

    @pytest.mark.parametrize("method", GENERATE_METHODS)
    def test_different_seeds_differ(self, method: str) -> None:
        assert not np.array_equal(
            getattr(_config(42), method)().ts_e, getattr(_config(43), method)().ts_e,
        ), f"{method}() ignores the seed -- 42 and 43 gave the same noise"

    @pytest.mark.parametrize("method", GENERATE_METHODS)
    def test_unseeded_still_draws_fresh_entropy(self, method: str) -> None:
        """The default must not have become deterministic."""
        assert not np.array_equal(
            getattr(_config(None), method)().ts_e,
            getattr(_config(None), method)().ts_e,
        ), f"{method}() is deterministic without a seed; default behaviour changed"

    def test_successive_calls_advance_the_stream(self) -> None:
        """A seed fixes the *sequence*, not every call.

        Otherwise a loop over one seeded config would silently draw the same
        realisation every iteration, which is the trap in re-seeding per call.
        """
        cfg = _config(42)
        first, second = cfg.generate_noise().ts_e, cfg.generate_noise().ts_e
        assert not np.array_equal(first, second), (
            "two successive draws from one seeded config are identical; "
            "the generator is being re-created per call"
        )
        replay = _config(42)
        np.testing.assert_array_equal(replay.generate_noise().ts_e, first)
        np.testing.assert_array_equal(replay.generate_noise().ts_e, second)

    def test_accepts_a_generator(self) -> None:
        first = _config(np.random.default_rng(9)).generate_noise().ts_e
        second = _config(np.random.default_rng(9)).generate_noise().ts_e
        np.testing.assert_array_equal(first, second)

    def test_seed_survives_get_updated(self) -> None:
        """`get_updated` rebuilds the config through the constructor."""
        updated = _config(42).get_updated({"snr": 30})
        assert updated.seed == 42
        assert updated.snr == 30

    def test_stored_generator_does_not_break_equality(self) -> None:
        """Generators compare by identity; that must not leak into config `==`."""
        assert _config(42) == _config(42)
        assert _config(None) == _config(None)


class TestThresholding:
    """`detection/thresholding.py` -- the three ladder generators."""

    def test_determine_scheme_same_seed_reproduces(self) -> None:
        first = thresholding.determine_scheme(
            1.0 / BRANCHING_PATTERN, seed=7, **SCHEME_KW,
        )
        second = thresholding.determine_scheme(
            1.0 / BRANCHING_PATTERN, seed=7, **SCHEME_KW,
        )
        np.testing.assert_array_equal(first.thresholds, second.thresholds)

    def test_determine_scheme_different_seeds_differ(self) -> None:
        first = thresholding.determine_scheme(
            1.0 / BRANCHING_PATTERN, seed=7, **SCHEME_KW,
        )
        second = thresholding.determine_scheme(
            1.0 / BRANCHING_PATTERN, seed=8, **SCHEME_KW,
        )
        assert not np.array_equal(first.thresholds, second.thresholds)

    def test_determine_scheme_unseeded_still_varies(self) -> None:
        first = thresholding.determine_scheme(1.0 / BRANCHING_PATTERN, **SCHEME_KW)
        second = thresholding.determine_scheme(1.0 / BRANCHING_PATTERN, **SCHEME_KW)
        assert not np.array_equal(first.thresholds, second.thresholds), (
            "unseeded ladders are now identical; default behaviour changed"
        )

    def test_evaluate_scheme_same_seed_reproduces(self) -> None:
        ladder = np.asarray(
            thresholding.determine_scheme(
                1.0 / BRANCHING_PATTERN, seed=7, **SCHEME_KW,
            ).thresholds,
            dtype=np.float64,
        )
        first = thresholding.evaluate_scheme(ladder, seed=3, **SCHEME_KW)
        second = thresholding.evaluate_scheme(ladder, seed=3, **SCHEME_KW)
        np.testing.assert_array_equal(_survival(first), _survival(second))

    def test_dynamic_threshold_scheme_seeds_its_generator(self) -> None:
        """The seed reaches `self.rng`. That is NOT the same as `run()` reproducing."""

        def draws(seed: int | None) -> np.ndarray:
            scheme = thresholding.DynamicThresholdScheme(
                nthresholds=20, seed=seed, **SCHEME_KW,
            )
            return scheme.rng.standard_normal(8)

        np.testing.assert_array_equal(draws(5), draws(5))
        assert not np.array_equal(draws(None), draws(None))

    def test_run_reproduces_single_threaded(self) -> None:
        """`run()` reproduces under a seed **only** when numba is single-threaded.

        `run_stage_legacy` is `@njit(parallel=True)` and draws from the shared
        generator inside a `prange`, so which thread takes which draw depends on
        scheduling. A seed fixes the stream, not the order it is consumed in, and the
        ladder therefore still varies under parallel execution -- measured at 2
        distinct ladders in 4 runs on 14 threads, against 1 in 4 on one thread.

        Seeding `__init__` is necessary but not sufficient at this site. Closing it
        needs per-iteration generators inside the kernel, which is a separate change.
        The single-threaded case is pinned because it is the part the seed does fix.
        The parallel case is deliberately not asserted: a test asserting that two runs
        *differ* would itself be flaky.
        """
        original = numba.get_num_threads()
        numba.set_num_threads(1)
        try:
            states = []
            for _ in range(2):
                scheme = thresholding.DynamicThresholdScheme(
                    np.array([4.0, 2.0, 2.0]),
                    ref_ducy=0.1,
                    nbins=32,
                    ntrials=256,
                    nprobs=8,
                    nthresholds=20,
                    snr_final=8.0,
                    seed=11,
                )
                scheme.run(thres_neigh=5)
                # The whole state record, not `threshold` alone: thresholds come off
                # a fixed linspace and are identical even unseeded, so comparing them
                # would pass vacuously. `success_h0` is the field that carries the
                # noise -- it differs between two unseeded single-threaded runs.
                states.append(scheme.states.copy())
        finally:
            numba.set_num_threads(original)
        np.testing.assert_array_equal(
            states[0]["success_h0"],
            states[1]["success_h0"],
            err_msg="seeded run() is not reproducible even on a single thread",
        )


class TestSimFFA:
    """`sensitivity/sim_ffa.py`."""

    def test_same_seed_reproduces(self) -> None:
        cfg = PulseSignalConfig(period=0.01, dt=64e-6, nsamps=NSAMPS, snr=20, ducy=0.1)
        limits = np.array([[100.0, 200.0], [-10.0, 10.0]])

        def draws(seed: int | None) -> np.ndarray:
            sim = FFASensitivitySim(cfg=cfg, param_limits=limits, seed=seed)
            return sim.rng.standard_normal(8)

        np.testing.assert_array_equal(draws(5), draws(5))
        assert not np.array_equal(draws(None), draws(None))


def test_no_bare_default_rng_in_library() -> None:
    """No `default_rng()` anywhere in `src/` may be called with no argument.

    Enumerated from the source so that a newly added unseeded generator fails
    here rather than reintroducing the defect unnoticed.
    """
    root = pathlib.Path(pyloki.__file__).parent
    offenders = []
    for path in sorted(root.rglob("*.py")):
        tree = ast.parse(path.read_text(), filename=str(path))
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            func = node.func
            if (
                isinstance(func, ast.Attribute)
                and func.attr == "default_rng"
                and not node.args
                and not node.keywords
            ):
                offenders.append(f"{path.relative_to(root)}:{node.lineno}")
    assert not offenders, (
        "unseeded np.random.default_rng() call(s) -- the caller cannot reproduce "
        f"this run: {', '.join(offenders)}"
    )
