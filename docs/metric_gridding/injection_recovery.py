"""Phase 3 step 3: inject a signal, prune, and ask whether it was recovered.

Recovery is judged in the **metric**, not in raw parameter distance: a candidate counts
as a recovery when its mismatch to the injected signal, in the metric of the full
accumulated baseline, is below `m_recover`. That is the only criterion that means the
same thing for both strategies -- a fixed tolerance in Hz would favour whichever grid
happens to be finer in frequency.
"""

from __future__ import annotations

import tempfile
import time
from dataclasses import dataclass

import numpy as np

from pyloki.core import metric
from pyloki.ffa import DynamicProgramming
from pyloki.periodogram import ScatteredPeriodogram
from pyloki.prune import prune_dyp_tree
from pyloki.simulation.pulse import PulseSignalConfig
from pyloki.utils.misc import C_VAL


@dataclass
class Injection:
    """Truth, in the leaf's own axis order [d_kmax .. d_2, d_1] plus f0."""

    freq: float
    accel: float = 0.0
    jerk: float = 0.0
    snap: float = 0.0

    def mod_kwargs(self, poly_order: int) -> dict[str, float]:
        kw = {"acc": self.accel, "jerk": self.jerk}
        if poly_order >= 4:
            kw["snap"] = self.snap
        return kw

    def truth_vector(self, poly_order: int) -> np.ndarray:
        """[d_kmax .. d_2, d_1] with d_1 the velocity offset (0 at the epoch)."""
        order = [self.snap, self.jerk, self.accel][-(poly_order - 1):]
        return np.array([*order, 0.0])


@dataclass
class Result:
    strategy: str
    snr: float
    recovered: bool
    best_mismatch: float
    best_score: float
    n_candidates: int
    seconds: float


def run_one(
    cfg,
    injection: Injection,
    snr: float,
    *,
    ducy: float = 0.1,
    n_runs: int = 1,
    max_sugg: int = 2**12,
    m_recover: float = 1.0,
    thresholds: np.ndarray | None = None,
) -> Result:
    """One injection, one strategy, start to finish."""
    t0 = time.perf_counter()
    sig = PulseSignalConfig(
        period=1.0 / injection.freq, dt=cfg.tsamp, nsamps=cfg.nsamps,
        snr=snr, ducy=ducy, mod_kwargs=injection.mod_kwargs(cfg.prune_poly_order),
    )
    tim = sig.generate(shape="gaussian")

    dyp = DynamicProgramming(tim, cfg)
    dyp.initialize()
    dyp.execute()

    ref = dyp.nsegments // 2
    if thresholds is None:
        # Placeholder ramp. Phase 3 step 2 replaces this with the Viterbi-optimised
        # scheme, recalibrated per strategy -- B(s) differs, so the schemes must too.
        n_levels = len(
            cfg.generate_branching_pattern(kind="poly_taylor_moving", ref_seg=ref),
        )
        thresholds = np.linspace(1.5, 8.0, n_levels)

    with tempfile.TemporaryDirectory() as td:
        rf = prune_dyp_tree(
            dyp, thresholds, n_runs=n_runs, max_sugg=max_sugg, outdir=td,
            file_prefix="inj", poly_basis="taylor", n_workers=1,
            use_moving_grid=True,
        )
        pgram = ScatteredPeriodogram.load(str(rf))
        df = pgram.data

    # Mismatch of every candidate to the truth, in the full-baseline metric.
    po = cfg.prune_poly_order
    tobs = cfg.nsamps * cfg.tsamp
    g = metric.poly_phase_metric(
        0.0, -tobs / 2, tobs / 2, po, injection.freq, cfg.nbins, cfg.ducy_max,
    )
    truth = injection.truth_vector(po)
    names = cfg.param_names
    cand = np.column_stack([df[n].to_numpy() for n in names[:-1]])
    # The last search axis is frequency; convert to a velocity offset like the truth.
    vel = C_VAL * (1.0 - df[names[-1]].to_numpy() / injection.freq)
    cand = np.column_stack([cand, vel])

    delta = cand - truth
    mism = np.einsum("ni,ij,nj->n", delta, g, delta)
    best = int(np.argmin(mism)) if len(mism) else -1
    return Result(
        strategy=cfg.tiling_strategy,
        snr=snr,
        recovered=bool(len(mism) and mism[best] <= m_recover),
        best_mismatch=float(mism[best]) if len(mism) else float("inf"),
        best_score=float(df["score"].to_numpy()[best]) if len(mism) else 0.0,
        n_candidates=int(len(df)),
        seconds=time.perf_counter() - t0,
    )
