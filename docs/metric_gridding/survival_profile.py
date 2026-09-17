"""Per-stage survival profile: is a covering leaf still alive at each prune level?

This is the measurement that D99 says nothing on either branch makes, and that both
open questions turn on. Two results rest on a *modelled* per-stage survival curve from
the single-leaf Gaussian family falsified in D94:

- my stage-split framing (D70) weights the geometric gain by "where detection is
  decided", and that window came from the falsified curve;
- the campaign's pair count spans 340-2840 (D98) depending on the stage weighting, which
  is the same curve.

So instead of modelling survival, record it. The pruning loop already keeps the
surviving candidates in `world_tree`, and `Pruning.execute()` is a loop over
`execute_iter()`, so the loop can be driven directly and the survivors inspected between
levels. Nothing in `src/` is modified -- other sessions depend on it.

At each level, leaves and truth are both expressed in that level's own frame (reference
`ref_cur`, half-span `t_half`) and the phase residual is evaluated over the accumulated
data, `tau` in `[-t_half, +t_half]`. What comes out per level is not a binary but the
distribution that a binary would be derived from: the minimum phase excursion to the
truth over all survivors, and how many survivors sit within a tolerance of it.
"""

from __future__ import annotations

import sys
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from injection_recovery import Injection  # noqa: E402

from pyloki.ffa import DynamicProgramming  # noqa: E402
from pyloki.prune import Pruning  # noqa: E402
from pyloki.simulation.pulse import PulseSignalConfig  # noqa: E402
from pyloki.utils import transforms  # noqa: E402
from pyloki.utils.misc import C_VAL  # noqa: E402
from sensitivity_loss import phase_excursion  # noqa: E402


@dataclass
class LevelRecord:
    """What survived at one prune level."""

    level: int
    n_leaves: int
    min_excursion: float          # nearest survivor to the truth, in eta/N_b
    n_within_1: int               # survivors within one tolerance
    n_within_3: int
    score_max: float
    t_half: float


@dataclass
class ProfileResult:
    strategy: str
    snr: float
    levels: list[LevelRecord] = field(default_factory=list)

    def alive(self, tol: float = 1.0) -> np.ndarray:
        """Boolean per level: was any survivor within `tol` of the truth?"""
        return np.array([r.min_excursion <= tol for r in self.levels])


def profile_one(
    cfg,
    injection: Injection,
    snr: float,
    thresholds: np.ndarray,
    *,
    ducy: float = 0.10,
    max_sugg: int = 2**18,
    ref_seg: int | None = None,
    seed: int | None = None,
) -> ProfileResult:
    """Run one real pruning search and record the survivor set at every level."""
    sig = PulseSignalConfig(
        period=1.0 / injection.freq, dt=cfg.tsamp, nsamps=cfg.nsamps,
        snr=snr, ducy=ducy,
        mod_kwargs=injection.mod_kwargs(cfg.prune_poly_order),
    )
    tim = sig.generate(shape="gaussian")

    dyp = DynamicProgramming(tim, cfg)
    dyp.initialize()
    dyp.execute()

    po = cfg.prune_poly_order
    tobs = cfg.nsamps * cfg.tsamp
    # Truth about the observation centre, in the leaf's own layout plus d_0.
    truth_c = np.zeros(po + 1)
    truth_c[:po] = injection.truth_vector(po)

    if ref_seg is None:
        ref_seg = dyp.nsegments // 2

    out = ProfileResult(strategy=cfg.tiling_strategy, snr=snr)
    with tempfile.TemporaryDirectory() as td:
        log_file = Path(td) / "log.txt"
        log_file.touch()
        prn = Pruning(
            dyp, thresholds, max_sugg=max_sugg, batch_size=1024,
            poly_basis="taylor", use_moving_grid=True,
        )
        prn.initialize(ref_seg, log_file)
        for _ in range(dyp.nsegments - 1):
            prn.execute_iter(log_file)
            lvl = prn.prune_level
            ref_cur, t_half = prn.scheme.get_current_coord(lvl, True)

            leaves = np.asarray(prn.world_tree.leaves)
            n = int(prn.world_tree.size)
            if n == 0 or leaves.shape[0] == 0:
                out.levels.append(LevelRecord(lvl, 0, np.inf, 0, 0, 0.0, t_half))
                continue
            leaves = leaves[:n]

            # Truth into this level's frame, then residual over the accumulated data.
            truth_l = transforms.shift_taylor_params(truth_c, ref_cur - tobs / 2)
            f0 = float(leaves[0, -1, 0])
            tol = cfg.eta / cfg.nbins
            tau = np.linspace(-t_half, t_half, 256)
            delta = leaves[:, :-2, 0] - truth_l[:po]
            exc = phase_excursion(delta, tau, f0, po) / tol
            scores = np.asarray(prn.world_tree.scores)[:n]
            out.levels.append(LevelRecord(
                level=lvl, n_leaves=n,
                min_excursion=float(exc.min()),
                n_within_1=int((exc <= 1.0).sum()),
                n_within_3=int((exc <= 3.0).sum()),
                score_max=float(scores.max()) if len(scores) else 0.0,
                t_half=float(t_half),
            ))
    return out
