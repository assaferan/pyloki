"""What a phase error costs in S/N, using the shipped scoring rather than a metric.

D44 (0.34%) was withdrawn twice over: it rested on a nearest-template number that was a
cap artefact, and it converted that number through `core/metric.py`, whose harmonic
weighting was an open question (O4) and whose mismatch is a mean-square quantity while
the template was selected by a sup-norm criterion. This module avoids the metric
entirely.

The phase residual `dPhi(t)` between a signal and a template does not attenuate the
pulse, it *smears* it: the folded profile is the true profile averaged over the phase
offsets visited during the observation,

    folded(phi) = <p(phi - dPhi(t))>_t ,

which in Fourier is exactly `P(k) * S(k)` with `S(k) = <exp(-2 pi i k dPhi(t))>_t`, the
characteristic function of the phase-residual distribution. No approximation and no
second-order truncation. Both profiles are then scored with `detection.scoring`'s own
boxcar matched filter over `generate_box_width_trials`, and since folding averages the
same data either way the noise level is identical, so the S/N ratio is the ratio of the
filter outputs.

One caveat that must travel with the numbers: the template is chosen by
`nearest_template.py`, which minimises the sup-norm of `dPhi`, not the S/N loss. The
best-scoring template need not be the sup-norm-nearest one, so evaluating the
sup-norm-nearest leaf gives an **upper bound** on the loss. Passing several near-optimal
leaves and taking the best tightens that bound, which is what `loss_for_leaves` does.
"""

from __future__ import annotations

import numpy as np

from pyloki.detection.scoring import boxcar_snr_1d, generate_box_width_trials
from pyloki.simulation.pulse import generate_folded_profile
from pyloki.utils.misc import C_VAL

from nearest_template import _basis


def smearing_factor(
    delta: np.ndarray, tau: np.ndarray, f0: float, poly_order: int, nbins: int,
    basis_fn=None,
) -> np.ndarray:
    """S(k) = <exp(-2 pi i k dPhi(t))>_t for k = 0..nbins//2.

    `basis_fn` selects the phase basis: the Taylor monomials by default, or
    `nearest_template_cheby.cheby_basis` for Chebyshev coefficients.
    """
    dphi = (f0 / C_VAL) * ((basis_fn or _basis)(tau, poly_order) @ delta)  # cycles
    k = np.arange(nbins // 2 + 1)
    return np.exp(-2j * np.pi * np.outer(k, dphi)).mean(axis=1)


def snr_ratio(
    delta: np.ndarray,
    tau: np.ndarray,
    f0: float,
    poly_order: int,
    nbins: int,
    ducy: float,
    ducy_max: float = 0.2,
    n_phase: int = 16,
    filt: str = "boxcar",
    basis_fn=None,
) -> float:
    """Recovered S/N divided by the S/N at zero phase error.

    Averaged over `n_phase` sub-bin positions of the pulse, because the signal's absolute
    phase is arbitrary and uniformly distributed. Without that average the answer carries
    a discretisation artefact of the same size as the effect being measured: a Gaussian
    smeared to a flat top can match a *boxcar* better than a sharp one does, which showed
    up as a negative loss.

    `filt="boxcar"` uses the pipeline's own filter bank; `filt="matched"` correlates
    against the true profile instead, which is a cleaner and provably monotone reference
    but not what the code scores with. Quoting both is the point -- if they disagree the
    number is an artefact of the filter, not a property of the grid.

    **This returns a profile overlap, and it is NOT a correlation between two templates'
    scores.** Both quantities weight the same smearing factor `S(k)`, but by different
    spectra:

        this function   sum_k |P(k)|^2 S(k) / sum_k |P(k)|^2      (the PULSE's spectrum)
        score corr.     sum_k |H(k)|^2 S(k) / sum_k |H(k)|^2      (the FILTER's)

    They coincide only when the filter is the pulse, and it is not: at `ducy = 0.10` with
    `N_b = 64` the bank selects a boxcar of width 3, whose `<k^2>` is about 4x the
    pulse's, so the smearing is weighted toward higher harmonics and the score
    decorrelates roughly 3x faster than the profile overlap suggests. Feeding a
    *between-template* residual `delta_A - delta_B` to this function and reading the
    result as a score correlation therefore understates the decorrelation substantially.
    Found by the injection-campaign session in its own extension of this module; recorded
    here so the next caller does not repeat it.
    """
    smear = smearing_factor(delta, tau, f0, poly_order, nbins, basis_fn)
    widths = generate_box_width_trials(nbins, ducy_max=ducy_max)
    num = den = 0.0
    for centre in 0.5 + np.arange(n_phase) / (n_phase * nbins):
        prof = generate_folded_profile(nbins=nbins, ducy=ducy, center=centre)
        smeared = np.fft.irfft(np.fft.rfft(prof) * smear, n=nbins)
        if filt == "matched":
            den += float(prof @ prof)
            num += float(prof @ smeared)
        else:
            den += float(boxcar_snr_1d(prof.astype(np.float32), widths, 1.0).max())
            num += float(boxcar_snr_1d(smeared.astype(np.float32), widths, 1.0).max())
    return num / den


def loss_for_leaves(
    leaves: list[tuple[float, np.ndarray]],
    tau: np.ndarray,
    f0: float,
    poly_order: int,
    nbins: int,
    ducy: float,
    ducy_max: float = 0.2,
    filt: str = "boxcar",
    basis_fn=None,
) -> float:
    """Smallest fractional S/N loss over a set of candidate leaves.

    An upper bound on the loss the strategy actually incurs, tighter the more
    near-optimal leaves are supplied.
    """
    return 1.0 - max(
        snr_ratio(d, tau, f0, poly_order, nbins, ducy, ducy_max, filt=filt,
                  basis_fn=basis_fn)
        for _, d in leaves
    )
