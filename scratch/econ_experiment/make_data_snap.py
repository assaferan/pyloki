"""Generate the poly_order=4 dataset: accel + jerk + SNAP.

The old data.npz has no snap at all and tobs = 2.1 s, so its resolve error is
~1e-7 cycles -- it can never show an SNR difference. Snap here is solved so the
NAIVE resolve error is exactly one phase bin at the real segment scale
(half_width_add = 0.5 * tseg_ffa), which is the precondition for any SNR gap.
"""

import numpy as np

import econ_metrics as em
from pyloki.simulation.pulse import PulseSignalConfig

# Config validated by size_experiment.py: all five preconditions pass.
# nseg is deliberately small -- the neglected 2nd-order Doppler term grows as
# nseg**7/(freq*tseg) relative to the 1st-order truncation error under test, so
# many short segments would swamp the effect. ACCEL_WIDEN decouples the accel grid
# width from snap (from_upper otherwise derives one from the other).
PERIOD, DT, NSAMPS_POW, NSEG, NBINS = 0.007, 64e-6, 20, 4, 64
BRUTE_DIV, ACCEL_WIDEN = 128, 20.0
ACCEL, JERK, SNR_IN = 500.0, 6.0, 12.0
TARGET_BINS = 1.0

nsamps = 2**NSAMPS_POW
freq = 1.0 / PERIOD
tseg_ffa = (nsamps // NSEG) * DT
snap = em.snap_for_phase_budget(JERK, ACCEL, freq, tseg_ffa, TARGET_BINS / NBINS)

print(f"tobs={nsamps * DT:.4f}s  freq={freq:.6f}Hz  tseg_ffa={tseg_ffa:.4f}s")
print(f"solved snap = {snap:.6g} m/s^4 for {TARGET_BINS} phase bin(s) of naive error")
e_n, e_e = em.phase_error_pair(np.array([snap, JERK, ACCEL, 0.0, 0.0]), freq, tseg_ffa)
print(f"naive={e_n:.6e} cyc  econ={e_e:.6e} cyc  ratio={e_e / e_n:.4f}")

cfg = PulseSignalConfig(
    period=PERIOD, dt=DT, nsamps=nsamps, snr=SNR_IN, ducy=0.1,
    mod_kwargs={"acc": ACCEL, "jerk": JERK, "snap": snap},
)
tim_data = cfg.generate(shape="gaussian")

out = "<repo>/scratch/econ_experiment/data_snap.npz"
np.savez(
    out, ts_e=tim_data.ts_e, ts_v=tim_data.ts_v, dt=tim_data.dt,
    period=PERIOD, accel=ACCEL, jerk=JERK, snap=snap, freq=cfg.freq,
    nsamps=nsamps, nseg=NSEG, nbins=NBINS, snr_in=SNR_IN,
    brute_div=BRUTE_DIV, accel_widen=ACCEL_WIDEN,
)
print(f"saved {out}")
