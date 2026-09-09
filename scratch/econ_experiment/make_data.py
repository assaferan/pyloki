import numpy as np
from pyloki.simulation.pulse import PulseSignalConfig

pulsar_period = 0.007
dt = 64e-6
accel = 500.0
jerk = 6.0
nsamps = 2**15

mod_kwargs = {"acc": accel, "jerk": jerk}
cfg = PulseSignalConfig(
    period=pulsar_period,
    dt=dt,
    nsamps=nsamps,
    snr=12.0,
    ducy=0.1,
    mod_kwargs=mod_kwargs,
)
print(f"tobs = {cfg.tobs:.3f} s, freq = {cfg.freq:.6f} Hz")
tim_data = cfg.generate(shape="gaussian")

np.savez(
    "/Users/assaferan/Documents/GitHub/pyloki/scratch/econ_experiment/data.npz",
    ts_e=tim_data.ts_e,
    ts_v=tim_data.ts_v,
    dt=tim_data.dt,
    period=pulsar_period,
    accel=accel,
    jerk=jerk,
    freq=cfg.freq,
    nsamps=nsamps,
)
print("Saved simulated dataset.")
