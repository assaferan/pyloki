import logging, warnings, time, sys
logging.disable(logging.CRITICAL); warnings.filterwarnings("ignore")
import numpy as np
from pyloki.detection.thresholding import DynamicThresholdScheme
KW = dict(branching_pattern=np.array([2.0]*32), ref_ducy=0.1, nbins=64,
          ntrials=1024, nprobs=10, nthresholds=100, snr_final=8.0)
mode = sys.argv[1]
d = DynamicThresholdScheme(mode=mode, **KW); d.run(thres_neigh=11)   # warm
ts = []
for _ in range(3):
    d = DynamicThresholdScheme(mode=mode, **KW)
    t = time.perf_counter(); d.run(thres_neigh=11); ts.append(time.perf_counter()-t)
print(f"{mode:8s} run() unseeded, 32 stages: {min(ts):.3f} s (best of 3)")
