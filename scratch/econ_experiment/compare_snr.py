"""Report the end-to-end SNR comparison once both variants have run."""

import pathlib

import numpy as np

BASE = pathlib.Path(__file__).resolve().parent

rows = {}
for variant in ("naive", "econ"):
    f = BASE / f"snr_{variant}.npz"
    if not f.exists():
        print(f"MISSING: {f.name} -- run_snr.py {variant} has not completed")
        continue
    rows[variant] = np.load(f)

if len(rows) < 2:
    raise SystemExit(1)

print("=" * 66)
print("End-to-end SNR: naive truncation vs Chebyshev economization")
print("=" * 66)
print(f"{'metric':<22} {'naive':>12} {'econ':>12} {'econ - naive':>14}")
print("-" * 66)
for key, label in (
    ("best_score", "best score"),
    ("best_score_ep", "best score_ep"),
    ("closest_score", "closest-to-true score"),
    ("closest_score_ep", "closest score_ep"),
    ("n_cands", "n candidates"),
):
    n, e = float(rows["naive"][key]), float(rows["econ"][key])
    print(f"{label:<22} {n:>12.4f} {e:>12.4f} {e - n:>+14.4f}")

print("-" * 66)
dn = float(rows["econ"]["closest_score"]) - float(rows["naive"]["closest_score"])
if abs(dn) < 1e-9:
    print("VERDICT: identical to floating point -- economization changed nothing.")
else:
    print(f"VERDICT: economization {'HELPS' if dn > 0 else 'HURTS'} the "
          f"closest-to-true candidate by {dn:+.4f} in score.")
