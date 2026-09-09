# Chebyshev-economization experiments

Scripts backing the analysis in `HANDOFF.md` (see §5-§7 there for what they
established and §6 for the per-script description). Nothing here is part of the
package; it is kept in git only so the reasoning is reproducible.

Generated artifacts (`*.npz`, `results_*/`, `variant_naive/`) are gitignored.
Regenerate them with the scripts below.

Run from the repo root with the venv active and this directory importable:

    PYTHONPATH=scratch/econ_experiment python scratch/econ_experiment/<script>

Typical order:

1. `parity_check.py` — fast, self-contained proof of the parity limitation.
2. `size_experiment.py --nsamps-pow 20 --nseg 4 --brute-div 128 \
       --accel-widen 20 --target-bins 1.0` — screens a candidate search config
   against the five preconditions (~6 s). Run this before anything expensive.
3. `make_data_snap.py` — writes `data_snap.npz` at the validated config.
4. `build_variant_naive.py` — builds `variant_naive/`. **Re-run after any change
   to `src/pyloki`**, or the naive arm silently runs stale code.
5. `run_snr.py econ` and `run_snr.py naive` — one end-to-end run each. The naive
   arm needs `variant_naive` first on `PYTHONPATH`:
       PYTHONPATH=scratch/econ_experiment/variant_naive:scratch/econ_experiment
6. `replicate.py N` then `analyze_replicates.py` — paired replicates and stats.

`resolve_probe2.py` is the corrected sup-norm probe plus the segment-scale sweep.
`make_data.py` / `run_search.py` are the earlier `poly_order=3` attempt, kept only
for reference: that data has no snap and `tobs = 2.1 s`, so it cannot show an
effect.
