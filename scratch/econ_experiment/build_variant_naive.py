"""Rebuild variant_naive/ as a copy of the CURRENT src/pyloki with the three
economize_taylor_params calls stripped.

Must be re-run after any change to src/pyloki, or the variant silently runs stale
code -- that is exactly how the first naive run picked up the pre-fix prune.py and
segfaulted while econ succeeded.
"""

import pathlib
import re
import shutil

REPO = pathlib.Path("/Users/assaferan/Documents/GitHub/pyloki")
DEST = REPO / "scratch/econ_experiment/variant_naive"

shutil.rmtree(DEST, ignore_errors=True)
DEST.mkdir(parents=True)
shutil.copytree(REPO / "src/pyloki", DEST / "pyloki")
for pc in (DEST / "pyloki").rglob("__pycache__"):
    shutil.rmtree(pc, ignore_errors=True)

p = DEST / "pyloki/core/taylor.py"
s = p.read_text()
before = s.count("economize_taylor_params")
# resolve only ever reads the low 3 coefficients, so removing the economization
# call IS naive truncation -- no other edit is needed.
s = s.replace("""    dvec_t_add = transforms.economize_taylor_params(
        dvec_t_add,
        half_width_add,
        n_keep=3,
    )
""", "")
s = s.replace("""        dvec_t_seg = transforms.economize_taylor_params(
            dvec_t_seg,
            half_width_seg,
            n_keep=3,
        )
""", "")
s = re.sub(r" *# NB: with a single order dropped.*\n *# identical to naive truncation.*\n", "", s)
if "economize_taylor_params" in s:
    msg = f"strip failed: {s.count('economize_taylor_params')} of {before} calls remain"
    raise SystemExit(msg)
p.write_text(s)
print(f"variant_naive rebuilt from current src ({before} economize calls removed)")
