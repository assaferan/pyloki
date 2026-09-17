"""Bind `04_upstream_report.md` to computed values, so prose cannot drift from data.

Two failures motivated this, both in one afternoon and both of the same kind: a correct
table with an incorrect sentence written over it. A number check alone would have caught
neither, so this does three things --

1. the cheap groups are recomputed live and must still match the committed JSON;
2. every headline figure quoted in the report must equal the JSON value it came from,
   at the precision the report prints it;
3. the `branch_max` verdicts are checked as *booleans* against the report's wording,
   because "which configurations build" is the claim that was got wrong in prose while
   the table beside it was right.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import pytest

DOCS = Path(__file__).resolve().parents[1] / "docs" / "metric_gridding"
sys.path.insert(0, str(DOCS))

import report_numbers as RN  # noqa: E402

REPORT = (DOCS / "04_upstream_report.md").read_text()
DATA = json.loads((DOCS / "report_numbers.json").read_text())
# Only the live report matters; the withdrawn draft is kept below a marker.
LIVE = REPORT.split("# ⚠ EVERYTHING BELOW THIS LINE")[0]


def _dig(path: str):
    node = DATA
    for part in path.split("."):
        node = node[part]
    return node


# --- 1. the cheap groups still reproduce ----------------------------------------

def test_corner_forms_reproduce():
    fresh = RN.corner_forms()
    assert {int(k): v for k, v in fresh["taylor"].items()} == {
        int(k): v for k, v in DATA["corner_forms"]["taylor"].items()}
    assert {int(k): v for k, v in fresh["chebyshev"].items()} == {
        int(k): v for k, v in DATA["corner_forms"]["chebyshev"].items()}


def test_corner_forms_match_their_closed_forms():
    """The point of the group: both grids follow an exact formula."""
    fresh = RN.corner_forms()
    for k in range(2, 9):
        assert fresh["taylor"][k] == pytest.approx(2.0 ** (k - 1) - 0.5, rel=1e-9)
        assert fresh["chebyshev"][k] == pytest.approx(k / 2.0, rel=1e-9)


def test_covering_radius_reproduces():
    fresh = RN.covering_radius()
    assert fresh == DATA["covering_radius"]


def test_branch_max_counts_reproduce():
    assert RN.branch_max_counts() == DATA["branch_max"]


# --- 2. the report quotes the computed values -----------------------------------

QUOTED = [
    ("nearest_template_taylor.aggressive.median", "{:.3f}"),
    ("nearest_template_taylor.quadrature.closer", "{:d}"),
    ("nearest_template_taylor.quadrature.median_proven_gain", "{:.2f}"),
    ("nearest_template_taylor.conservative.closer", "{:d}"),
    ("nearest_template_taylor.conservative.unresolved", "{:d}"),
    ("nearest_template_chebyshev.aggressive.median", "{:.3f}"),
    ("nearest_template_chebyshev.quadrature.closer", "{:d}"),
    ("nearest_template_chebyshev.quadrature.median_proven_gain", "{:.2f}"),
    ("corner_forms.taylor.4", "{:.1f}"),
    ("corner_forms.chebyshev.4", "{:.1f}"),
]


@pytest.mark.parametrize(("path", "fmt"), QUOTED)
def test_report_quotes_computed_value(path, fmt):
    text = fmt.format(_dig(path))
    assert text in LIVE, f"{path} = {text} is not quoted in the live report"


@pytest.mark.parametrize("basis", ["taylor", "chebyshev"])
@pytest.mark.parametrize("strategy", ["aggressive", "quadrature"])
def test_report_quotes_amplitude_loss(basis, strategy):
    """The report prints these as percentages to 2 dp."""
    value = _dig(f"amplitude_{basis}.{strategy}.ducy_0.10.boxcar")
    assert f"{100 * value:.2f}" in LIVE


# --- 3. the branch_max claim, as booleans ---------------------------------------

@pytest.mark.parametrize("basis", ["taylor", "chebyshev"])
@pytest.mark.parametrize("strategy", ["aggressive", "quadrature", "conservative"])
def test_report_wording_matches_build_verdict(basis, strategy):
    """`num_points > branch_max` is strict, so equality builds.

    The report must mark a configuration as raising exactly when it does. Getting this
    backwards for Chebyshev+quadrature (16 of 16, which builds) denied the only
    configuration with a proven gain that needs no config change.
    """
    entry = _dig(f"branch_max.{basis}.{strategy}")
    # Guard the match against a digit to its left: a bare substring test makes
    # "7 — raises" match inside "27 — raises".
    marker = rf"(?<!\d){entry['max_per_axis']} — raises"
    found = re.search(marker, LIVE) is not None
    if entry["builds"]:
        assert not found, (
            f"{basis}/{strategy} builds ({entry['max_per_axis']} <= "
            f"{entry['branch_max']}) but the report says it raises")
    else:
        assert found, (
            f"{basis}/{strategy} raises ({entry['max_per_axis']} > "
            f"{entry['branch_max']}) but the report does not say so")


def test_build_verdict_is_consistent_with_the_guard():
    for basis, row in DATA["branch_max"].items():
        for strategy, entry in row.items():
            assert entry["builds"] == (
                entry["max_per_axis"] <= entry["branch_max"]), (basis, strategy)


def test_only_one_nonaggressive_configuration_builds():
    """Stated as a test so it cannot be over-generalised in prose again."""
    builds = {(b, s) for b, row in DATA["branch_max"].items()
              for s, e in row.items() if e["builds"] and s != "aggressive"}
    assert builds == {("chebyshev", "quadrature")}
