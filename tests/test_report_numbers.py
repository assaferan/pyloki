"""Bind `04_upstream_report.md` to computed values, so prose cannot drift from data.

Two failures motivated this, both in one afternoon and both of the same kind: a correct
table with an incorrect sentence written over it. A number check alone would have caught
neither, so this does three things --

1. the cheap groups are recomputed live and must still match the committed JSON;
2. every headline figure quoted in the report must equal the JSON value it came from,
   at the precision the report prints it;
3. the withdrawn `branch_max` verdicts stay withdrawn -- the report must not re-assert
   which configurations build.

Section 3 used to do the opposite. It checked the `branch_max` verdicts as booleans
against the report's wording, on the reasoning that "which configurations build" was the
claim got wrong in prose while the table beside it was right. The table was not right: it
came from a proxy that does not measure the guarded quantity (D120), so section 3 spent
five assertions pinning a falsified number to the prose that quoted it, and the suite that
gave the table its authority became the thing that would have blocked its removal.

The lesson, and the reason this docstring is longer than the tests below: asserting that a
figure still reproduces is not asserting that it measures the claim, and only the first is
cheap to automate. A test suite cannot tell you that you computed the wrong quantity. It
can only stop you from quietly re-adopting one you already withdrew, which is all
section 3 now attempts.
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
    """Walk DATA by a `/`-separated path.

    `/` rather than `.` because several keys contain a dot (`ducy_0.10`), which a
    dot-separated walk silently splits into two missing keys.
    """
    node = DATA
    for part in path.split("/"):
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
    """Pins a SUPERSEDED artifact (D120), not a result.

    `branch_max_counts` does not measure the guarded quantity; its output survives in the
    JSON only as an audit trail. This test keeps that artifact stable so the withdrawal
    stays legible -- it asserts nothing about `branch_max` in the shipped code.
    """
    assert RN.branch_max_counts() == DATA["branch_max"]


# --- 2. the report quotes the computed values -----------------------------------

QUOTED = [
    ("nearest_template_taylor/aggressive/median", "{:.3f}"),
    ("nearest_template_taylor/quadrature/closer", "{:d}"),
    ("nearest_template_taylor/quadrature/median_proven_gain", "{:.2f}"),
    ("nearest_template_taylor/conservative/closer", "{:d}"),
    ("nearest_template_taylor/conservative/unresolved", "{:d}"),
    ("nearest_template_chebyshev/aggressive/median", "{:.3f}"),
    ("nearest_template_chebyshev/quadrature/closer", "{:d}"),
    ("nearest_template_chebyshev/quadrature/median_proven_gain", "{:.2f}"),
    ("corner_forms/taylor/4", "{:.1f}"),
    ("corner_forms/chebyshev/4", "{:.1f}"),
]


@pytest.mark.parametrize(("path", "fmt"), QUOTED)
def test_report_quotes_computed_value(path, fmt):
    text = fmt.format(_dig(path))
    assert text in LIVE, f"{path} = {text} is not quoted in the live report"


@pytest.mark.parametrize("basis", ["taylor", "chebyshev"])
@pytest.mark.parametrize("strategy", ["aggressive", "quadrature"])
def test_report_quotes_amplitude_loss(basis, strategy):
    """The report prints these as percentages to 2 dp."""
    value = _dig(f"amplitude_{basis}/{strategy}/ducy_0.10/boxcar")
    assert f"{100 * value:.2f}" in LIVE


@pytest.mark.parametrize("basis", ["taylor", "chebyshev"])
def test_report_quotes_the_paired_advantage(basis):
    """And quotes the PAIRED statistic, not the ratio of the two medians.

    They differ by about a factor of two (+0.47% against +1.00% at 10% duty in the
    Taylor basis), and an earlier draft printed the paired figure beside the two
    medians as though dividing them would reproduce it.
    """
    adv = _dig(f"amplitude_{basis}/quadrature_advantage/ducy_0.10/boxcar")
    assert f"{100 * adv['per_cell_median']:.2f}" in LIVE
    assert adv["per_cell_median"] < adv["ratio_of_medians"]


def test_taylor_advantage_statistics_really_do_differ():
    """Guards the caveat above: if they ever coincide, the warning is noise."""
    adv = _dig("amplitude_taylor/quadrature_advantage/ducy_0.10/boxcar")
    assert adv["ratio_of_medians"] > 1.5 * adv["per_cell_median"]


# --- 3. the withdrawn branch_max verdicts stay withdrawn -------------------------

# Phrasings that state a per-configuration verdict. The guard's own mechanics ("16
# builds and 17 raises") are NOT verdicts and must stay sayable -- that sentence is the
# only part of the section that survived D120.
VERDICT_MARKERS = ("— raises", "— at the limit", "configuration that builds",
                   "configurations that build", "unreachable at the default")


def _verdict_lines(text: str) -> list[str]:
    """Paragraphs stating a build verdict OUTSIDE a withdrawal context.

    A withdrawn claim has to be quotable in order to be withdrawn -- the README's
    withdrawn-table row names it verbatim, and the report's withdrawal paragraph restates
    it before refuting it -- so a passage is only a violation when it asserts the verdict
    without marking it as retracted.

    Paragraph granularity, not line: prose here is hard-wrapped, so a claim and the
    "withdrawn" that governs it routinely sit on different lines.
    """
    blocks = re.split(r"\n\s*\n", text)
    return [b.strip() for b in blocks
            if any(m in b for m in VERDICT_MARKERS)
            and not re.search(r"withdraw|WITHDRAWN|D120|D121|NOT established"
                              r"|not supported|Do not quote", b)]


@pytest.mark.parametrize("doc", ["04_upstream_report.md", "README.md"])
def test_no_live_branch_max_verdict(doc):
    """No document may state which configurations build (D120)."""
    offenders = _verdict_lines((DOCS / doc).read_text())
    assert offenders == [], (
        f"{doc} states a withdrawn branch_max verdict:\n" + "\n".join(offenders))


def test_report_carries_the_withdrawal():
    """So the section cannot be quietly dropped instead of marked withdrawn."""
    assert "WITHDRAWN IN FULL" in LIVE
    assert "D120" in LIVE
    # The refutation is the load-bearing part: a successor who deletes it loses the
    # reason and is free to rebuild the same proxy.
    assert "40 full 63-level prunes" in LIVE


# --- 4. the measured replacement (D121) -----------------------------------------

MEASURED = DATA["branch_max_measured"]


def _sections(text: str) -> list[str]:
    """Split markdown on headings, so a table stays with the prose that frames it."""
    parts, cur = [], []
    for line in text.splitlines():
        if line.startswith("#") and cur:
            parts.append("\n".join(cur))
            cur = []
        cur.append(line)
    if cur:
        parts.append("\n".join(cur))
    return parts


def test_withdrawn_counts_never_appear_as_live_figures():
    """The proxy's numbers may be named as withdrawn, never quoted as the answer.

    Kept as a distinct check from section 3 because the two fail differently: section 3
    catches a reinstated *verdict*, this catches a reinstated *number* -- e.g. a table
    that quietly restores 28 for Taylor+quadrature without the word "raises" anywhere.
    """
    withdrawn = {(b, s): e["max_per_axis"]
                 for b, row in DATA["branch_max"].items() for s, e in row.items()}
    for (basis, strategy), old in withdrawn.items():
        new = MEASURED[basis][strategy]
        if old == new:
            continue        # 7 is both the old Taylor/aggressive and the new conservative
        for doc in ("04_upstream_report.md", "README.md"):
            for section in _sections((DOCS / doc).read_text()):
                # Scope by SECTION, not paragraph. Bare numbers recur everywhere --
                # "| 20 |" is a prune level in the survival table -- so an unscoped
                # search reports those. Paragraph scoping is too tight in the other
                # direction: a table one blank line below its own heading would escape.
                if "branch_max" not in section and "num_points" not in section:
                    continue
                if not re.search(rf"\|\s*\**{old}\**\s*(—|\|)", section):
                    continue
                assert re.search(r"withdraw|WITHDRAWN|D120", section), (
                    f"{doc} quotes the withdrawn count {old} for {basis}/{strategy} "
                    f"in a live table; the measured value is {new} (D121)")


@pytest.mark.parametrize("basis", ["taylor", "chebyshev"])
def test_docs_quote_the_measured_row(basis):
    """Both documents must carry the measured row exactly as recorded."""
    row = [MEASURED[basis][s] for s in ("aggressive", "quadrature", "conservative")]
    for doc in ("04_upstream_report.md", "README.md"):
        text = (DOCS / doc).read_text()
        pattern = r"\|\s*" + basis.capitalize() + r"\s*\|" + "".join(
            rf"\s*\**{v}\**\s*\|" for v in row)
        assert re.search(pattern, text), (
            f"{doc} does not carry the measured {basis} row {row} (D121)")


def test_measured_maxima_all_fit_the_shipped_default():
    """The claim the docs make, as a boolean over the recorded numbers."""
    worst = max(v for b in ("taylor", "chebyshev")
                for v in MEASURED[b].values())
    assert worst <= MEASURED["shipped_branch_max"]
    assert MEASURED["all_build"] is True


def test_measured_maxima_are_not_buffer_limited():
    """A maximum that moves with max_sugg would make the table conditional."""
    for cell, row in MEASURED["buffer_dependence"].items():
        if cell.startswith("_"):
            continue
        assert len(set(row)) == 1, (
            f"{cell} moves with the buffer: {row}; the table is then conditional "
            "on max_sugg and must say so")
