"""Edge cases the base-pair scorer must get right (PI Round 4e).

Each case is checked against the independent implementation in
`scripts/independent_scorer.py`, so a change in one that is not mirrored in the
other fails here rather than silently changing every downstream number.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rff.analysis.scoring import prf  # noqa: E402

_spec = importlib.util.spec_from_file_location(
    "independent_scorer", ROOT / "scripts" / "independent_scorer.py"
)
independent = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(independent)


def both(pred, truth):
    a = prf(*independent.counts(pred, truth))
    b = independent.counts(pred, truth)
    f1_b = independent.f1_from_counts(*b)
    return a, b, f1_b


def test_empty_truth_and_empty_prediction():
    a, b, f1_b = both(set(), set())
    assert b == (0, 0, 0)
    assert a["f1"] == 1.0 and f1_b == 1


def test_empty_truth_with_predictions_is_all_false_positive():
    a, b, f1_b = both({(0, 5)}, set())
    assert b == (0, 1, 0)
    assert a["precision"] == 0.0
    assert a["recall"] == 1.0  # vacuous: no true pairs to miss
    assert a["f1"] == 0.0 and f1_b == 0


def test_empty_prediction_with_truth_is_all_false_negative():
    a, b, f1_b = both(set(), {(0, 5)})
    assert b == (0, 0, 1)
    assert a["precision"] == 1.0
    assert a["recall"] == 0.0
    assert a["f1"] == 0.0 and f1_b == 0


def test_pseudoknotted_truth_pairs_are_scored_as_given():
    """{(0,7), (3,10)} cross. One of them predicted -> TP 1, FN 1."""
    truth = {(0, 7), (3, 10)}
    pred = {(0, 7)}
    a, b, f1_b = both(pred, truth)
    assert b == (1, 0, 1)
    assert a["f1"] == pytest.approx(2 / 3)
    assert f1_b == independent.f1_from_counts(1, 0, 1) == pytest.approx(2 / 3)


def test_non_canonical_truth_pairs_are_not_filtered_by_the_scorer():
    """Scoring never looks at the sequence; filtering is a build-time decision."""
    truth = {(0, 4), (10, 14)}  # e.g. A-A pairs after a bad label source
    pred = {(0, 4)}
    a, b, _ = both(pred, truth)
    assert b == (1, 0, 1)
    assert a["f1"] == pytest.approx(2 / 3)


def test_pm1_slip_is_detectable_and_reported_separately():
    truth = {(10, 20), (30, 40)}
    pred = set(truth)
    assert independent.f1_from_counts(*independent.counts(pred, truth)) == 1
    for delta in (1, -1):
        slipped = independent.shift_pairs(pred, delta)
        tp, fp, fn = independent.counts(slipped, truth)
        assert (tp, fp, fn) == (0, 2, 2)
        assert independent.f1_from_counts(tp, fp, fn) == 0


def test_shift_drops_pairs_that_fall_off_the_sequence():
    assert independent.shift_pairs({(0, 5)}, -1) == set()
