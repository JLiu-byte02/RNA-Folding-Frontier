"""The frozen H2.4 operating-point order (Amendment 4 §A4.1), tested directly.

These cases are deliberately synthetic. The point of the unit test is the
*rule*, not the biology: whichever operating point the rule would pick must be
the one the report uses, and no third fallback may appear.
"""

from __future__ import annotations

import pytest

from rff.analysis.operating_point import Scan, select_operating_point


def test_exact_recall_match_wins_even_when_pair_budget_also_matches():
    scan = [
        Scan(0.50, 0.400, 40),
        Scan(0.60, 0.600, 30),
        Scan(0.70, 0.800, 20),
    ]
    picked = select_operating_point(scan, dp_recall=0.600, dp_pairs=30)
    assert picked["basis"] == "recall_matched"
    assert picked["recall_matched"]["threshold"] == 0.60
    assert picked["recall_matched"]["matches"] == "exact"


def test_nearest_point_above_recall_is_used_when_no_exact_match():
    """A4.1: greedy must not be evaluated at the weaker operating point."""
    scan = [
        Scan(0.50, 0.500, 40),
        Scan(0.60, 0.700, 30),
        Scan(0.70, 0.900, 10),
    ]
    picked = select_operating_point(scan, dp_recall=0.800, dp_pairs=12)
    assert picked["basis"] == "recall_matched"
    assert picked["recall_matched"]["threshold"] == 0.70
    assert picked["recall_matched"]["matches"] == "nearest_above"
    assert picked["recall_matched"]["gap"] == pytest.approx(0.1)


def test_pair_budget_is_the_only_fallback_when_no_threshold_reaches_recall():
    scan = [
        Scan(0.50, 0.300, 25),
        Scan(0.60, 0.400, 20),
        Scan(0.70, 0.500, 15),
    ]
    picked = select_operating_point(scan, dp_recall=0.900, dp_pairs=21)
    assert picked["basis"] == "pair_budget"
    assert picked["pair_budget"]["threshold"] == 0.60
    assert picked["pair_budget"]["gap_pairs"] == 1


def test_no_match_returns_none_rather_than_inventing_a_fallback():
    scan = [
        Scan(0.50, 0.300, 10),
        Scan(0.60, 0.400, 5),
    ]
    picked = select_operating_point(scan, dp_recall=0.900, dp_pairs=60)
    assert picked["basis"] is None
    assert picked["pair_budget"]["feasible"] is False


def test_empty_scan_is_an_error_not_a_silent_verdict():
    with pytest.raises(ValueError):
        select_operating_point([], dp_recall=0.5, dp_pairs=10)
