"""The primary-scoring transform frozen in Amendment 6.

A nested decoder cannot emit a crossing pair, so scoring it against a
pseudoknotted truth without reduction charges it for something it is structurally
forbidden to produce. These tests pin the reduction, its determinism, and the
fact that it never increases the pair count.
"""

from __future__ import annotations

import random
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rff.analysis.labels import (  # noqa: E402
    crosses,
    max_nested_subset,
    pseudoknot_pairs,
)


def is_nested(pairs) -> bool:
    pl = sorted(pairs)
    return not any(crosses(pl[a], pl[b]) for a in range(len(pl)) for b in range(a + 1, len(pl)))


def test_nested_input_is_returned_unchanged():
    pairs = {(0, 10), (1, 9), (2, 8), (12, 20)}
    assert max_nested_subset(pairs) == pairs


def test_single_pseudoknot_keeps_the_larger_side():
    # (0, 12) crosses (4, 20); both sides have one pair, so the tie-break decides
    pairs = {(0, 12), (4, 20), (30, 40)}
    kept = max_nested_subset(pairs)
    assert (30, 40) in kept
    assert len(kept) == 2
    assert kept in ({(0, 12), (30, 40)}, {(4, 20), (30, 40)})
    assert is_nested(kept)


def test_maximum_is_two_of_three_chain():
    # (0,10) crosses (5,15); (5,15) crosses (12,20); (0,10) and (12,20) are disjoint
    pairs = {(0, 10), (5, 15), (12, 20)}
    kept = max_nested_subset(pairs)
    assert kept == {(0, 10), (12, 20)}


def test_pseudoknot_pairs_are_exactly_the_crossing_ones():
    nested = {(0, 10), (1, 9), (12, 20)}
    assert pseudoknot_pairs(nested) == set()
    # (4,14) interleaves with (0,10) and (1,9), and (12,20) interleaves with (4,14)
    assert pseudoknot_pairs(nested | {(4, 14)}) == {(0, 10), (1, 9), (4, 14), (12, 20)}


def test_result_is_deterministic_and_nested_on_random_inputs():
    rng = random.Random(7)
    for _ in range(40):
        n = 60
        pairs = set()
        for _ in range(12):
            i, j = sorted(rng.sample(range(n), 2))
            pairs.add((i, j))
        first = max_nested_subset(pairs)
        assert first == max_nested_subset(pairs)
        assert is_nested(first)
        assert len(first) <= len(pairs)


def test_shortest_span_tie_break_is_lexicographic_not_random():
    # four mutually crossing pairs, all size 1: keep the lexicographically first
    pairs = {(0, 9), (1, 8), (2, 7), (3, 6)}
    # these are nested, not crossing, so all four survive
    assert len(max_nested_subset(pairs)) == 4
    crossing = {(0, 4), (1, 5), (2, 6)}
    kept = max_nested_subset(crossing)
    assert len(kept) == 1
    assert kept == {(0, 4)}
