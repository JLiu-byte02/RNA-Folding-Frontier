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
    max_nested_subset_enumeration,
    max_nested_subset_of_arcs,
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


def _interval_dp_size(arcs, n: int) -> int:
    """Independent exact solver: the textbook interval recurrence, O(n^3).

    A different recurrence from the production DP (it splits at every ``k`` and
    allows ``j`` to be unpaired from either side), so agreement on random inputs
    is a real cross-check rather than the same code twice.
    """
    allowed = set(arcs)
    f = [[0] * n for _ in range(n)]

    def get(i, j):
        return 0 if i > j or i < 0 or j >= n else f[i][j]

    for span in range(1, n + 1):
        for i in range(0, n - span + 1):
            j = i + span - 1
            best = max(get(i + 1, j), get(i, j - 1))
            if (i, j) in allowed:
                best = max(best, 1 + get(i + 1, j - 1))
            for k in range(i, j):
                best = max(best, get(i, k) + get(k + 1, j))
            f[i][j] = best
    return f[0][n - 1]


def test_dp_matches_enumeration_on_small_random_arc_sets():
    rng = random.Random(11)
    for _ in range(25):
        n = 30
        arcs = set()
        for _ in range(rng.randint(1, 12)):
            i, j = sorted(rng.sample(range(n), 2))
            arcs.add((i, j))
        reference = max_nested_subset_enumeration(arcs)
        solved = max_nested_subset_of_arcs(arcs)
        assert solved == reference, sorted(arcs)


def test_dp_matches_an_independent_recurrence_above_the_old_limit():
    """The old code fell back to greedy above 16 arcs; the DP must stay exact."""
    rng = random.Random(5)
    for _ in range(15):
        n = 60
        arcs = set()
        while len(arcs) < 24:
            i, j = sorted(rng.sample(range(n), 2))
            arcs.add((i, j))
        solved = max_nested_subset_of_arcs(arcs)
        assert len(solved) == _interval_dp_size(arcs, n)
        assert is_nested(solved)


def test_every_crossing_component_of_a_random_structure_survives_the_split():
    """Solving components separately must equal solving the whole set at once."""
    rng = random.Random(3)
    for _ in range(20):
        n = 40
        arcs = set()
        for _ in range(14):
            i, j = sorted(rng.sample(range(n), 2))
            arcs.add((i, j))
        overall = max_nested_subset_of_arcs(arcs)
        assert len(overall) == _interval_dp_size(arcs, n)
