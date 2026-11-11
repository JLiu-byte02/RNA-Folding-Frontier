"""Verify the DP decoders against exhaustive enumeration."""

from __future__ import annotations

import itertools
import sys
from functools import lru_cache
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rff.models.dp import nussinov, pair_map_to_pairs  # noqa: E402


def brute_force(scores: np.ndarray, min_loop: int = 3) -> float:
    """Best total weight over all non-crossing pairings (recursive, small n)."""
    n = scores.shape[0]

    @lru_cache(maxsize=None)
    def best(i: int, j: int) -> float:
        if i >= j:
            return 0.0
        val = best(i, j - 1)
        for t in range(i + min_loop + 1, j + 1):
            val = max(val, scores[i, t] + best(i + 1, t - 1) + best(t + 1, j))
        return val

    return best(0, n - 1)


def score_of(pairs, scores):
    return float(sum(scores[i, j] for i, j in pairs))


def test_matches_brute_force():
    rng = np.random.default_rng(0)
    for n in range(0, 13):
        for _ in range(5):
            scores = rng.normal(size=(n, n)) * 2.0
            scores = (scores + scores.T) / 2
            expected = brute_force(scores)
            pairs, val = nussinov(scores)
            assert abs(val - expected) < 1e-6, (n, val, expected)
            assert abs(score_of(pairs, scores) - val) < 1e-6


def test_pairs_are_nested_and_non_crossing():
    rng = np.random.default_rng(1)
    n = 60
    scores = rng.normal(size=(n, n))
    scores = (scores + scores.T) / 2
    pairs, _val = nussinov(scores, min_loop=3)
    positions = [x for p in pairs for x in p]
    assert len(positions) == len(set(positions))
    for (i, j), (a, b) in itertools.combinations(pairs, 2):
        assert j - i > 3
        assert not (i < a < j < b) and not (a < i < b < j)


def test_positive_scores_are_used():
    n = 20
    scores = np.zeros((n, n))
    scores[0, 7] = 5.0
    pairs, val = nussinov(scores)
    assert (0, 7) in pairs and val == 5.0


def test_greedy_decode_stays_non_crossing():
    rng = np.random.default_rng(2)
    n = 40
    prob = rng.random((n, n)) * 0.6
    pairs = pair_map_to_pairs(prob, threshold=0.5)
    positions = [x for p in pairs for x in p]
    assert len(positions) == len(set(positions))
    for (i, j), (a, b) in itertools.combinations(pairs, 2):
        assert not (i < a < j < b) and not (a < i < b < j)
