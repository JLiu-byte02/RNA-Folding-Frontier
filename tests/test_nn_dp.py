"""Correctness tests for decoder A (nearest-neighbour DP).

The DP is checked against exhaustive enumeration of every secondary structure on
short sequences: for random pair scores and loop penalties the DP optimum must
equal the brute-force optimum, and the traceback must reproduce that optimum.
That establishes the recursion and the traceback independently of any energy
parameters, which is the part a ViennaRNA comparison cannot give us.
"""

from __future__ import annotations

import itertools
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rff.decoders.nn_dp import (  # noqa: E402
    NEG,
    NNDecoder,
    PairScoreModel,
    logit_score,
    mask_allowed,
    pairs_to_dotbracket,
)


def all_structures(n: int, min_loop: int = 3):
    """Every non-crossing pair set with ``j - i > min_loop``, as a tuple of pairs."""

    def rec(i: int, j: int):
        if i > j:
            yield ()
            return
        for rest in rec(i + 1, j):
            yield rest
        for k in range(i + min_loop + 1, j + 1):
            for inner in rec(i + 1, k - 1):
                for rest in rec(k + 1, j):
                    yield tuple(sorted(((i, k),) + inner + rest))

    yield from rec(0, n - 1)


def children_of(pairs, i, j):
    inside = [p for p in pairs if i < p[0] and p[1] < j]
    out = []
    for p in inside:
        if not any(
            q is not p and q[0] < p[0] and p[1] < q[1] for q in inside
        ):
            out.append(p)
    return sorted(out)


def structure_score(pairs, length: int, model: PairScoreModel, S: np.ndarray) -> float:
    """Score a structure exactly as the DP decomposes it."""
    pairs = sorted(pairs)
    total = 0.0

    def score_pair(i, j):
        kids = children_of(pairs, i, j)
        cands = []
        if not kids:
            cands.append(model.hairpin_penalty(j - i - 1))
        if len(kids) == 1:
            k, l = kids[0]
            left, right = k - i - 1, j - l - 1
            cands.append(
                model.interior_penalty(left + right, abs(left - right))
                + score_pair(k, l)
            )
        if kids:
            unpaired = (j - i - 1) - sum(l - k + 1 for k, l in kids)
            cands.append(
                model.ml_closing
                + unpaired * model.ml_unpaired
                + sum(model.ml_branch + score_pair(k, l) for k, l in kids)
            )
        return float(S[i, j]) + max(cands)

    top = [p for p in pairs if not any(q[0] < p[0] and p[1] < q[1] for q in pairs)]
    total = sum(score_pair(i, j) for i, j in top)
    return total


def random_model(n: int, rng, with_penalties: bool) -> tuple[PairScoreModel, np.ndarray]:
    S = np.full((n, n), NEG)
    for i in range(n):
        for j in range(i + 4, n):
            S[i, j] = rng.normal()
    if not with_penalties:
        return PairScoreModel(pair_score=S), S
    model = PairScoreModel(
        pair_score=S,
        hairpin=np.array([0.0, -0.5, -1.0, -1.5, -2.0, -2.5]),
        interior=np.array([[-0.2 * k - 0.1 * m for m in range(8)] for k in range(31)]),
        ml_closing=-2.0,
        ml_branch=-1.0,
        ml_unpaired=-0.3,
    )
    return model, S


@pytest.mark.parametrize("with_penalties", [False, True])
@pytest.mark.parametrize("seed", [0, 1, 2])
def test_dp_optimal_value_matches_bruteforce(with_penalties, seed):
    rng = np.random.default_rng(seed)
    n = 11
    model, S = random_model(n, rng, with_penalties)
    dec = NNDecoder(model)
    pairs, value = dec.decode("A" * n)

    best = max(structure_score(p, n, model, S) for p in all_structures(n))
    assert value == pytest.approx(best, abs=1e-6)
    assert structure_score(pairs, n, model, S) == pytest.approx(value, abs=1e-6)


@pytest.mark.parametrize("seed", [3, 4])
def test_dp_traceback_is_valid_structure(seed):
    rng = np.random.default_rng(seed)
    n = 14
    model, S = random_model(n, rng, with_penalties=True)
    pairs, _ = NNDecoder(model).decode("A" * n)
    used = [p for pair in pairs for p in pair]
    assert len(used) == len(set(used)), "a position was paired twice"
    for i, j in pairs:
        assert j - i > model.min_loop
    for (a, b), (c, d) in itertools.combinations(pairs, 2):
        crossing = (a < c < b < d) or (c < a < d < b)
        assert not crossing, f"crossing pairs {(a, b)} {(c, d)}"


def test_zero_scores_give_empty_structure():
    n = 12
    S = np.zeros((n, n))
    pairs, value = NNDecoder(PairScoreModel(pair_score=S)).decode("A" * n)
    assert pairs == []
    assert value == 0.0


def test_strong_negative_pair_scores_are_avoided():
    n = 12
    S = np.full((n, n), NEG)
    S[0, 11] = -5.0
    S[1, 10] = -5.0
    pairs, value = NNDecoder(PairScoreModel(pair_score=S)).decode("A" * n)
    # empty structure scores 0, which beats any negative pair
    assert pairs == []
    assert value == 0.0


def test_canonical_mask_restricts_pairs():
    seq = "GGGAAACCC"
    mask = mask_allowed(seq)
    assert mask[0, 8] and mask[1, 7] and mask[2, 6]
    assert not mask[3, 5]
    S = np.where(mask, 1.0, NEG)
    pairs, _ = NNDecoder(PairScoreModel(pair_score=S)).decode(seq)
    for i, j in pairs:
        assert mask[i, j]
        assert seq[i] + seq[j] in {"AU", "UA", "GC", "CG", "GU", "UG"}


def test_logit_score_is_symmetric_and_monotone():
    p = np.array([[0.0, 0.2, 0.9], [0.05, 0.0, 0.5], [0.2, 0.4, 0.0]])
    s = logit_score(p)
    assert np.allclose(s, s.T)
    assert s[0, 2] == pytest.approx(s[2, 0])
    assert s[0, 1] < 0 < s[0, 2]


def test_dotbracket_roundtrip():
    assert pairs_to_dotbracket(6, [(0, 5), (1, 4)]) == "((..))"


def test_min_loop_is_enforced_on_short_sequences():
    # four nucleotides cannot form a pair with a hairpin loop of >= 3
    S = np.full((4, 4), 1.0)
    pairs, _ = NNDecoder(PairScoreModel(pair_score=S)).decode("AAAA")
    assert pairs == []
    # five can: (0, 4) closes a loop of exactly three
    S = np.full((5, 5), 1.0)
    pairs, _ = NNDecoder(PairScoreModel(pair_score=S)).decode("AAAAA")
    assert pairs == [(0, 4)]
