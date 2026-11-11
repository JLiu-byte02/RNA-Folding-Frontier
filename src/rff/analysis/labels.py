"""Label-agreement analysis for sequences annotated by two sources.

Given two annotations of the *same* sequence, the pair-level agreement between
them bounds how well any predictor can score against either one. This module
computes that agreement and the distribution of disagreement.
"""

from __future__ import annotations

from collections import Counter
from typing import Iterable, Sequence


def pair_prf(
    truth: Iterable[tuple[int, int]], pred: Iterable[tuple[int, int]]
) -> dict[str, float]:
    """Precision/recall/F1 of ``pred`` against ``truth`` at the base-pair level."""
    t = set(truth)
    p = set(pred)
    tp = len(t & p)
    fp = len(p - t)
    fn = len(t - p)
    precision = tp / (tp + fp) if tp + fp else 1.0
    recall = tp / (tp + fn) if tp + fn else 1.0
    f1 = (
        2 * precision * recall / (precision + recall)
        if precision + recall
        else 1.0
    )
    return {
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "precision": precision,
        "recall": recall,
        "f1": f1,
    }


def classify_disagreement(
    a: Sequence[tuple[int, int]],
    b: Sequence[tuple[int, int]],
) -> dict[str, int]:
    """Split disagreements into categories that matter mechanistically.

    - ``shared``: the same pair in both annotations.
    - ``shifted``: an ``a`` pair ``(i, j)`` whose partner ``j`` is paired to a
      neighbouring position in ``b`` (register/offset ambiguity).
    - ``a_only`` / ``b_only``: pairs that exist in one annotation only.
    """
    sa = {(min(i, j), max(i, j)) for i, j in a}
    sb = {(min(i, j), max(i, j)) for i, j in b}
    partner_b: dict[int, int] = {}
    for i, j in sb:
        partner_b[i] = j
        partner_b[j] = i

    counts = Counter(shared=len(sa & sb), shifted=0, a_only=0, b_only=0)
    paired_in_b = set(sa & sb)
    for i, j in sa - sb:
        nb = partner_b.get(i)
        if nb is not None and (min(i, nb), max(i, nb)) not in paired_in_b:
            counts["shifted"] += 1
            paired_in_b.add((min(i, nb), max(i, nb)))
        else:
            counts["a_only"] += 1
    counts["b_only"] = len(sb - sa) - counts["shifted"]
    return dict(counts)


def sequence_composition(seq: str) -> dict[str, float]:
    n = max(1, len(seq))
    counts = Counter(seq)
    return {base: counts.get(base, 0) / n for base in "ACGU"}
