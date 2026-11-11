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


def crosses(a: tuple[int, int], b: tuple[int, int]) -> bool:
    """True when two pairs interleave: a.i < b.i < a.j < b.j, or the mirror."""
    i, j = a
    k, l = b
    return (i < k < j < l) or (k < i < l < j)


def pseudoknot_pairs(pairs: Iterable[tuple[int, int]]) -> set[tuple[int, int]]:
    """Pairs that cross at least one other pair in the same structure."""
    normalised = sorted({(min(i, j), max(i, j)) for i, j in pairs if i != j})
    out = set()
    for idx, a in enumerate(normalised):
        for b in normalised[idx + 1:]:
            if b[0] >= a[1]:
                break  # sorted by first index: no later pair can start before a ends
            if crosses(a, b):
                out.add(a)
                out.add(b)
    return out


def max_nested_subset(pairs: Iterable[tuple[int, int]]) -> set[tuple[int, int]]:
    """Largest non-crossing subset of ``pairs``, ties broken deterministically.

    Amendment 6 makes this the primary scoring transform: every decoder in our
    comparison emits nested structures, so a truth set containing crossing
    (pseudoknotted) pairs is reduced before the primary F1. Otherwise a nested
    predictor is penalised for pairs it is structurally forbidden to emit.

    Method: pairs that cross nothing are always compatible and are kept. The
    remaining crossing pairs are split into connected components of the crossing
    graph; each component is solved exactly by enumeration when it has at most
    ``EXACT_LIMIT`` pairs (real pseudoknot components are tiny), otherwise by a
    documented greedy fallback so the function never explodes on a pathological
    label. Tie-break: among equal-sized optima, prefer the lexicographically
    smallest sorted pair list, so the result is deterministic.
    """
    normalised = sorted({(min(i, j), max(i, j)) for i, j in pairs if i != j})
    if not normalised:
        return set()
    crossing = pseudoknot_pairs(normalised)
    if not crossing:
        return set(normalised)

    free = [p for p in normalised if p not in crossing]
    node_list = sorted(crossing)
    index = {p: k for k, p in enumerate(node_list)}
    adjacency = {p: set() for p in node_list}
    for a_i, a in enumerate(node_list):
        for b in node_list[a_i + 1:]:
            if b[0] >= a[1]:
                break
            if crosses(a, b):
                adjacency[a].add(b)
                adjacency[b].add(a)

    seen: set[tuple[int, int]] = set()
    chosen = set(free)
    exact_limit = 16
    for start in node_list:
        if start in seen:
            continue
        stack = [start]
        component: list[tuple[int, int]] = []
        seen.add(start)
        while stack:
            node = stack.pop()
            component.append(node)
            for nxt in adjacency[node]:
                if nxt not in seen:
                    seen.add(nxt)
                    stack.append(nxt)
        component.sort()

        if len(component) <= exact_limit:
            best_choice: list[tuple[int, int]] = []
            for mask in range(1 << len(component)):
                picked = [component[k] for k in range(len(component)) if mask >> k & 1]
                if any(
                    crosses(picked[a], picked[b])
                    for a in range(len(picked))
                    for b in range(a + 1, len(picked))
                ):
                    continue
                if (
                    len(picked) > len(best_choice)
                    or (
                        len(picked) == len(best_choice)
                        and picked and best_choice
                        and picked < best_choice
                    )
                ):
                    best_choice = picked
            chosen.update(best_choice)
        else:
            # Pathological label: greedy by shortest span, then lexicographic.
            for p in sorted(component, key=lambda q: (q[1] - q[0], q)):
                if not any(crosses(p, kept) for kept in chosen if kept in crossing):
                    chosen.add(p)
    return chosen
