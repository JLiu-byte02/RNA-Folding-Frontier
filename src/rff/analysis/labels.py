"""Label-agreement analysis for sequences annotated by two sources.

Given two annotations of the *same* sequence, the pair-level agreement between
them bounds how well any predictor can score against either one. This module
computes that agreement and the distribution of disagreement.
"""

from __future__ import annotations

from collections import Counter, defaultdict
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


def conflicts(a: tuple[int, int], b: tuple[int, int]) -> bool:
    """True when both pairs cannot appear in one structure.

    Two reasons: they interleave (crossing, i.e. a pseudoknot), or they share an
    endpoint (a base with two partners, which no secondary structure allows).
    ``crosses`` alone misses the second case, and a solver that only checks
    crossing happily returns (17, 27) together with (18, 27).
    """
    if a[0] in b or a[1] in b:  # any shared endpoint, including a.j == b.i
        return True
    return crosses(a, b)


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


def conflict_components(
    pairs: Iterable[tuple[int, int]],
) -> tuple[set[tuple[int, int]], list[list[tuple[int, int]]]]:
    """Split ``pairs`` into conflict-free pairs and conflicting components.

    A pair that conflicts with nothing can never be dropped, so it is always
    kept. The rest form a conflict graph whose connected components are
    independent: no arc of one component conflicts with an arc of another,
    because any such conflict would be an edge. Solving each component
    separately is therefore exact, and it keeps the quadratic DP confined to the
    small span of a pseudoknot instead of the whole molecule.

    Conflicts are crossing *and* shared endpoints (see :func:`conflicts`).
    """
    normalised = sorted({(min(i, j), max(i, j)) for i, j in pairs if i != j})
    in_conflict: set[tuple[int, int]] = set()
    for a_i, a in enumerate(normalised):
        for b in normalised[a_i + 1:]:
            if b[0] > a[1]:
                break
            if conflicts(a, b):
                in_conflict.add(a)
                in_conflict.add(b)
    free = {p for p in normalised if p not in in_conflict}
    if not in_conflict:
        return free, []

    node_list = sorted(in_conflict)
    adjacency: dict[tuple[int, int], set[tuple[int, int]]] = {p: set() for p in node_list}
    for a_i, a in enumerate(node_list):
        for b in node_list[a_i + 1:]:
            if b[0] > a[1]:
                break
            if conflicts(a, b):
                adjacency[a].add(b)
                adjacency[b].add(a)

    seen: set[tuple[int, int]] = set()
    components: list[list[tuple[int, int]]] = []
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
        components.append(sorted(component))
    return free, components


def max_nested_subset_of_arcs(arcs: Iterable[tuple[int, int]]) -> set[tuple[int, int]]:
    """Maximum non-crossing subset of a fixed arc set — polynomial, exact.

    PI Round 4g §2: the previous implementation enumerated all subsets of a
    crossing component and fell back to a greedy rule above 16 pairs. The
    problem is Nussinov DP restricted to the given arcs: either ``j`` is
    unpaired, or it closes some allowed arc ``(k, j)``::

        f[i][j] = max(f[i][j-1],
                      max_{(k,j) in A, k >= i} f[i][k-1] + 1 + f[k+1][j-1])

    which is O(n^2 + n|A|) cells and exact for any component size. Tie-break is
    unchanged: among maximum-size solutions, take the lexicographically smallest
    sorted arc list. Reconstruction walks the arcs in ascending order and keeps
    an arc when it lies inside a still-free interval and the interval's optimum
    is preserved — the standard greedy proof that this yields the
    lexicographic minimum.
    """
    started = sorted({(min(i, j), max(i, j)) for i, j in arcs if i != j})
    if not started:
        return set()
    lo = min(i for i, _ in started)
    hi = max(j for _, j in started)
    m = hi - lo + 1
    local = [(i - lo, j - lo) for i, j in started]

    by_end: dict[int, list[int]] = defaultdict(list)
    for i, j in local:
        by_end[j].append(i)
    for ends in by_end.values():
        ends.sort()

    f = [[0] * m for _ in range(m)]

    def best(i: int, j: int) -> int:
        return 0 if i > j else f[i][j]

    for j in range(m):
        ends = by_end.get(j, ())
        for i in range(j, -1, -1):
            value = best(i, j - 1)
            for k in ends:
                if k < i:
                    continue
                candidate = best(i, k - 1) + 1 + best(k + 1, j - 1)
                if candidate > value:
                    value = candidate
            f[i][j] = value

    optimum = f[0][m - 1]
    chosen: list[tuple[int, int]] = []
    intervals: list[tuple[int, int]] = [(0, m - 1)]
    for i, j in local:
        for slot, (left, right) in enumerate(intervals):
            if not (left <= i and j <= right):
                continue
            with_arc = (
                best(left, i - 1) + 1 + best(i + 1, j - 1) + best(j + 1, right)
            )
            if with_arc == f[left][right]:
                chosen.append((i, j))
                replacement = []
                if left <= i - 1:
                    replacement.append((left, i - 1))
                if i + 1 <= j - 1:
                    replacement.append((i + 1, j - 1))
                if j + 1 <= right:
                    replacement.append((j + 1, right))
                intervals[slot:slot + 1] = replacement
            break

    if len(chosen) != optimum:  # pragma: no cover - would be a DP bug
        raise AssertionError(
            f"nested-subset reconstruction kept {len(chosen)} of {optimum} arcs"
        )
    return {(i + lo, j + lo) for i, j in chosen}


def max_nested_subset_enumeration(
    arcs: Iterable[tuple[int, int]], limit: int = 16
) -> set[tuple[int, int]]:
    """Reference solver: exhaustive search, used only to test the DP.

    Exported because the polynomial solver's correctness is checked against it on
    every crossing component of the real label sets (Amendment 5 §A5.4 style:
    a fast path is only trusted once it agrees with a slow, obviously-correct
    one). Refuses sets larger than ``limit`` rather than hanging.
    """
    items = sorted({(min(i, j), max(i, j)) for i, j in arcs if i != j})
    if len(items) > limit:
        raise ValueError(f"{len(items)} arcs exceeds the enumeration limit {limit}")
    k = len(items)
    # Bitmask form keeps the exhaustive search usable on the thousands of
    # 10-16 arc components in the real label sets: validity is one AND per set
    # bit instead of a quadratic pair loop.
    blocking = [0] * k
    for a in range(k):
        for b in range(k):
            if a != b and conflicts(items[a], items[b]):
                blocking[a] |= 1 << b

    best: list[tuple[int, int]] = []
    for mask in range(1 << k):
        size = mask.bit_count()
        if size < len(best):
            continue
        rest = mask
        valid = True
        while rest:
            low = rest & -rest
            if mask & blocking[low.bit_length() - 1]:
                valid = False
                break
            rest ^= low
        if not valid:
            continue
        if size > len(best) or (size == len(best) and best):
            picked = [items[i] for i in range(k) if mask >> i & 1]
            if size > len(best) or picked < best:
                best = picked
    return set(best)


def max_nested_subset(pairs: Iterable[tuple[int, int]]) -> set[tuple[int, int]]:
    """Largest non-crossing subset of ``pairs``, ties broken deterministically.

    Amendment 6 makes this the primary scoring transform: every decoder in our
    comparison emits nested structures, so a truth set containing crossing
    (pseudoknotted) pairs is reduced before the primary F1. Otherwise a nested
    predictor is penalised for pairs it is structurally forbidden to emit.

    Pairs that cross nothing are kept; each crossing component is solved exactly
    by :func:`max_nested_subset_of_arcs`. The tie-break is lexicographic on the
    sorted pair list, and solving components independently reproduces it: the
    first position where two candidates differ lies in whichever component
    contributes the smaller arc, and each component returns its own minimum.
    """
    free, components = conflict_components(pairs)
    chosen = set(free)
    for component in components:
        chosen |= max_nested_subset_of_arcs(component)
    return chosen
