"""Exact dynamic programming over non-crossing base-pair sets.

The structured decoder family in E1' needs the same object three times: a
maximum-weight non-crossing pairing (MFE decode), its partition function (for
thermodynamic quantities and likelihood training), and base-pair marginals
(for MEA decode and calibrated confidence). This module provides the max-plus
version first because it can be verified exhaustively on short sequences.

A structure is a set of pairs (i, j) with j - i > ``min_loop`` and no two pairs
crossing. The score of a structure is the sum of the pair scores, which is the
Nussinov parameterisation; decoders that add stacking or loop terms must express
them as pair-local scores or extend the recursion explicitly.
"""

from __future__ import annotations

import numpy as np


def nussinov(
    scores: np.ndarray, min_loop: int = 3, return_matrix: bool = False
) -> tuple[list[tuple[int, int]], float]:
    """Maximum-weight non-crossing pairing.

    ``scores[i, j]`` is the additive weight of pairing i with j (higher is
    better). Unpaired positions contribute nothing. Ties are broken towards the
    shorter partner, which keeps the output deterministic.
    """
    s = np.asarray(scores, dtype=np.float64)
    if s.ndim != 2 or s.shape[0] != s.shape[1]:
        raise ValueError("scores must be square")
    n = s.shape[0]
    m = np.zeros((n + 1, n + 1), dtype=np.float64)
    if n == 0:
        return [], 0.0
    arg = np.full((n, n), -1, dtype=np.int32)

    for d in range(1, n):
        js = np.arange(d, n)
        is_ = js - d
        for i, j in zip(is_, js):
            best = m[i, j - 1]
            arg[i, j] = -1
            lo = i + min_loop + 1
            if lo <= j:
                ts = np.arange(lo, j + 1)
                vals = s[i, ts] + m[i + 1, ts - 1] + m[ts + 1, j]
                k = int(np.argmax(vals))
                if vals[k] > best:
                    best = float(vals[k])
                    arg[i, j] = int(ts[k])
            m[i, j] = best

    pairs: list[tuple[int, int]] = []
    stack = [(0, n - 1)]
    while stack:
        i, j = stack.pop()
        if i >= j:
            continue
        t = arg[i, j]
        if t < 0:
            stack.append((i, j - 1))
        else:
            pairs.append((int(i), int(t)))
            stack.append((i + 1, int(t) - 1))
            stack.append((int(t) + 1, j))
    pairs.sort()
    if return_matrix:
        return pairs, float(m[0, n - 1]), m
    return pairs, float(m[0, n - 1])


def pair_map_to_pairs(
    prob: np.ndarray, min_loop: int = 3, threshold: float = 0.5
) -> list[tuple[int, int]]:
    """Greedy non-crossing decode of a pair-probability map.

    Pairs are accepted in decreasing probability order while they stay
    compatible with everything already accepted. This is the standard low-cost
    decoder; ``nussinov`` on log-odds is the structured alternative.
    """
    n = prob.shape[0]
    order = [
        (float(prob[i, j]), i, j)
        for i in range(n)
        for j in range(i + min_loop + 1, n)
        if prob[i, j] >= threshold
    ]
    order.sort(reverse=True)
    used = np.zeros(n, dtype=bool)
    pairs: list[tuple[int, int]] = []
    for _p, i, j in order:
        if used[i] or used[j]:
            continue
        if any((i < a < j < b) or (a < i < b < j) for a, b in pairs):
            continue
        used[i] = used[j] = True
        pairs.append((i, j))
    return sorted(pairs)
