"""Scoring primitives for the P2 diagnostics.

The pre-registered metric is a single F1 at a fixed 0.5 threshold. That number
cannot separate two very different situations: a predictor that got worse, and a
predictor whose logits shifted scale so the fixed cut now behaves like a
different operating point. These helpers compute the threshold-free version
(average precision) and an oracle threshold, plus the stratifications the PI
asked for.

Everything here is plain numpy so it can be tested without a GPU or a model.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class Region:
    """The kept (untouched) interval, half-open, plus its complement."""

    lo: int
    hi: int

    def where(self, i: int, j: int) -> str:
        a = self.lo <= i < self.hi
        b = self.lo <= j < self.hi
        if a and b:
            return "both_kept"
        if a or b:
            return "one_kept"
        return "neither_kept"


def distance_bin(i: int, j: int) -> str:
    d = abs(i - j)
    if d < 32:
        return "d<32"
    if d <= 128:
        return "d32-128"
    return "d>128"


def length_band(length: int) -> str:
    """Same bands as the P2 table, so D1 and P2 can be read side by side."""
    for upper, name in ((80, "L<80"), (150, "80-149"), (300, "150-299"), (10**9, "L>=300")):
        if length < upper:
            return name
    return "L>=300"


def decode_pairs(
    scored: dict[tuple[int, int], float],
    threshold: float,
    n: int,
    min_loop: int = 3,
    partner: bool = True,
    nesting: bool = True,
) -> list[tuple[int, int]]:
    """Greedy decode of a sparse probability map, in three nested modes.

    ``partner=False`` is a bare threshold. ``partner=True, nesting=False`` accepts
    in decreasing probability under only the one-partner constraint. Both true is
    the production decoder, matching ``pair_map_to_pairs`` exactly (same tie-break:
    larger probability, then larger ``i``, then larger ``j``).
    """
    cands = [
        (p, i, j)
        for (i, j), p in scored.items()
        if p >= threshold and j - i > min_loop
    ]
    if not partner:
        return sorted((i, j) for _p, i, j in cands)
    cands.sort(reverse=True)
    used = [False] * n
    out: list[tuple[int, int]] = []
    for _p, i, j in cands:
        if used[i] or used[j]:
            continue
        if nesting and any((i < a < j < b) or (a < i < b < j) for a, b in out):
            continue
        used[i] = used[j] = True
        out.append((i, j))
    return sorted(out)


def normalise_scored(scored: dict) -> dict[tuple[int, int], float]:
    """Accept both ``{(i, j): p}`` and ``{"i,j": p}`` (the dump format)."""
    out: dict[tuple[int, int], float] = {}
    for key, value in scored.items():
        if isinstance(key, str):
            i, j = key.split(",")
            out[(int(i), int(j))] = float(value)
        else:
            i, j = key
            out[(int(i), int(j))] = float(value)
    return out


def increment_sweep(
    rows: list[dict],
    cuts: list[float],
    partner: bool = True,
    nesting: bool = True,
) -> list[list[tuple[int, int]]]:
    """Decode every row at every cut in one pass over its candidates.

    A candidate only affects the cuts at or below its own probability, so one
    descending walk maintains all decoders simultaneously and costs
    ``sum_candidates (#cuts <= p)`` instead of ``len(cuts) * n_candidates``.
    """
    per_row: list[list[list[tuple[int, int]]]] = []
    for row in rows:
        scored = normalise_scored(row["_scored"])
        n = row["length"]
        cands = sorted(
            ((p, i, j) for (i, j), p in scored.items() if j - i > 3),
            reverse=True,
        )
        used = [[False] * n for _ in cuts]
        out: list[list[tuple[int, int]]] = [[] for _ in cuts]
        for p, i, j in cands:
            for c, cut in enumerate(cuts):
                if p < cut:
                    break
                if not partner:
                    out[c].append((i, j))
                    continue
                if used[c][i] or used[c][j]:
                    continue
                if nesting and any(
                    (i < a < j < b) or (a < i < b < j) for a, b in out[c]
                ):
                    continue
                used[c][i] = used[c][j] = True
                out[c].append((i, j))
        per_row.append([sorted(o) for o in out])
    return per_row


def prf(tp: int, fp: int, fn: int) -> dict:
    precision = tp / (tp + fp) if tp + fp else 1.0
    recall = tp / (tp + fn) if tp + fn else 1.0
    if tp + fp + fn == 0:
        # Vacuous case: nothing true and nothing predicted.
        f1 = 1.0
    elif precision + recall == 0:
        f1 = 0.0
    else:
        f1 = 2 * precision * recall / (precision + recall)
    return {
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "precision": precision,
        "recall": recall,
        "f1": f1,
    }


def pair_scores(truth: set[tuple[int, int]], scored: dict[tuple[int, int], float]) -> tuple[np.ndarray, np.ndarray]:
    """Return (labels, scores) over the union of truth and scored candidates.

    Cells that are neither true nor scored are excluded: the head only emits a
    candidate where its probability exceeded the dump threshold, and treating
    every unlisted cell as a confident negative would inflate AUPRC.
    """
    keys = set(truth) | set(scored)
    y = np.array([1 if k in truth else 0 for k in keys], dtype=np.int8)
    s = np.array([scored.get(k, 0.0) for k in keys], dtype=np.float64)
    return y, s


def average_precision(y: np.ndarray, s: np.ndarray) -> float:
    """Step-wise AP, the same definition sklearn uses for AUPRC."""
    if y.sum() == 0:
        return float("nan")
    order = np.argsort(-s, kind="mergesort")
    y = y[order]
    tp = np.cumsum(y)
    precision = tp / np.arange(1, len(y) + 1)
    recall = tp / y.sum()
    ap = np.sum(np.diff(np.concatenate([[0.0], recall])) * precision)
    return float(ap)


def best_threshold(y: np.ndarray, s: np.ndarray, steps: int = 200) -> tuple[float, float]:
    """Oracle threshold: the cut on this very set that maximises F1."""
    if len(s) == 0:
        return 0.5, float("nan")
    cuts = np.unique(np.quantile(s, np.linspace(0, 1, steps)))
    best = (0.5, -1.0)
    for c in cuts:
        pred = s >= c
        tp = int((pred & (y == 1)).sum())
        fp = int((pred & (y == 0)).sum())
        fn = int((~pred & (y == 1)).sum())
        f1 = prf(tp, fp, fn)["f1"]
        if f1 > best[1]:
            best = (float(c), f1)
    return best


def summarise(
    rows: list[dict],
    kept_lo: int,
    kept_hi: int,
    threshold: float = 0.5,
) -> dict:
    """Aggregate one condition's rows into the D1/D2 table entries.

    Each row is ``{"truth": [[i,j],...], "scored": {"i,j": prob}, "length": L}``
    with pairs sorted and i<j.
    """
    region = Region(kept_lo, kept_hi)
    tot = {"tp": 0, "fp": 0, "fn": 0}
    by_region: dict[str, dict[str, int]] = {}
    by_dist: dict[str, dict[str, int]] = {}
    ys, ss = [], []
    pred_counts, true_counts = [], []

    for row in rows:
        truth = {tuple(p) for p in row["truth"]}
        scored = {tuple(int(x) for x in k.split(",")): float(v) for k, v in row["scored"].items()}
        pred = {k for k, v in scored.items() if v >= threshold}
        pred_counts.append(len(pred))
        true_counts.append(len(truth))

        for name, tp_set in (("tp", truth & pred), ("fp", pred - truth), ("fn", truth - pred)):
            tot[name] += len(tp_set)
            for i, j in tp_set:
                key = region.where(i, j)
                by_region.setdefault(key, {"tp": 0, "fp": 0, "fn": 0})[name] += 1
                d = distance_bin(i, j)
                by_dist.setdefault(d, {"tp": 0, "fp": 0, "fn": 0})[name] += 1

        y, s = pair_scores(truth, scored)
        ys.append(y)
        ss.append(s)

    y_all = np.concatenate(ys)
    s_all = np.concatenate(ss)
    oracle_cut, oracle_f1 = best_threshold(y_all, s_all)
    return {
        "n_sequences": len(rows),
        "at_threshold": round(threshold, 3),
        **{k: int(v) for k, v in tot.items()},
        **{k: (round(v, 4) if isinstance(v, float) else v) for k, v in prf(**tot).items()},
        "mean_predicted_pairs": round(float(np.mean(pred_counts)), 2),
        "mean_true_pairs": round(float(np.mean(true_counts)), 2),
        "auprc": round(average_precision(y_all, s_all), 4),
        "oracle_threshold": round(oracle_cut, 4),
        "oracle_f1": round(oracle_f1, 4),
        "by_region": {k: {kk: (round(vv, 4) if isinstance(vv, float) else vv) for kk, vv in prf(**v).items()} for k, v in sorted(by_region.items())},
        "by_distance": {k: {kk: (round(vv, 4) if isinstance(vv, float) else vv) for kk, vv in prf(**v).items()} for k, v in sorted(by_dist.items())},
    }
