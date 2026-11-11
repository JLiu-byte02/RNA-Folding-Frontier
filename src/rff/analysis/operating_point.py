"""Operating-point selection for the H2.4 comparison (Amendment 4 §A4.1).

The rule is deliberately small and pure: it takes the scanned greedy
``(threshold, recall, pair_count)`` triples and decoder A's pooled numbers, and
returns which matched comparison applies. Keeping it out of the scoring script
means the frozen order can be unit-tested without inventing a toy RNA whose
greedy and DP pair counts happen to differ by the right amount.

Order, frozen before any data existed:

1. ``recall_matched`` — a threshold whose pooled recall equals A's within
   ``recall_tol``; otherwise the closest point **at or above** A's recall.
2. ``pair_budget`` — only if (1) is infeasible: the threshold whose pair count is
   closest to A's, feasible within ``pair_tol_frac`` of A's pair count.
3. ``None`` — neither is feasible; the caller reports ``inconclusive`` and adds
   no third fallback.
"""

from __future__ import annotations

from typing import Iterable, NamedTuple


class Scan(NamedTuple):
    threshold: float
    recall: float
    pairs: int


def select_operating_point(
    scan: Iterable[Scan],
    dp_recall: float,
    dp_pairs: int,
    recall_tol: float = 0.002,
    pair_tol_frac: float = 0.01,
) -> dict:
    rows = list(scan)
    if not rows:
        raise ValueError("empty threshold scan")

    primary = None
    for row in rows:
        if abs(row.recall - dp_recall) <= recall_tol:
            gap = abs(row.recall - dp_recall)
            if primary is None or gap < primary["gap"]:
                primary = {
                    "gap": gap,
                    "threshold": row.threshold,
                    "recall": row.recall,
                    "pairs": row.pairs,
                    "matches": "exact",
                }
    if primary is None:
        above = [r for r in rows if r.recall >= dp_recall]
        if above:
            row = min(above, key=lambda r: r.recall)
            primary = {
                "gap": row.recall - dp_recall,
                "threshold": row.threshold,
                "recall": row.recall,
                "pairs": row.pairs,
                "matches": "nearest_above",
            }

    budget_row = min(rows, key=lambda r: abs(r.pairs - dp_pairs))
    pair_tol = max(1.0, pair_tol_frac * dp_pairs)
    budget_gap = abs(budget_row.pairs - dp_pairs)
    budget_feasible = bool(budget_gap <= pair_tol)

    if primary is not None:
        basis = "recall_matched"
    elif budget_feasible:
        basis = "pair_budget"
    else:
        basis = None

    return {
        "basis": basis,
        "recall_matched": primary,
        "pair_budget": {
            "threshold": budget_row.threshold,
            "recall": budget_row.recall,
            "pairs": budget_row.pairs,
            "gap_pairs": budget_gap,
            "tolerance_pairs": pair_tol,
            "feasible": budget_feasible,
        },
    }


__all__ = ["Scan", "select_operating_point"]
