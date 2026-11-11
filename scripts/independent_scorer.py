#!/usr/bin/env python3
"""An independent base-pair scorer, written without reusing project code.

PI Round 4e: a scorer bug would silently corrupt every number from here on, so
the scorer needs a second implementation that agrees per sequence. This file
deliberately shares nothing with `rff.analysis.scoring`:

* the dot-bracket parser scans backwards for each closing bracket instead of
  using a forward stack;
* counts are maintained as exact `Fraction` arithmetic and F1 is computed from
  the confusion counts directly.

Only the *interface* is shared: `pairs_from_db`, `counts`, `f1_from_counts`.
`scripts/scorer_crosscheck.py` compares the two implementations per sequence.
"""

from __future__ import annotations

from fractions import Fraction


def pairs_from_db(db: str) -> set[tuple[int, int]]:
    """Parse dot-bracket by matching each ')' with the nearest unclosed '('."""
    opens: list[int] = []
    pairs: set[tuple[int, int]] = set()
    for pos, sym in enumerate(db):
        if sym == "(":
            opens.append(pos)
        elif sym == ")":
            if not opens:
                raise ValueError("closing bracket without an opener")
            pairs.add((opens.pop(), pos))
        elif sym not in ".-_":
            raise ValueError(f"unexpected symbol {sym!r} at {pos}")
    if opens:
        raise ValueError("unclosed opening bracket")
    return pairs


def counts(pred: set[tuple[int, int]], truth: set[tuple[int, int]]) -> tuple[int, int, int]:
    tp = len(pred.intersection(truth))
    fp = len(pred.difference(truth))
    fn = len(truth.difference(pred))
    return tp, fp, fn


def f1_from_counts(tp: int, fp: int, fn: int) -> Fraction:
    if tp + fp + fn == 0:
        return Fraction(1)
    denominator = 2 * tp + fp + fn
    if denominator == 0:
        return Fraction(0)
    return Fraction(2 * tp, denominator)


def family_macro_f1(per_family: dict[str, list[int]]) -> Fraction:
    if not per_family:
        return Fraction(0)
    total = sum(
        (f1_from_counts(*triple) for triple in per_family.values()),
        start=Fraction(0),
    )
    return total / len(per_family)


def shift_pairs(pairs: set[tuple[int, int]], delta: int) -> set[tuple[int, int]]:
    """The ±1 slip diagnostic: move every pair index by ``delta``."""
    out = set()
    for i, j in pairs:
        a, b = i + delta, j + delta
        if a >= 0 and b > a:
            out.add((a, b))
    return out
