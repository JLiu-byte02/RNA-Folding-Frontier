"""Shared scoring for the P2-block and P3 train-and-test tables.

The dumps contain per-sequence probabilities; everything here is plain numpy so
the confirmatory tables can be recomputed without a GPU. Family-macro F1 is the
mean of per-family mean per-sequence F1, matching the rest of the project.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np

from rff.analysis.labels import max_nested_subset
from rff.analysis.scoring import decode_pairs, normalise_scored


@dataclass(frozen=True)
class SeqScore:
    id: str
    family: str
    tp: int
    fp: int
    fn: int

    @property
    def f1(self) -> float:
        if self.tp + self.fp + self.fn == 0:
            return 1.0
        if self.tp == 0:
            return 0.0
        return 2 * self.tp / (2 * self.tp + self.fp + self.fn)


def read_jsonl(path: Path) -> list[dict]:
    with open(path, encoding="utf-8") as fh:
        return [json.loads(line) for line in fh if line.strip()]


def score_dump(
    path: Path,
    *,
    threshold: float = 0.5,
    min_true_pairs: int = 4,
    truth_key: str = "truth_v3",
    nested: bool = False,
) -> dict[str, SeqScore]:
    """Decode one dump and return per-sequence counts keyed by raw id."""
    out: dict[str, SeqScore] = {}
    for row in read_jsonl(path):
        raw = row.get(truth_key)
        if raw is None or len(raw) < min_true_pairs:
            continue
        truth = {tuple(int(x) for x in pair) for pair in raw}
        if nested:
            truth = max_nested_subset(truth)
        scored = normalise_scored(row["scored"])
        pred = set(
            decode_pairs(
                scored,
                threshold,
                int(row["length"]),
                partner=True,
                nesting=True,
            )
        )
        out[row["id"]] = SeqScore(
            id=row["id"],
            family=row["family"],
            tp=len(truth & pred),
            fp=len(pred - truth),
            fn=len(truth - pred),
        )
    return out


def common_ids(*tables: dict[str, SeqScore]) -> set[str]:
    if not tables:
        return set()
    ids = set(tables[0])
    for table in tables[1:]:
        ids &= set(table)
    return ids


def family_values(
    scores: dict[str, SeqScore], ids: Iterable[str] | None = None
) -> dict[str, list[float]]:
    by_family: dict[str, list[float]] = {}
    for rid in ids if ids is not None else scores:
        score = scores[rid]
        by_family.setdefault(score.family, []).append(score.f1)
    return by_family


def family_means(
    scores: dict[str, SeqScore], ids: Iterable[str] | None = None
) -> dict[str, float]:
    return {
        family: float(np.mean(values))
        for family, values in family_values(scores, ids).items()
    }


def family_macro(
    scores: dict[str, SeqScore], ids: Iterable[str] | None = None
) -> float:
    values = family_means(scores, ids)
    return float(np.mean(list(values.values()))) if values else float("nan")


def micro(scores: dict[str, SeqScore], ids: Iterable[str] | None = None) -> float:
    selected = list(ids if ids is not None else scores)
    tp = sum(scores[rid].tp for rid in selected)
    fp = sum(scores[rid].fp for rid in selected)
    fn = sum(scores[rid].fn for rid in selected)
    if tp + fp + fn == 0:
        return 1.0
    if tp == 0:
        return 0.0
    return 2 * tp / (2 * tp + fp + fn)


def seed_average_family_means(
    per_seed: list[dict[str, SeqScore]], ids: Iterable[str]
) -> dict[str, float]:
    """Average each family's per-sequence-mean F1 across seeds.

    The id set must be common to every seed table; the caller is responsible for
    passing ``common_ids(*per_seed)``. Families absent from a seed are omitted
    rather than silently treated as zero.
    """
    ids = list(ids)
    by_seed = [family_means(table, ids) for table in per_seed]
    families = set(by_seed[0]) if by_seed else set()
    for table in by_seed[1:]:
        families &= set(table)
    return {
        family: float(np.mean([table[family] for table in by_seed]))
        for family in families
    }


def bootstrap_mean(
    values: dict[str, float],
    *,
    n_boot: int = 10_000,
    rng_seed: int = 20261112,
) -> dict[str, float]:
    if not values:
        return {"mean": float("nan"), "lo": float("nan"), "hi": float("nan"), "n": 0}
    families = sorted(values)
    arr = np.asarray([values[f] for f in families], dtype=float)
    rng = np.random.default_rng(rng_seed)
    draws = rng.integers(0, len(arr), size=(n_boot, len(arr)))
    means = arr[draws].mean(axis=1)
    return {
        "mean": float(arr.mean()),
        "lo": float(np.percentile(means, 2.5)),
        "hi": float(np.percentile(means, 97.5)),
        "n": int(len(arr)),
    }


def bootstrap_ratio(
    condition_values: dict[str, dict[str, float]],
    *,
    numerator: str,
    denominator_a: str,
    denominator_e: str,
    n_boot: int = 10_000,
    rng_seed: int = 20261112,
) -> dict[str, float | bool]:
    """Bootstrap ``(num - e) / (a - e)`` over families.

    The same resampled families are used for all three conditions. A ratio is
    undefined when the resampled ``a - e`` is non-positive; the fraction of such
    draws is reported instead of hiding it.
    """
    families = (
        set(condition_values[numerator])
        & set(condition_values[denominator_a])
        & set(condition_values[denominator_e])
    )
    if not families:
        return {
            "ratio": float("nan"),
            "lo": float("nan"),
            "hi": float("nan"),
            "n_families": 0,
            "undefined_fraction": 1.0,
        }
    families = sorted(families)
    num = np.asarray([condition_values[numerator][f] for f in families])
    a = np.asarray([condition_values[denominator_a][f] for f in families])
    e = np.asarray([condition_values[denominator_e][f] for f in families])
    point_den = float(a.mean() - e.mean())
    point = (
        float((num.mean() - e.mean()) / point_den) if point_den > 0 else float("nan")
    )
    rng = np.random.default_rng(rng_seed)
    draws = rng.integers(0, len(families), size=(n_boot, len(families)))
    num_m = num[draws].mean(axis=1)
    a_m = a[draws].mean(axis=1)
    e_m = e[draws].mean(axis=1)
    den = a_m - e_m
    valid = den > 0
    ratios = np.full(n_boot, np.nan)
    ratios[valid] = (num_m[valid] - e_m[valid]) / den[valid]
    return {
        "ratio": point,
        "lo": float(np.nanpercentile(ratios, 2.5)) if valid.any() else float("nan"),
        "hi": float(np.nanpercentile(ratios, 97.5)) if valid.any() else float("nan"),
        "n_families": int(len(families)),
        "undefined_fraction": float(1.0 - valid.mean()),
    }


def bootstrap_delta(
    block_values: dict[str, float],
    a_values: dict[str, float],
    *,
    n_boot: int = 10_000,
    rng_seed: int = 20261112,
) -> dict[str, float]:
    families = sorted(set(block_values) & set(a_values))
    if not families:
        return {"delta": float("nan"), "lo": float("nan"), "hi": float("nan"), "n": 0}
    delta = np.asarray([block_values[f] - a_values[f] for f in families])
    rng = np.random.default_rng(rng_seed)
    draws = rng.integers(0, len(delta), size=(n_boot, len(delta)))
    means = delta[draws].mean(axis=1)
    return {
        "delta": float(delta.mean()),
        "lo": float(np.percentile(means, 2.5)),
        "hi": float(np.percentile(means, 97.5)),
        "n": int(len(families)),
    }
