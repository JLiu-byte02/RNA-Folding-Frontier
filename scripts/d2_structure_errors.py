#!/usr/bin/env python3
"""D2 — H2.1 and H2.2: are long-range false positives "spurious partners"?

H2 was pre-registered in Amendment 2 §A2.4. This script computes its first two
tests exactly as written, on the native head (condition a) dumps:

**H2.1** raw-threshold false-positive rate at ``|i-j| > 128`` **per sequence**;
passes if the mean rate for ``L >= 150`` is at least twice the mean rate for
``80 <= L < 150``.

**H2.2** among those long-range false positives, the fraction whose pair lies in
a stretch of at least 4 contiguous stacked canonical or GU pairs, compared with
a null of random canonical, non-true pairs matched on distance. Passes if the
observed fraction is at least twice the null.

Both are reported with a family-block bootstrap interval, because sequences
inside a family are not independent — the same reason the P2 table resamples
families. The ratio is the statistic; the interval is on the ratio, not on each
arm separately.

The decoder here is the bare threshold, because H2's wording is "raw-threshold".
That is not the production decoder and must not be quoted as its error profile;
the production-decoder splits are printed separately and labelled exploratory.

H2.3 (the block-boundary test) needs P2-block caches and is not computed here.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rff.analysis.scoring import normalise_scored, require_finite  # noqa: E402

CANONICAL = set("AU UA GC CG GU UG".split())
LONG_RANGE = 128


def read_jsonl(path: Path) -> list[dict]:
    with open(path) as fh:
        return [json.loads(line) for line in fh if line.strip()]


def helix_length(seq: str, i: int, j: int, cap: int = 40) -> int:
    """Length of the contiguous canonical stack containing (i, j)."""
    if seq[i] + seq[j] not in CANONICAL:
        return 0
    back = 0
    while (
        back < cap
        and i - back - 1 >= 0
        and j + back + 1 < len(seq)
        and seq[i - back - 1] + seq[j + back + 1] in CANONICAL
    ):
        back += 1
    fwd = 0
    while (
        fwd < cap
        and i + fwd + 1 < j - fwd - 1
        and seq[i + fwd + 1] + seq[j - fwd - 1] in CANONICAL
    ):
        fwd += 1
    return 1 + back + fwd


def null_pair(
    seq: str, d: int, truth: set[tuple[int, int]], rng: np.random.Generator, tries: int = 20
) -> tuple[int, int] | None:
    """Random canonical, non-true pair with the same distance."""
    if d <= 0 or d >= len(seq):
        return None
    for _ in range(tries):
        i = int(rng.integers(0, len(seq) - d))
        j = i + d
        if (i, j) in truth or (j, i) in truth:
            continue
        if seq[i] + seq[j] in CANONICAL:
            return i, j
    return None


def bootstrap_ratio(
    per_family: dict[str, dict[str, float]], n_boot: int, seed: int
) -> dict:
    """Bootstrap over families: ratio of pooled observed to pooled null rates.

    ``per_family[f] = {"hits_obs", "n_obs", "hits_null", "n_null"}``.
    """
    fams = sorted(per_family)
    if len(fams) < 3:
        return {"ratio": None, "lo": None, "hi": None, "n_families": len(fams)}

    def ratio_of(sel) -> float:
        obs_n = sum(per_family[f]["n_obs"] for f in sel)
        obs_h = sum(per_family[f]["hits_obs"] for f in sel)
        null_n = sum(per_family[f]["n_null"] for f in sel)
        null_h = sum(per_family[f]["hits_null"] for f in sel)
        if obs_n == 0 or null_n == 0 or null_h == 0:
            return float("nan")
        return (obs_h / obs_n) / (null_h / null_n)

    point = ratio_of(fams)
    rng = np.random.default_rng(seed)
    draws = []
    for _ in range(n_boot):
        idx = rng.integers(0, len(fams), len(fams))
        draws.append(ratio_of([fams[k] for k in idx]))
    draws = np.array([d for d in draws if np.isfinite(d)])
    if draws.size == 0:
        return {"ratio": None, "lo": None, "hi": None, "n_families": len(fams)}
    return {
        "ratio": round(float(point), 4) if np.isfinite(point) else None,
        "lo": round(float(np.percentile(draws, 2.5)), 4),
        "hi": round(float(np.percentile(draws, 97.5)), 4),
        "n_families": len(fams),
    }


def bootstrap_mean_ratio(
    groups: dict[str, dict[str, float]], n_boot: int, seed: int
) -> dict:
    """Bootstrap over families: ratio of mean(long) to mean(short)."""
    fams = sorted(set(groups["short"]) | set(groups["long"]))
    if len(fams) < 3:
        return {"ratio": None, "lo": None, "hi": None, "n_families": len(fams)}

    def ratio_of(sel) -> float:
        s = [groups["short"][f] for f in sel if f in groups["short"]]
        l = [groups["long"][f] for f in sel if f in groups["long"]]
        if not s or not l or np.mean(s) == 0:
            return float("nan")
        return float(np.mean(l) / np.mean(s))

    point = ratio_of(fams)
    rng = np.random.default_rng(seed)
    draws = np.array(
        [ratio_of([fams[k] for k in rng.integers(0, len(fams), len(fams))]) for _ in range(n_boot)]
    )
    draws = draws[np.isfinite(draws)]
    return {
        "ratio": round(float(point), 4) if np.isfinite(point) else None,
        "lo": round(float(np.percentile(draws, 2.5)), 4) if draws.size else None,
        "hi": round(float(np.percentile(draws, 97.5)), 4) if draws.size else None,
        "n_families": len(fams),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump-dir", required=True)
    ap.add_argument("--names", nargs="*", default=[])
    ap.add_argument("--split", required=True, help="jsonl with id -> seq")
    ap.add_argument("--threshold", type=float, default=0.5)
    ap.add_argument("--min-helix", type=int, default=4)
    ap.add_argument("--n-boot", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--null-reps", type=int, default=1)
    ap.add_argument("--label-version", choices=("v1", "v3"), default="v1")
    ap.add_argument(
        "--min-true-pairs",
        type=int,
        default=0,
        help="exclude rows with fewer true pairs; 4 is Amendment 8 primary",
    )
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    seqs = {r["id"]: r["seq"] for r in read_jsonl(ROOT / args.split)}
    dump_dir = Path(args.dump_dir)
    files = (
        [dump_dir / f"{n}.jsonl" for n in args.names]
        if args.names
        else sorted(dump_dir.glob("*.jsonl"))
    )
    rng = np.random.default_rng(args.seed)

    report: dict[str, dict] = {}
    for f in files:
        rows = read_jsonl(f)
        if not rows or "truth" not in rows[0] or "scored" not in rows[0]:
            continue
        short: dict[str, list[float]] = defaultdict(list)  # 80-149
        long: dict[str, list[float]] = defaultdict(list)  # >=150
        short_per_100: dict[str, list[float]] = defaultdict(list)
        long_per_100: dict[str, list[float]] = defaultdict(list)
        short_per_1k_cells: dict[str, list[float]] = defaultdict(list)
        long_per_1k_cells: dict[str, list[float]] = defaultdict(list)
        h22: dict[str, dict[str, float]] = defaultdict(
            lambda: {"hits_obs": 0, "n_obs": 0, "hits_null": 0, "n_null": 0}
        )
        fp_dist: dict[str, int] = defaultdict(int)
        region_counts: dict[str, dict[str, int]] = defaultdict(
            lambda: {"tp": 0, "fp": 0, "fn": 0}
        )
        n_skipped = 0
        n_excluded_below_min = 0
        n_excluded_parse = 0
        truth_key = "truth" if args.label_version == "v1" else "truth_v3"
        if truth_key not in rows[0]:
            raise ValueError(f"{f.name}: missing {truth_key}; regenerate the dump")
        for row in rows:
            seq = seqs.get(row["id"])
            if seq is None:
                n_skipped += 1
                continue
            raw_truth = row[truth_key]
            if raw_truth is None:
                n_excluded_parse += 1
                continue
            if len(raw_truth) < args.min_true_pairs:
                n_excluded_below_min += 1
                continue
            truth = {tuple(p) for p in raw_truth}
            scored = normalise_scored(row["scored"])
            require_finite(scored, where=f"{row['id']}")
            pred = {k for k, v in scored.items() if v >= args.threshold}
            fam = row["family"]
            long_fps = [p for p in (pred - truth) if abs(p[0] - p[1]) > LONG_RANGE]
            if 80 <= len(seq) < 150:
                short[fam].append(len(long_fps))
                short_per_100[fam].append(100.0 * len(long_fps) / len(seq))
                n_long_cells = max(0, (len(seq) - 128) * (len(seq) - 129) // 2)
                short_per_1k_cells[fam].append(
                    1000.0 * len(long_fps) / max(1, n_long_cells)
                )
            elif len(seq) >= 150:
                long[fam].append(len(long_fps))
                long_per_100[fam].append(100.0 * len(long_fps) / len(seq))
                n_long_cells = max(0, (len(seq) - 128) * (len(seq) - 129) // 2)
                long_per_1k_cells[fam].append(
                    1000.0 * len(long_fps) / max(1, n_long_cells)
                )

            kept_lo, kept_hi = row["kept"]
            for i, j in pred - truth:
                d = abs(i - j)
                key = "d<32" if d < 32 else ("d32-128" if d <= 128 else "d>128")
                fp_dist[key] += 1
                in_lo = kept_lo <= i < kept_hi
                in_hi = kept_lo <= j < kept_hi
                where = (
                    "both_kept"
                    if in_lo and in_hi
                    else ("one_kept" if in_lo or in_hi else "neither_kept")
                )
                region_counts[where]["fp"] += 1
                if d > LONG_RANGE:
                    h22[fam]["n_obs"] += 1
                    if helix_length(seq, i, j) >= args.min_helix:
                        h22[fam]["hits_obs"] += 1
                    for _ in range(args.null_reps):
                        pair = null_pair(seq, d, truth, rng)
                        if pair is None:
                            continue
                        h22[fam]["n_null"] += 1
                        if helix_length(seq, *pair) >= args.min_helix:
                            h22[fam]["hits_null"] += 1
            for i, j in truth & pred:
                in_lo = kept_lo <= i < kept_hi
                in_hi = kept_lo <= j < kept_hi
                where = (
                    "both_kept"
                    if in_lo and in_hi
                    else ("one_kept" if in_lo or in_hi else "neither_kept")
                )
                region_counts[where]["tp"] += 1
            for i, j in truth - pred:
                in_lo = kept_lo <= i < kept_hi
                in_hi = kept_lo <= j < kept_hi
                where = (
                    "both_kept"
                    if in_lo and in_hi
                    else ("one_kept" if in_lo or in_hi else "neither_kept")
                )
                region_counts[where]["fn"] += 1

        def mean_rate(group: dict[str, list[float]]) -> float:
            vals = [float(np.mean(v)) for v in group.values() if v]
            return round(float(np.mean(vals)), 4) if vals else float("nan")

        short_rate = mean_rate(short)
        long_rate = mean_rate(long)
        h21_ratio = (
            round(long_rate / short_rate, 4)
            if short_rate and np.isfinite(short_rate) and short_rate > 0
            else None
        )
        h22_ratio = bootstrap_ratio(h22, args.n_boot, args.seed)
        n_obs = sum(v["n_obs"] for v in h22.values())
        n_null = sum(v["n_null"] for v in h22.values())
        short_per_100_rate = mean_rate(short_per_100)
        long_per_100_rate = mean_rate(long_per_100)
        short_cell_rate = mean_rate(short_per_1k_cells)
        long_cell_rate = mean_rate(long_per_1k_cells)

        report[f.stem] = {
            "n_sequences": len(rows) - n_excluded_below_min - n_excluded_parse - n_skipped,
            "label_policy": {
                "version": args.label_version,
                "min_true_pairs": args.min_true_pairs,
                "n_input": len(rows),
                "n_excluded_below_min": n_excluded_below_min,
                "n_excluded_label_parse_failure": n_excluded_parse,
            },
            "n_sequences_without_seq": n_skipped,
            "decoder": f"bare threshold {args.threshold}",
            "H2_1": {
                "long_range_definition": f"|i-j| > {LONG_RANGE}",
                "mean_fp_per_sequence_80_149": short_rate,
                "n_sequences_80_149": sum(len(v) for v in short.values()),
                "mean_fp_per_sequence_ge150": long_rate,
                "n_sequences_ge150": sum(len(v) for v in long.values()),
                "ratio_ge150_over_80_149": h21_ratio,
                "bootstrap_ratio": bootstrap_mean_ratio(
                    {"short": {k: float(np.mean(v)) for k, v in short.items() if v},
                     "long": {k: float(np.mean(v)) for k, v in long.items() if v}},
                    args.n_boot,
                    args.seed,
                ),
                "passes_2x_rule": bool(h21_ratio is not None and h21_ratio >= 2.0),
            },
            "exploratory_H2_1_normalised": {
                "note": (
                    "not a pre-registered gate; the raw per-sequence count grows "
                    "with the number of candidate pairs, so these normalisations "
                    "test whether the length trend survives"
                ),
                "mean_fp_per_100nt_80_149": short_per_100_rate,
                "mean_fp_per_100nt_ge150": long_per_100_rate,
                "ratio_per_100nt": (
                    round(long_per_100_rate / short_per_100_rate, 4)
                    if short_per_100_rate and np.isfinite(short_per_100_rate)
                    else None
                ),
                "mean_fp_per_1000_long_range_cells_80_149": short_cell_rate,
                "mean_fp_per_1000_long_range_cells_ge150": long_cell_rate,
                "ratio_per_1000_long_range_cells": (
                    round(long_cell_rate / short_cell_rate, 4)
                    if short_cell_rate and np.isfinite(short_cell_rate)
                    else None
                ),
            },
            "H2_2": {
                "min_helix": args.min_helix,
                "n_long_range_fp": n_obs,
                "observed_frac": round(
                    sum(v["hits_obs"] for v in h22.values()) / n_obs, 4
                )
                if n_obs
                else None,
                "n_null_pairs": n_null,
                "null_frac": round(
                    sum(v["hits_null"] for v in h22.values()) / n_null, 4
                )
                if n_null
                else None,
                "bootstrap_ratio": h22_ratio,
                "ratio": h22_ratio["ratio"],
                "passes_2x_rule": bool(
                    h22_ratio["ratio"] is not None and h22_ratio["ratio"] >= 2.0
                ),
                "null_definition": (
                    "random canonical pair with the same |i-j| inside the same "
                    "sequence, not present in the true structure"
                ),
            },
            "exploratory_fp_by_distance": dict(sorted(fp_dist.items())),
            "exploratory_region_counts": {
                k: dict(v) for k, v in sorted(region_counts.items())
            },
        }
        print(json.dumps({f.stem: {"H2_1": report[f.stem]["H2_1"]["ratio_ge150_over_80_149"],
                                   "H2_2": report[f.stem]["H2_2"]["ratio"]}}), flush=True)

    out = ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
