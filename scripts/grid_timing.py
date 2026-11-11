#!/usr/bin/env python3
"""How long does one soft-constrained fold take, by length bin?

Amendment 7 §A7.4 pre-registers the `(s, b)` search but also makes its cost a
decision point: if the full training set is unaffordable, a fixed subsample is
used. That decision needs a measured number, not an estimate, so this times the
production decoder (`ViennaSoftDecoder`, which calls ``sc_add_bp`` once per
candidate cell) on real sequences from the local ArchiveII mirror.

Output: seconds per fold by bin, the implied cost of the whole grid over the
real training set, and the number of ViennaRNA constraints per call. Local
ArchiveII lengths are not the training split's lengths, so the projection is
labelled as an estimate and the bin table is the reported measurement.
"""

from __future__ import annotations

import argparse
import json
import statistics
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402

from rff.decoders.nn_dp import NEG  # noqa: E402

BINS = [(0, 80), (80, 150), (150, 300), (300, 600), (600, 2000)]
GRID_POINTS = 50  # A7.4: 5x5 coarse + 5x5 refinement


def band(n: int) -> str:
    for lo, hi in BINS:
        if lo <= n < hi:
            return f"{lo}-{hi - 1}"
    return f">={BINS[-1][1]}"


def score_matrix(seq: str, rng: np.random.Generator) -> np.ndarray:
    """A dense log-odds-like map: every canonical cell is a candidate."""
    n = len(seq)
    raw = rng.normal(0.0, 1.0, size=(n, n))
    scores = np.triu(raw, k=1)
    scores = scores + scores.T
    allowed = {("A", "U"), ("U", "A"), ("G", "C"), ("C", "G"), ("G", "U"), ("U", "G")}
    mask = np.full((n, n), NEG, dtype=float)
    for i in range(n):
        for j in range(i + 1, n):
            if (seq[i], seq[j]) in allowed:
                mask[i, j] = scores[i, j]
    return mask


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--archive", default="data/raw/archiveII/data.json")
    ap.add_argument("--per-bin", type=int, default=8)
    ap.add_argument("--limit-scan", type=int, default=0, help="0 = scan the whole file")
    ap.add_argument("--train-sequences", type=int, default=17_700)
    ap.add_argument("--seed", type=int, default=17)
    ap.add_argument("--out", default="experiments/009_grid_timing/grid_timing.json")
    args = ap.parse_args()

    from rff.decoders.vienna_soft import ViennaSoftDecoder

    data = json.loads((ROOT / args.archive).read_text(encoding="utf-8"))
    items = list(data.items())
    if args.limit_scan:
        items = items[: args.limit_scan]

    counts: dict[str, int] = {}
    for rid, rec in items:
        seq = rec.get("sequence") or ""
        if seq:
            counts[band(len(seq))] = counts.get(band(len(seq)), 0) + 1

    rng = np.random.default_rng(args.seed)
    by_bin: dict[str, list[tuple[str, str]]] = {}
    for rid, rec in items:
        seq = rec.get("sequence") or ""
        if not seq:
            continue
        by_bin.setdefault(band(len(seq)), []).append((rid, seq))
    picked: dict[str, list[tuple[str, str]]] = {}
    for name, bucket in by_bin.items():
        idx = rng.choice(len(bucket), size=min(args.per_bin, len(bucket)), replace=False)
        picked[name] = [bucket[i] for i in sorted(idx)]

    decoder = ViennaSoftDecoder(scale=1.0, offset=0.0)
    rows = []
    for name in [band(lo) for lo, _ in BINS]:
        samples = picked.get(name, [])
        if not samples:
            continue
        times = []
        constraints = []
        for rid, seq in samples:
            scores = score_matrix(seq, rng)
            start = time.perf_counter()
            decoder.decode(seq, scores)
            times.append(time.perf_counter() - start)
            constraints.append(decoder.last["n_constraints"])
        rows.append(
            {
                "band": name,
                "n_sequences_in_file": counts[name],
                "n_timed": len(times),
                "median_s": round(statistics.median(times), 4),
                "mean_s": round(statistics.fmean(times), 4),
                "p90_s": round(sorted(times)[max(0, int(0.9 * len(times)) - 1)], 4),
                "median_constraints": int(statistics.median(constraints)),
            }
        )

    projected = sum(r["median_s"] * r["n_sequences_in_file"] for r in rows)
    report = {
        "purpose": "cost of one soft-constrained fold, by length bin (A7.4 decision input)",
        "archive": args.archive,
        "grid_points": GRID_POINTS,
        "assumed_training_sequences": args.train_sequences,
        "rows": rows,
        "implied_seconds_per_full_pass_over_this_file": round(projected, 1),
        "implied_cpu_hours_per_full_pass_over_this_file": round(projected / 3600, 2),
        "implied_cpu_hours_for_the_grid_assuming_train_size": round(
            projected / max(1, sum(counts.values())) * args.train_sequences * GRID_POINTS / 3600, 2
        ),
        "caveat": (
            "ArchiveII lengths are not the training split's lengths; the bin table is the "
            "measurement and the grid projection is an estimate that must be recomputed on "
            "the cluster from the real split before the A7.4 subsample decision"
        ),
        "decoder": "rff.decoders.vienna_soft.ViennaSoftDecoder, one sc_add_bp per candidate cell",
    }
    out = ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: report[k] for k in ("implied_cpu_hours_per_full_pass_over_this_file",
                                             "implied_cpu_hours_for_the_grid_assuming_train_size")}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
