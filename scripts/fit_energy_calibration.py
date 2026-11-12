#!/usr/bin/env python3
"""Fit the A-energy-fixed pseudo-energy calibration (Amendment 4 §A4.2).

The decoder turns a head probability into a ViennaRNA soft constraint as

    dG_ij = -(s * logit_ij + b)

and H2.4 is confirmatory only if ``(s, b)`` are fitted on the **frozen train
split** and then frozen. This script does that fit by 5-fold cross-validation,
maximising pooled family-macro F1 (the same metric used everywhere else), and
writes the chosen values with enough provenance to be quoted.

The CV estimate is used as a hyper-parameter selector, which is what it is here:
nothing is trained inside the folds because the head is already trained. For
each ``(s, b)`` on the grid, every sequence is decoded exactly once; its fold
assignment is fixed by ``--seed`` and by the record id, so the folds are
reproducible without storing them.

Cost: ``n_sequences * len(grid)`` ViennaRNA decodes. At ~0.2 s per decode and
17.7 k train sequences, a 5 x 5 grid is about 12 CPU-hours, so run it as a
batch job with ``--workers`` equal to the allocation.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import multiprocessing
import os
import sys
import zlib
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rff.analysis.labels import max_nested_subset  # noqa: E402
from rff.analysis.scoring import normalise_scored  # noqa: E402
from rff.decoders.nn_dp import NEG, logit_score, mask_allowed  # noqa: E402


def read_jsonl(path: Path) -> list[dict]:
    with open(path) as fh:
        return [json.loads(line) for line in fh if line.strip()]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fold_of(record_id: str, k: int) -> int:
    return zlib.crc32(record_id.encode()) % k


def dense_scores(scored: dict, seq: str) -> np.ndarray:
    n = len(seq)
    prob = np.zeros((n, n))
    present = np.zeros((n, n), dtype=bool)
    for (i, j), p in scored.items():
        if i < j:
            prob[i, j] = prob[j, i] = p
            present[i, j] = present[j, i] = True
    score = logit_score(prob)
    allowed = mask_allowed(seq)
    return np.where(present & allowed & np.triu(np.ones((n, n), dtype=bool), 1), score, NEG)


_ITEMS: list[dict] = []


def _init(items: list[dict]) -> None:
    global _ITEMS
    _ITEMS = items


def _decode_one(task: tuple[int, float, float]) -> tuple[int, float, set[tuple[int, int]]]:
    from rff.decoders.vienna_soft import ViennaSoftDecoder

    idx, scale, offset = task
    item = _ITEMS[idx]
    scores = dense_scores(item["scored"], item["seq"])
    pairs, _ = ViennaSoftDecoder(scale=scale, offset=offset).decode(item["seq"], scores)
    return idx, item["family"], set(pairs)


def family_macro_f1(preds: dict[int, set], items: list[dict], idxs: list[int]) -> float:
    """Mean of per-family mean per-sequence F1 (the project's primary metric).

    Pooling tp/fp/fn per family would weight long sequences more heavily and is
    not the metric used by the D1/P2 tables. The earlier implementation pooled
    counts; it is corrected here before any H2.4 number is produced.
    """
    per_family: dict[str, list[float]] = {}
    for idx in idxs:
        item = items[idx]
        pred = preds[idx]
        truth = item["truth"]
        tp = len(pred & truth)
        fp = len(pred - truth)
        fn = len(truth - pred)
        if tp + fp + fn == 0:
            f1 = 1.0
        elif tp == 0:
            f1 = 0.0
        else:
            f1 = 2 * tp / (2 * tp + fp + fn)
        per_family.setdefault(item["family"], []).append(f1)
    means = [float(np.mean(values)) for values in per_family.values()]
    return float(np.mean(means)) if means else float("nan")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--train-dump", required=True, help="per-sequence dump on the frozen train split")
    ap.add_argument("--train-split", required=True, help="jsonl with id -> seq")
    ap.add_argument("--grid-s", default="0.25,0.5,1.0,2.0,4.0")
    ap.add_argument(
        "--grid-b",
        default="-6,-3,0,3,6",
        help="coarse offsets, Amendment 7 section A7.4 (±2 ln 20 around zero)",
    )
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--max-seqs", type=int, default=0, help="0 = all records")
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument("--label-version", choices=("v1", "v3"), default="v3")
    ap.add_argument("--min-true-pairs", type=int, default=4)
    ap.add_argument("--truth-mode", choices=("full", "nested"), default="nested")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    seqs = {r["id"]: r["seq"] for r in read_jsonl(ROOT / args.train_split)}
    rows = read_jsonl(Path(args.train_dump))
    if args.max_seqs and len(rows) > args.max_seqs:
        keep = rng.choice(len(rows), size=args.max_seqs, replace=False)
        rows = [rows[i] for i in sorted(keep)]

    items = []
    n_below_min = 0
    n_parse_failed = 0
    truth_key = "truth" if args.label_version == "v1" else "truth_v3"
    for row in rows:
        seq = seqs.get(row["id"])
        if seq is None:
            continue
        raw_truth = row[truth_key]
        if raw_truth is None:
            n_parse_failed += 1
            continue
        if len(raw_truth) < args.min_true_pairs:
            n_below_min += 1
            continue
        truth = {tuple(p) for p in raw_truth}
        items.append(
            {
                "id": row["id"],
                "family": row["family"],
                "seq": seq,
                "scored": normalise_scored(row["scored"]),
                "truth": max_nested_subset(truth) if args.truth_mode == "nested" else truth,
            }
        )
    if not items:
        raise SystemExit("no usable records: check --train-dump / --train-split")

    folds = [[i for i, it in enumerate(items) if fold_of(it["id"], args.folds) == k]
             for k in range(args.folds)]

    grid = [(float(s), float(b)) for s in args.grid_s.split(",") for b in args.grid_b.split(",")]
    results = []
    global _ITEMS
    _ITEMS = items
    for scale, offset in grid:
        preds: dict[int, set] = {}
        # Linux fork lets workers inherit the read-only item table. Passing it
        # through initializer arguments pickled the full table for every grid
        # point and drove MaxRSS into the hundreds of GB; this avoids that.
        try:
            ctx = multiprocessing.get_context("fork")
            executor = ProcessPoolExecutor(max_workers=args.workers, mp_context=ctx)
        except ValueError:
            # Windows/local fallback: spawn workers cannot inherit _ITEMS.
            executor = ProcessPoolExecutor(
                max_workers=args.workers, initializer=_init, initargs=(items,)
            )
        with executor as ex:
            for idx, _family, pairs in ex.map(
                _decode_one, [(i, scale, offset) for i in range(len(items))], chunksize=16
            ):
                preds[idx] = pairs
        fold_f1 = [family_macro_f1(preds, items, idxs) for idxs in folds if idxs]
        results.append(
            {
                "s": scale,
                "b": offset,
                "folds": [round(f, 4) for f in fold_f1],
                "mean_f1": round(float(np.mean(fold_f1)), 4),
                "sd_f1": round(float(np.std(fold_f1)), 4),
                "pooled_f1": round(
                    family_macro_f1(preds, items, list(range(len(items)))), 4
                ),
            }
        )

    best = max(results, key=lambda r: r["mean_f1"])
    report = {
        "status": "fitted; freeze before scoring id_test or withheld (A4.2)",
        "n_sequences": len(items),
        "label_policy": {
            "version": args.label_version,
            "min_true_pairs": args.min_true_pairs,
            "truth_mode": args.truth_mode,
            "n_excluded_below_min": n_below_min,
            "n_excluded_label_parse_failure": n_parse_failed,
        },
        "folds": args.folds,
        "seed": args.seed,
        "max_seqs": args.max_seqs or None,
        "grid": {"s": [float(x) for x in args.grid_s.split(",")],
                 "b": [float(x) for x in args.grid_b.split(",")]},
        "metric": (
            "family-macro F1 (mean of per-family mean per-sequence F1), "
            "held-out folds; selection by mean over folds"
        ),
        "selected": {"s": best["s"], "b": best["b"], "mean_f1": best["mean_f1"]},
        "all": results,
        "inputs": {
            "train_dump_sha256": sha256_file(Path(args.train_dump)),
            "train_split_sha256": sha256_file(ROOT / args.train_split),
        },
        "note": (
            "A4.2 requires these values to be frozen and hashed before id_test or "
            "withheld is scored; the offset b absorbs the pos_weight-20 shift of "
            "the logits so that s is not absorbing calibration"
        ),
    }
    out = ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"selected": report["selected"], "out": str(out)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
