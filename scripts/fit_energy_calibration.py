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
import os
import sys
import zlib
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

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
    pairs, _ = ViennaSoftDecoder(scale=scale, offset=offset).decode(item["seq"], item["scores"])
    return idx, item["family"], set(pairs)


def family_macro_f1(preds: dict[int, set], items: list[dict], idxs: list[int]) -> float:
    per_family: dict[str, list[int]] = {}
    for idx in idxs:
        item = items[idx]
        pred = preds[idx]
        truth = item["truth"]
        tp = len(pred & truth)
        fp = len(pred - truth)
        fn = len(truth - pred)
        rec = per_family.setdefault(item["family"], [0, 0, 0])
        rec[0] += tp
        rec[1] += fp
        rec[2] += fn
    f1s = []
    for tp, fp, fn in per_family.values():
        denom = 2 * tp + fp + fn
        f1s.append(2 * tp / denom if denom else 0.0)
    return float(np.mean(f1s)) if f1s else float("nan")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--train-dump", required=True, help="per-sequence dump on the frozen train split")
    ap.add_argument("--train-split", required=True, help="jsonl with id -> seq")
    ap.add_argument("--grid-s", default="0.25,0.5,1.0,2.0,4.0")
    ap.add_argument("--grid-b", default="-3,-1.5,0,1.5,3")
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--max-seqs", type=int, default=0, help="0 = all records")
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    seqs = {r["id"]: r["seq"] for r in read_jsonl(ROOT / args.train_split)}
    rows = read_jsonl(Path(args.train_dump))
    if args.max_seqs and len(rows) > args.max_seqs:
        keep = rng.choice(len(rows), size=args.max_seqs, replace=False)
        rows = [rows[i] for i in sorted(keep)]

    items = []
    for row in rows:
        seq = seqs.get(row["id"])
        if seq is None:
            continue
        items.append(
            {
                "id": row["id"],
                "family": row["family"],
                "seq": seq,
                "scores": dense_scores(normalise_scored(row["scored"]), seq),
                "truth": {tuple(p) for p in row["truth"]},
            }
        )
    if not items:
        raise SystemExit("no usable records: check --train-dump / --train-split")

    folds = [[i for i, it in enumerate(items) if fold_of(it["id"], args.folds) == k]
             for k in range(args.folds)]

    grid = [(float(s), float(b)) for s in args.grid_s.split(",") for b in args.grid_b.split(",")]
    results = []
    for scale, offset in grid:
        preds: dict[int, set] = {}
        with ProcessPoolExecutor(max_workers=args.workers, initializer=_init, initargs=(items,)) as ex:
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
        "folds": args.folds,
        "seed": args.seed,
        "max_seqs": args.max_seqs or None,
        "grid": {"s": [float(x) for x in args.grid_s.split(",")],
                 "b": [float(x) for x in args.grid_b.split(",")]},
        "metric": "pooled family-macro F1, held-out folds; selection by mean over folds",
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
