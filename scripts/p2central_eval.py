#!/usr/bin/env python3
"""Evaluate corrected P2-central dumps with three readouts and two cuts."""

from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
CONDITIONS = ("a", "b", "c", "d", "e")
SPLITS = ("id_test", "test_withheld")
GRID = tuple(round(0.20 + 0.05 * k, 2) for k in range(15))


def read_jsonl(path: Path) -> list[dict]:
    with open(path, encoding="utf-8") as fh:
        return [json.loads(line) for line in fh if line.strip()]


def f1(tp: int, fp: int, fn: int) -> float:
    if tp + fp + fn == 0:
        return 1.0
    if tp == 0:
        return 0.0
    return 2 * tp / (2 * tp + fp + fn)


def family_macro(per_seq: list[tuple[str, float]]) -> float:
    by_family: dict[str, list[float]] = defaultdict(list)
    for family, value in per_seq:
        by_family[family].append(value)
    return float(np.mean([np.mean(v) for v in by_family.values()])) if by_family else float("nan")


def metrics(rows: list[dict], cut: float, kept_only: bool) -> dict:
    per_seq = []
    total = [0, 0, 0]
    n_excluded_empty_kept = 0
    for row in rows:
        truth = row["_truth"]
        pred = row["_pred"].get(cut)
        if pred is None:
            raise KeyError(f"cut {cut} not decoded")
        kept = row["kept"]
        if kept_only:
            truth = {p for p in truth if kept[0] <= p[0] < kept[1] and kept[0] <= p[1] < kept[1]}
            pred = {p for p in pred if kept[0] <= p[0] < kept[1] and kept[0] <= p[1] < kept[1]}
            if not truth:
                n_excluded_empty_kept += 1
                continue
        tp = len(truth & pred)
        fp = len(pred - truth)
        fn = len(truth - pred)
        total[0] += tp
        total[1] += fp
        total[2] += fn
        per_seq.append((row["family"], f1(tp, fp, fn)))
    return {
        "family_macro_f1": round(family_macro(per_seq), 6),
        "micro_f1": round(f1(*total), 6),
        "tp": total[0],
        "fp": total[1],
        "fn": total[2],
        "n_sequences": len(per_seq),
        "n_excluded_empty_kept_truth": n_excluded_empty_kept,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump-prefix", default="experiments/002_e1prime_runs/dump_p2fix_seed")
    ap.add_argument("--seeds", nargs="+", type=int, default=[17, 23, 41])
    ap.add_argument("--threshold", type=float, default=0.5)
    ap.add_argument("--min-true-pairs", type=int, default=4)
    ap.add_argument("--out", default="experiments/002_e1prime_runs/p2central_fixed.json")
    args = ap.parse_args()

    from rff.analysis.scoring import decode_pairs

    result: dict[str, dict] = {}
    for split in SPLITS:
        for cond in CONDITIONS:
            per_seed = []
            for seed in args.seeds:
                path = ROOT / f"{args.dump_prefix}{seed}/{cond}_{split}.jsonl"
                if not path.exists():
                    continue
                rows = []
                for row in read_jsonl(path):
                    raw = row.get("truth_v3")
                    if raw is None or len(raw) < args.min_true_pairs:
                        continue
                    truth = {tuple(p) for p in raw}
                    scored = {
                        tuple(int(x) for x in key.split(",")): float(value)
                        for key, value in row["scored"].items()
                    }
                    preds = {
                        cut: set(
                            decode_pairs(
                                scored, cut, int(row["length"]), partner=True, nesting=True
                            )
                        )
                        for cut in GRID
                    }
                    rows.append(
                        {
                            "family": row["family"],
                            "kept": row["kept"],
                            "length": int(row["length"]),
                            "_truth": truth,
                            "_pred": preds,
                        }
                    )
                if not rows:
                    continue
                # Oracle cut: maximize family-macro all-pairs F1 on 80-299,
                # the same primary band used by D1.
                band = [r for r in rows if 80 <= r["length"] < 300]
                oracle_rows = band or rows
                curve = []
                for cut in GRID:
                    vals = []
                    for r in oracle_rows:
                        truth = r["_truth"]
                        pred = r["_pred"][cut]
                        tp, fp, fn = len(truth & pred), len(pred - truth), len(truth - pred)
                        vals.append((r["family"], f1(tp, fp, fn)))
                    curve.append(family_macro(vals))
                best = int(np.argmax(curve))
                oracle = GRID[best]
                per_seed.append(
                    {
                        "seed": seed,
                        "fixed": {
                            "all": metrics(rows, args.threshold, False),
                            "kept_kept": metrics(rows, args.threshold, True),
                        },
                        "oracle": {
                            "cut": oracle,
                            "all": metrics(rows, oracle, False),
                            "kept_kept": metrics(rows, oracle, True),
                        },
                    }
                )
            result[f"{split}|{cond}"] = {
                "per_seed": per_seed,
                "n_seeds": len(per_seed),
            }
    out = ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"out": args.out, "blocks": len(result)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
