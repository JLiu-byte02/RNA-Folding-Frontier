#!/usr/bin/env python3
"""Fast fixed-threshold scoring for all-N / dinucleotide-shuffle controls."""

from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rff.analysis.labels import pair_prf  # noqa: E402
from rff.analysis.scoring import decode_pairs  # noqa: E402


def read_jsonl(path: Path) -> list[dict]:
    with open(path, encoding="utf-8") as fh:
        return [json.loads(line) for line in fh if line.strip()]


def score_dump(path: Path, threshold: float, min_true: int) -> dict:
    per_seq = []
    total = {"tp": 0, "fp": 0, "fn": 0}
    skipped = 0
    for row in read_jsonl(path):
        truth_raw = row.get("truth_v3")
        if truth_raw is None or len(truth_raw) < min_true:
            skipped += 1
            continue
        truth = {tuple(p) for p in truth_raw}
        scored = {
            tuple(int(x) for x in key.split(",")): float(value)
            for key, value in row["scored"].items()
        }
        pred = set(
            decode_pairs(scored, threshold, int(row["length"]), partner=True, nesting=True)
        )
        prf = pair_prf(truth, pred)
        per_seq.append((row["family"], prf["f1"]))
        total["tp"] += prf["tp"]
        total["fp"] += prf["fp"]
        total["fn"] += prf["fn"]

    by_family: dict[str, list[float]] = defaultdict(list)
    for family, f1 in per_seq:
        by_family[family].append(f1)
    macro = float(np.mean([np.mean(v) for v in by_family.values()])) if by_family else None
    denom = 2 * total["tp"] + total["fp"] + total["fn"]
    micro_f1 = 2 * total["tp"] / denom if denom else None
    return {
        "n_scored": len(per_seq),
        "n_excluded": skipped,
        "n_families": len(by_family),
        "family_macro_f1": round(macro, 6) if macro is not None else None,
        "micro_f1": round(micro_f1, 6) if micro_f1 is not None else None,
        "tp": total["tp"],
        "fp": total["fp"],
        "fn": total["fn"],
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default="experiments/002_e1prime_runs")
    ap.add_argument("--threshold", type=float, default=0.5)
    ap.add_argument("--min-true-pairs", type=int, default=4)
    ap.add_argument("--out", default="experiments/002_e1prime_runs/template_controls.json")
    args = ap.parse_args()

    root = ROOT / args.root
    results: dict[str, dict] = {}
    for seed in (17, 23, 41):
        dump_dir = root / f"dump_seqfree_seed{seed}"
        if not dump_dir.exists():
            continue
        for name in (
            "all_n_id_test",
            "all_n_test_withheld",
            "shuffle_id_test",
            "shuffle_test_withheld",
        ):
            path = dump_dir / f"{name}.jsonl"
            if path.exists():
                results[f"seed{seed}|{name}"] = score_dump(
                    path, args.threshold, args.min_true_pairs
                )
    if not results:
        raise SystemExit("no template-control dumps found")

    means: dict[str, dict] = {}
    for name in (
        "all_n_id_test",
        "all_n_test_withheld",
        "shuffle_id_test",
        "shuffle_test_withheld",
    ):
        vals = [
            v["family_macro_f1"]
            for k, v in results.items()
            if k.endswith("|" + name) and v["family_macro_f1"] is not None
        ]
        means[name] = {
            "n_seeds": len(vals),
            "mean_family_macro_f1": round(float(np.mean(vals)), 6) if vals else None,
        }

    out = {
        "description": "fixed 0.5 nested decoder; v3 primary truth; all-N/shuffle controls",
        "threshold": args.threshold,
        "min_true_pairs": args.min_true_pairs,
        "per_seed": results,
        "means": means,
    }
    path = ROOT / args.out
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"out": args.out, "means": means}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
