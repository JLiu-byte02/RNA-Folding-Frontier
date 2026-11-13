#!/usr/bin/env python3
"""Freeze the E-C1 family list before any coverage proxy is computed."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def read_rows(path: Path, min_true: int) -> dict[str, dict]:
    out = {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if not line.strip():
                continue
            row = json.loads(line)
            truth = row.get("truth_v3")
            if truth is None or len(truth) < min_true:
                continue
            out[row["id"]] = {
                "family": row["family"],
                "source": row.get("source", "?"),
            }
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump-root", default="experiments/002_e1prime_runs")
    ap.add_argument("--split", default="test_withheld")
    ap.add_argument("--seeds", nargs="+", type=int, default=[17, 23, 41])
    ap.add_argument("--min-true-pairs", type=int, default=4)
    ap.add_argument(
        "--out", default="experiments/002_e1prime_runs/ec1_families.json"
    )
    args = ap.parse_args()

    tables = []
    for seed in args.seeds:
        for arm in ("a", "b", "e"):
            path = (
                ROOT
                / args.dump_root
                / f"dump_p3_{arm}_seed{seed}"
                / f"{args.split}.jsonl"
            )
            if not path.exists():
                raise SystemExit(f"missing dump: {path}")
            tables.append(read_rows(path, args.min_true_pairs))

    ids = set(tables[0])
    for table in tables[1:]:
        ids &= set(table)

    families: dict[str, list[str]] = {}
    source = {}
    for rid in sorted(ids):
        meta = tables[0][rid]
        for table in tables[1:]:
            if table[rid]["family"] != meta["family"]:
                raise SystemExit(f"family mismatch for {rid}")
        families.setdefault(meta["family"], []).append(rid)
        source[rid] = meta["source"]

    result = {
        "split": args.split,
        "seeds": args.seeds,
        "min_true_pairs": args.min_true_pairs,
        "n_sequences": len(ids),
        "n_families": len(families),
        "families": {
            family: {
                "ids": members,
                "n": len(members),
                "sources": sorted({source[r] for r in members}),
            }
            for family, members in sorted(families.items())
        },
    }
    out = ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"out": args.out, "n_families": len(families), "n_sequences": len(ids)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
