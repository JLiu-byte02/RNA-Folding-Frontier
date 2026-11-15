"""Anatomy of false positives: near-miss or structural?

For every false-positive pair a method emits, compute the Manhattan distance in
pair coordinates to the nearest true pair:

    d(i,j) = min over truth (i',j') of |i-i'| + |j-j'|

and bin d. A method whose FPs are mostly d<=2 locates the helix but offsets it;
a method whose FPs are mostly d>5 invents structure elsewhere. The comparison is
descriptive and uses the same cohort, labels, and threshold as the main table.

    python scripts/fp_anatomy.py --split test_withheld --dump ... --preds ... --out ...
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from rff.analysis.labels import max_nested_subset  # noqa: E402
import score_external_baseline as seb  # noqa: E402

BINS = [(1, 2), (3, 5), (6, 10), (11, 10**9)]
LABELS = ["1-2", "3-5", "6-10", ">10"]


def _resolve(path: str) -> Path:
    p = Path(path)
    return p if p.is_absolute() else ROOT / p


def parse_dump(path: Path, threshold: float) -> dict[str, set[tuple[int, int]]]:
    out = {}
    for row in seb.read_jsonl(path):
        pairs = set()
        for key, prob in (row.get("scored") or {}).items():
            if prob >= threshold:
                i, j = key.split(",")
                pairs.add((int(i), int(j)))
        out[row["id"]] = pairs
    return out


def parse_preds(path: Path) -> dict[str, set[tuple[int, int]]]:
    out = {}
    for row in seb.read_jsonl(path):
        if "pairs" in row:
            out[row["id"]] = {(int(a), int(b)) for a, b in row["pairs"]}
        elif "dbn" in row:
            from rff.data.io import pairs_from_dbn

            out[row["id"]] = set(pairs_from_dbn(row["dbn"]))
        else:
            raise SystemExit(f"{path}: row has neither dbn nor pairs")
    return out


def bin_of(d: int) -> str:
    for (lo, hi), name in zip(BINS, LABELS):
        if lo <= d <= hi:
            return name
    return "exact"


def summarise(fp_dist: list[int], n_fp: int) -> dict:
    counts = Counter(bin_of(d) for d in fp_dist)
    total = sum(counts.values()) or 1
    return {
        "n_fp": n_fp,
        "frac": {k: round(counts.get(k, 0) / total, 4) for k in LABELS},
        "near_miss_frac": round((counts.get("1-2", 0)) / total, 4),
        "mean_distance": round(sum(fp_dist) / total, 3) if fp_dist else 0.0,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--split", required=True)
    ap.add_argument("--cohort", required=True)
    ap.add_argument("--dump", required=True)
    ap.add_argument("--preds", required=True, help="baseline prediction JSONL to compare")
    ap.add_argument("--labels", default=seb.LABELS_V3)
    ap.add_argument("--threshold", type=float, default=0.5)
    ap.add_argument("--min-true-pairs", type=int, default=4)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    split_rows = {r["id"]: r for r in seb.read_jsonl(_resolve(args.split))}
    cohort = seb.read_jsonl(_resolve(args.cohort))
    labels = seb.load_labels_v3(_resolve(args.labels))
    head = parse_dump(_resolve(args.dump), args.threshold)
    base = parse_preds(_resolve(args.preds))

    per_method = {"head": {"fp_dist": [], "n_fp": 0, "n_seq_with_fp": 0},
                  "baseline": {"fp_dist": [], "n_fp": 0, "n_seq_with_fp": 0}}
    n_seq = 0
    for row in cohort:
        rid = row["id"]
        src = split_rows.get(rid)
        key = rid[5:] if rid.startswith("rfam:") else rid
        if key in labels:
            raw = labels[key]
            if raw is None:
                continue
        elif src is not None:
            raw = {(int(a), int(b)) for a, b in src["pairs"]}
        else:
            continue
        if len(raw) < args.min_true_pairs:
            continue
        if rid not in head or rid not in base:
            continue
        truth = max_nested_subset(raw)
        n_seq += 1
        for name, pred in (("head", head[rid]), ("baseline", base[rid])):
            fps = pred - truth
            if fps:
                per_method[name]["n_seq_with_fp"] += 1
            per_method[name]["n_fp"] += len(fps)
            for (i, j) in fps:
                d = min((abs(i - a) + abs(j - b) for (a, b) in truth), default=10**6)
                if d == 0:
                    continue
                per_method[name]["fp_dist"].append(d)

    out = {"split": args.split, "n_sequences": n_seq, "threshold": args.threshold,
           "methods": {}}
    for name, d in per_method.items():
        s = summarise(d["fp_dist"], d["n_fp"])
        s["n_sequences_with_fp"] = d["n_seq_with_fp"]
        s["frac_sequences_with_fp"] = round(d["n_seq_with_fp"] / n_seq, 4) if n_seq else 0
        out["methods"][name] = s
    dest = _resolve(args.out)
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps(out, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
