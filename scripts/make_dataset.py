#!/usr/bin/env python3
"""Turn a labelled source into the JSONL form the trainers read.

Pairs are filtered to the project policy: canonical AU/UA/GC/CG plus GU/UG.
Everything else (non-canonical pairs, pseudoknots expressed with other bracket
types) is dropped, which keeps the label space identical to the bpRNA split
builder.

``--source archiveii`` splits by family prefix (ArchiveII record ids start with
the family name) so a held-out split is family-level, not record-level.
"""

from __future__ import annotations

import argparse
import json
import random
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

CANONICAL = {("A", "U"), ("U", "A"), ("G", "C"), ("C", "G"), ("G", "U"), ("U", "G")}


def filter_pairs(seq: str, pairs) -> list[list[int]]:
    keep = set()
    for pair in pairs:
        i, j = int(pair[0]), int(pair[1])
        if i > j:
            i, j = j, i
        if 0 <= i < j < len(seq) and (seq[i], seq[j]) in CANONICAL:
            keep.add((i, j))
    return [list(p) for p in sorted(keep)]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", default="archiveii", choices=["archiveii"])
    ap.add_argument("--data-json", default="data/raw/archiveII/data.json")
    ap.add_argument("--out", default="data/processed/archiveii")
    ap.add_argument("--val-fraction", type=float, default=0.15)
    ap.add_argument("--seed", type=int, default=17)
    args = ap.parse_args()

    rng = random.Random(args.seed)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    data = json.loads((ROOT / args.data_json).read_text())
    records = []
    dropped = Counter()
    for rid, rec in data.items():
        seq = rec["sequence"].upper().replace("T", "U")
        pairs = filter_pairs(seq, rec.get("structure", []))
        if not pairs:
            dropped["no_pairs"] += 1
        records.append(
            {
                "id": rid,
                "seq": seq,
                "pairs": pairs,
                "family": rid.split("_", 1)[0],
                "source": "archiveii",
            }
        )

    fams = sorted({r["family"] for r in records})
    val_fams = set(rng.sample(fams, max(1, int(round(len(fams) * args.val_fraction)))))
    train = [r for r in records if r["family"] not in val_fams]
    val = [r for r in records if r["family"] in val_fams]

    for name, recs in (("train", train), ("val", val)):
        path = out / f"{name}.jsonl"
        with open(path, "w") as fh:
            for rec in recs:
                fh.write(json.dumps(rec) + "\n")
    summary = {
        "n": len(records),
        "val_families": sorted(val_fams),
        "counts": {"train": len(train), "val": len(val)},
        "dropped": dict(dropped),
        "mean_pairs": round(
            sum(len(r["pairs"]) for r in records) / max(1, len(records)), 2
        ),
        "mean_len": round(sum(len(r["seq"]) for r in records) / max(1, len(records)), 1),
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
