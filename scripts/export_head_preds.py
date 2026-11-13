#!/usr/bin/env python3
"""Export decoded head predictions from a p2central dump as a predictions JSONL.

The external-baseline table and the head must be scored by exactly the same
function (``scripts/score_external_baseline.py score``), so this script turns a
dump into ``{"id": ..., "pairs": [[i, j], ...]}`` lines at a fixed cut.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rff.analysis.scoring import decode_pairs  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump", required=True)
    ap.add_argument("--cut", type=float, default=0.5)
    ap.add_argument("--out", required=True)
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()

    src = ROOT / args.dump
    out = ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with src.open("r", encoding="utf-8") as fin, out.open("w", encoding="utf-8") as fout:
        for line in fin:
            if not line.strip():
                continue
            row = json.loads(line)
            scored = {
                tuple(int(x) for x in key.split(",")): float(value)
                for key, value in row["scored"].items()
            }
            preds = decode_pairs(
                scored, args.cut, int(row["length"]), partner=True, nesting=True
            )
            fout.write(
                json.dumps({"id": row["id"], "pairs": [list(p) for p in sorted(preds)]})
                + "\n"
            )
            n += 1
            if args.limit and n >= args.limit:
                break
    print(json.dumps({"out": args.out, "rows": n, "cut": args.cut}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
