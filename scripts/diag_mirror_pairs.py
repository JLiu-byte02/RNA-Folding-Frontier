#!/usr/bin/env python3
"""Diagnose why mirror labels differ from bpseq labels for the same sequence."""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rff.data.io import load_path  # noqa: E402


def main() -> int:
    pool: dict[str, list] = {}
    for line in open(ROOT / "data/processed/bprna/pool.jsonl"):
        rec = json.loads(line)
        pool.setdefault(rec["seq"], []).append((rec["id"], rec["pairs"]))

    recs = load_path(ROOT / "data/raw/bpfold/BPfold_data/bpRNA/TR0")
    shown = 0
    for rec in recs:
        cands = pool.get(rec.seq)
        if not cands:
            continue
        rid, pairs = cands[0]
        a = {tuple(sorted(p)) for p in rec.pairs}
        b = {tuple(sorted(p)) for p in pairs}
        if a == b:
            continue
        print("ID", rec.rid, "|", rid, "| len", len(rec.seq))
        print("  pairs bpseq/mirror/inter:", len(a), len(b), len(a & b))
        print("  bpseq-only:", sorted(a - b)[:8])
        print("  mirror-only:", sorted(b - a)[:8])
        print(
            "  bases:",
            [(i, j, rec.seq[i], rec.seq[j]) for i, j in sorted(b - a)[:6]],
        )
        # offset hypothesis: shift mirror indices by +1
        shifted = {(i - 1, j - 1) for i, j in b if i > 0 and j > 0}
        print("  match if mirror is 1-based:", len(a & shifted))
        shown += 1
        if shown >= 3:
            break
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
