#!/usr/bin/env python3
"""Build all-N and dinucleotide-shuffled inputs for template controls."""

from __future__ import annotations

import argparse
import json
import random
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def dinucleotide_shuffle(seq: str, rng: random.Random) -> str:
    """Altschul-Erickson shuffle preserving dinucleotide counts.

    The original sequence is an Euler trail in the directed graph whose nodes
    are bases and whose edges are consecutive dinucleotides. Randomising the
    adjacency order and running Hierholzer's algorithm samples another trail
    with the same edge multiset, hence the same dinucleotide counts and the same
    first/last base.
    """
    if len(seq) < 2:
        return seq
    adjacency: dict[str, list[str]] = defaultdict(list)
    for a, b in zip(seq, seq[1:]):
        adjacency[a].append(b)
    for edges in adjacency.values():
        rng.shuffle(edges)

    stack = [seq[0]]
    path: list[str] = []
    while stack:
        node = stack[-1]
        if adjacency.get(node):
            stack.append(adjacency[node].pop())
        else:
            path.append(stack.pop())
    path.reverse()
    if len(path) != len(seq) or path[0] != seq[0] or path[-1] != seq[-1]:
        raise ValueError("dinucleotide Euler trail reconstruction failed")
    return "".join(path)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--mode", choices=("all_n", "shuffle"), required=True)
    ap.add_argument("--seed", type=int, default=17)
    args = ap.parse_args()

    rng = random.Random(args.seed)
    records = []
    with open(ROOT / args.input, encoding="utf-8") as fh:
        for line in fh:
            if not line.strip():
                continue
            rec = json.loads(line)
            seq = rec["seq"]
            new_seq = "N" * len(seq) if args.mode == "all_n" else dinucleotide_shuffle(seq, rng)
            if len(new_seq) != len(seq):
                raise ValueError(f"{rec['id']}: length changed")
            out = dict(rec)
            out["seq"] = new_seq
            out["control"] = args.mode
            records.append(out)

    dest = ROOT / args.out
    dest.parent.mkdir(parents=True, exist_ok=True)
    with open(dest, "w", encoding="utf-8") as fh:
        for rec in records:
            fh.write(json.dumps(rec) + "\n")
    print(json.dumps({"out": args.out, "n": len(records), "mode": args.mode, "seed": args.seed}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
