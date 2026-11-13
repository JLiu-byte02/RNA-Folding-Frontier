#!/usr/bin/env python3
"""Build learned one-hot nucleotide caches for the E-C1b control."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import torch

ROOT = Path(__file__).resolve().parents[1]
INDEX = {"A": 0, "C": 1, "G": 2, "U": 3}


def one_hot(seq: str) -> torch.Tensor:
    out = torch.zeros((len(seq), 4), dtype=torch.float16)
    for i, char in enumerate(seq):
        idx = INDEX.get(char.upper())
        if idx is not None:
            out[i, idx] = 1.0
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--split", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--max-len", type=int, default=1024)
    ap.add_argument("--shard-size", type=int, default=256)
    args = ap.parse_args()

    records = [
        json.loads(line)
        for line in open(ROOT / args.split, encoding="utf-8")
        if line.strip()
    ]
    out = ROOT / args.out
    out.mkdir(parents=True, exist_ok=True)
    ids, embs, n_shard = [], [], 0

    def flush() -> None:
        nonlocal ids, embs, n_shard
        if not ids:
            return
        torch.save(
            {
                "ids": ids,
                "lengths": torch.tensor([e.shape[0] for e in embs]),
                "embs": embs,
            },
            out / f"shard_{n_shard:04d}.pt",
        )
        n_shard += 1
        ids, embs = [], []

    for rec in records:
        seq = rec["seq"][: args.max_len]
        ids.append(rec["id"])
        embs.append(one_hot(seq))
        if len(ids) >= args.shard_size:
            flush()
    flush()
    (out / "manifest.json").write_text(
        json.dumps(
            {
                "kind": "one-hot",
                "data": args.split,
                "records_cached": len(records),
                "shards": n_shard,
                "max_len": args.max_len,
                "dim": 4,
            },
            indent=2,
        )
        + "\n"
    )
    print(json.dumps({"out": args.out, "records": len(records), "shards": n_shard}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
