#!/usr/bin/env python3
"""Split a FASTA into fixed-size shards, keeping records intact."""

from __future__ import annotations

import argparse
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("input")
    ap.add_argument("outdir")
    ap.add_argument("--size", type=int, default=500)
    args = ap.parse_args()

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    for old in out.glob("shard_*.fa"):
        old.unlink()

    records: list[tuple[str, str]] = []
    name, buf = None, []
    with open(args.input) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(buf)))
                name, buf = line[1:].strip().split()[0], []
            else:
                buf.append(line.strip())
    if name is not None:
        records.append((name, "".join(buf)))

    n_shards = 0
    for start in range(0, len(records), args.size):
        chunk = records[start : start + args.size]
        path = out / f"shard_{start // args.size:05d}.fa"
        with open(path, "w") as fh:
            for name, seq in chunk:
                fh.write(f">{name}\n{seq}\n")
        n_shards += 1
    print(f"records: {len(records)}, shards: {n_shards}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
