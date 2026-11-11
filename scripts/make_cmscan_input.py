#!/usr/bin/env python3
"""Assemble the sequences that need an Rfam family assignment.

Rfam SEED sequences do not need cmscan: their family is known from the seed file
they came from. That leaves the bpRNA training pool, ArchiveII, bpRNA-new and
bpRNA-TS0, roughly 20k sequences instead of 170k.

Headers carry the source so the parser can attribute every assignment:
    >bprna:<id>            training pool
    >archiveii:<id>        ArchiveII test
    >bprnanew:<id>         bpRNA-new test
    >ts0:<id>              bpRNA TS0 test
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rff.data.io import load_path  # noqa: E402


def write_records(fh, tag: str, records) -> int:
    n = 0
    for rec in records:
        if not rec.seq or set(rec.seq) - set("ACGU"):
            continue
        fh.write(f">{tag}:{rec.rid}\n{rec.seq}\n")
        n += 1
    return n


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bprna-fasta", default="data/processed/bprna/cluster90/pool90")
    ap.add_argument("--out", default="data/processed/cmscan/input.fasta")
    ap.add_argument("--with-tests", action="store_true", default=True)
    args = ap.parse_args()

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    counts: dict[str, int] = {}

    with open(out, "w") as fh:
        n = 0
        name = None
        seq_parts: list[str] = []
        with open(ROOT / args.bprna_fasta) as src:
            for line in src:
                if line.startswith(">"):
                    if name is not None:
                        fh.write(f">{name}\n{''.join(seq_parts)}\n")
                        n += 1
                    rid = line[1:].strip().split()[0]
                    name = f"bprna:{rid}"
                    seq_parts = []
                else:
                    seq_parts.append(line.strip())
        if name is not None:
            fh.write(f">{name}\n{''.join(seq_parts)}\n")
            n += 1
        counts["bprna"] = n

        if args.with_tests:
            for tag, path in (
                ("archiveii", "data/raw/bpfold/BPfold_data/archiveII/archiveII"),
                ("ts0", "data/raw/bpfold/BPfold_data/bpRNA/TS0"),
                ("bprnanew", "data/raw/bpfold/BPfold_data/bpRNAnew/bpRNAnew.nr500.canonicals"),
            ):
                records = load_path(ROOT / path)
                counts[tag] = write_records(fh, tag, records)

    meta = {"out": str(args.out), "counts": counts, "total": sum(counts.values())}
    (out.parent / "input_manifest.json").write_text(json.dumps(meta, indent=2))
    print(json.dumps(meta, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
