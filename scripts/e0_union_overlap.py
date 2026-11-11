#!/usr/bin/env python3
"""ArchiveII vs the union of its two standard training corpora.

BPfold-class models train on RNAStrAlign and bpRNA-1m. Auditing each trainer
separately understates contamination, so this script reports the headline number
against the union, per ArchiveII family, plus the redundancy *inside* ArchiveII
that inflates per-sequence means.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rff.data.audit import near_duplicates, seq_hash, summarize  # noqa: E402
from rff.data.io import load_path  # noqa: E402


def load_archiveii(path: Path):
    records = load_path(path)
    for rec in records:
        rec.family = rec.rid.split("_", 1)[0]
    return records


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tr0", default="data/raw/bpfold/BPfold_data/bpRNA/TR0")
    ap.add_argument(
        "--rnastralign", default="data/raw/bpfold/BPfold_data/RNAStrAlign_bpseq"
    )
    ap.add_argument(
        "--archiveii", default="data/raw/bpfold/BPfold_data/archiveII/archiveII"
    )
    ap.add_argument("--out", required=True)
    ap.add_argument("--k", type=int, default=15)
    ap.add_argument("--min-containment", type=float, default=0.5)
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    tr0 = load_path(ROOT / args.tr0)
    rst = load_path(ROOT / args.rnastralign)
    arch = load_archiveii(ROOT / args.archiveii)
    union = tr0 + rst
    print(f"TR0 {len(tr0)}  RNAStrAlign {len(rst)}  union {len(union)}  ArchiveII {len(arch)}")

    tr0_hashes = {seq_hash(r.seq) for r in tr0}
    rst_hashes = {seq_hash(r.seq) for r in rst}
    union_hashes = tr0_hashes | rst_hashes

    per_family: dict[str, Counter] = defaultdict(Counter)
    hit_rows = []
    for rec in arch:
        h = seq_hash(rec.seq)
        in_tr0 = h in tr0_hashes
        in_rst = h in rst_hashes
        fam = rec.family or "unknown"
        per_family[fam]["n"] += 1
        if in_tr0:
            per_family[fam]["exact_tr0"] += 1
        if in_rst:
            per_family[fam]["exact_rnastralign"] += 1
        if in_tr0 or in_rst:
            per_family[fam]["exact_union"] += 1
            hit_rows.append(
                {"id": rec.rid, "family": fam, "in_tr0": in_tr0, "in_rnastralign": in_rst}
            )

    exact_union = sum(1 for r in arch if seq_hash(r.seq) in union_hashes)

    print("computing near-duplicates vs union ...", flush=True)
    near = near_duplicates(union, arch, k=args.k, min_containment=args.min_containment)
    near_ids = {h["query"] for h in near["hits"]}
    for rec in arch:
        if rec.rid in near_ids:
            per_family[rec.family or "unknown"]["near_union"] += 1

    # redundancy inside ArchiveII
    within = Counter(seq_hash(r.seq) for r in arch)
    within_dupes = {h: c for h, c in within.items() if c > 1}

    result = {
        "n_tr0": len(tr0),
        "n_rnastralign": len(rst),
        "n_archiveii": len(arch),
        "exact_tr0": sum(1 for r in arch if seq_hash(r.seq) in tr0_hashes),
        "exact_rnastralign": sum(1 for r in arch if seq_hash(r.seq) in rst_hashes),
        "exact_union": exact_union,
        "exact_union_rate": round(exact_union / len(arch), 4),
        "near_union_k": args.k,
        "near_union_min_containment": args.min_containment,
        "near_union_pairs": near["n_hits"],
        "near_union": len(near_ids),
        "near_union_rate": round(len(near_ids) / len(arch), 4),
        "within_archiveii_duplicate_sequences": len(within_dupes),
        "within_archiveii_extra_copies": sum(within_dupes.values()) - len(within_dupes),
        "per_family": {fam: dict(c) for fam, c in sorted(per_family.items())},
        "union_summary": summarize(arch),
    }
    (out / "summary.json").write_text(json.dumps(result, indent=2))
    (out / "exact_hits.json").write_text(json.dumps(hit_rows, indent=1))

    lines = [
        "# ArchiveII vs RNAStrAlign + bpRNA TR0",
        "",
        f"- ArchiveII sequences: {result['n_archiveii']}",
        f"- Exact duplicates of TR0: {result['exact_tr0']} "
        f"({result['exact_tr0'] / len(arch):.2%})",
        f"- Exact duplicates of RNAStrAlign: {result['exact_rnastralign']} "
        f"({result['exact_rnastralign'] / len(arch):.2%})",
        f"- **Exact duplicates of the union: {exact_union} "
        f"({result['exact_union_rate']:.2%})**",
        f"- Near-duplicates of the union (k={args.k}, containment >= "
        f"{args.min_containment}): {result['near_union']} "
        f"({result['near_union_rate']:.2%})",
        f"- Duplicate sequences inside ArchiveII: "
        f"{result['within_archiveii_duplicate_sequences']} "
        f"({result['within_archiveii_extra_copies']} extra copies)",
        "",
        "| family | n | exact TR0 | exact RNAStrAlign | exact union | near union |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for fam, c in sorted(per_family.items()):
        lines.append(
            f"| {fam} | {c['n']} | {c['exact_tr0']} | {c['exact_rnastralign']} | "
            f"{c['exact_union']} | {c['near_union']} |"
        )
    (out / "report.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
