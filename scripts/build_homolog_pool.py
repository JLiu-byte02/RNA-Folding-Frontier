#!/usr/bin/env python3
"""Find a low-identity same-family homolog for each training sequence.

This builds condition (c) of E1': replace a sequence's distal context with the
corresponding region of a *homolog* that is still the same family but below 60%
identity. Same-family means the fold is informative, low identity means the
model cannot copy the answer from a near-duplicate.

Procedure per family: sample same-family candidates (length-restricted, so the
alignment stays cheap), align each globally, and keep the first whose identity
falls in [min, max). Identity is the fraction of alignment columns that are
identical matches, the usual definition; the verified value is stored so the
condition can be described rather than assumed.

An earlier version clustered with CD-HIT at 0.60 and failed: CD-HIT 4.8.1
rejects nucleotide thresholds below 0.8 ("invalid clstr threshold"). Clustering
is therefore not used; the alignment band is the filter, which is the property
the condition actually needs.
"""

from __future__ import annotations

import argparse
import json
import random
from collections import defaultdict
from pathlib import Path

from Bio import Align


def read_jsonl(path: Path) -> list[dict]:
    with open(path) as fh:
        return [json.loads(line) for line in fh if line.strip()]


def identity(aligner: Align.PairwiseAligner, a: str, b: str) -> float:
    """Identical columns / alignment columns."""
    aln = aligner.align(a, b)[0]
    matches = columns = 0
    for (s1, e1), (s2, e2) in zip(aln.aligned[0], aln.aligned[1]):
        matches += sum(x == y for x, y in zip(a[s1:e1], b[s2:e2]))
        columns += e1 - s1
    # Ungapped blocks are the aligned columns; gaps only lengthen the alignment,
    # so using aligned columns keeps the denominator the matched span and makes
    # the threshold slightly stricter than a gapped-column count.
    return matches / max(1, columns)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default="splits/s1prime_filtered/train.jsonl")
    ap.add_argument("--out", default="data/processed/homologs/homolog_map.json")
    ap.add_argument("--identity", type=float, default=0.60)
    ap.add_argument("--min-identity", type=float, default=0.40)
    ap.add_argument("--max-verify", type=int, default=12)
    ap.add_argument("--len-ratio", type=float, default=0.7)
    ap.add_argument("--seed", type=int, default=17)
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()

    rng = random.Random(args.seed)
    recs = read_jsonl(Path(args.input))
    by_family: dict[str, list[dict]] = defaultdict(list)
    for rec in recs:
        by_family[rec["family"]].append(rec)

    aligner = Align.PairwiseAligner()
    # Scoring by match/mismatch rather than a substitution matrix: bpRNA and
    # Rfam seed rows contain IUPAC ambiguity codes that NUC.4.4 rejects.
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -6.0
    aligner.extend_gap_score = -1.0
    aligner.mode = "global"

    mapping: dict[str, dict] = {}
    stats = {"families": 0, "with_candidate": 0, "families_skipped": 0, "alignments": 0}
    for fam, group in sorted(by_family.items()):
        if len(group) < 2:
            stats["families_skipped"] += 1
            continue
        stats["families"] += 1
        for rec in group:
            L = len(rec["seq"])
            pool = [
                other
                for other in group
                if other["id"] != rec["id"]
                and args.len_ratio * L <= len(other["seq"]) <= L / args.len_ratio
            ]
            if not pool:
                continue
            rng.shuffle(pool)
            chosen = None
            for cand in pool[: args.max_verify]:
                ident = identity(aligner, rec["seq"], cand["seq"])
                stats["alignments"] += 1
                if args.min_identity <= ident < args.identity:
                    chosen = {**cand, "identity": round(float(ident), 3)}
                    break
            if chosen:
                mapping[rec["id"]] = {
                    "homolog": chosen["id"],
                    "identity": chosen["identity"],
                    "homolog_len": len(chosen["seq"]),
                    "family": fam,
                }
                stats["with_candidate"] += 1
        if args.limit and stats["with_candidate"] >= args.limit:
            break

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(mapping, indent=1))
    report = {
        "input": args.input,
        "n_records": len(recs),
        "n_mapped": len(mapping),
        "coverage": round(len(mapping) / max(1, len(recs)), 4),
        "identity_range": [args.min_identity, args.identity],
        "len_ratio": args.len_ratio,
        "max_verify": args.max_verify,
        "per_family": {
            fam: {
                "n": len(group),
                "mapped": sum(1 for r in group if r["id"] in mapping),
            }
            for fam, group in sorted(by_family.items())
        },
        **stats,
    }
    (out.parent / "homolog_report.json").write_text(json.dumps(report, indent=2))
    print(json.dumps({k: v for k, v in report.items() if k != "per_family"}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
