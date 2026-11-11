#!/usr/bin/env python3
"""Freeze the test_novel evaluation subset (PREREGISTRATION Amendment 1.6).

Primary set: novel clans whose SEED ``#=GF SS`` source is *published* and whose
SS_cons has at least four base pairs.  Within each family, sequences are
deduplicated at 90% identity and capped at ten, drawn in a seeded random order,
so the cap cannot be gamed by taking the easiest sequences first.

The output is written with a sha256 over both the sequence file and the
catalogue, so "no model is scored on test_novel until this is committed" is
checkable by anyone: the hash in the repo must match the file the scores were
computed from.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import random
from collections import defaultdict
from pathlib import Path

from Bio.Align import PairwiseAligner

ROOT = Path(__file__).resolve().parents[1]


def read_jsonl(path: Path) -> list[dict]:
    with open(path) as fh:
        return [json.loads(line) for line in fh if line.strip()]


def type_bucket(rfam_type: str) -> str:
    t = rfam_type.lower()
    if "cis-reg" in t:
        return "cis-regulatory/riboswitch"
    if "ribozyme" in t:
        return "ribozyme"
    if "snorna" in t:
        return "snoRNA"
    if "mirna" in t:
        return "miRNA"
    if "srna" in t:
        return "sRNA"
    return "other"


def has_multiloop(pairs: list[tuple[int, int]]) -> bool:
    """True if any loop is closed by a pair that directly contains >=2 pairs.

    A multiloop is a loop with at least two directly enclosed helices; counting
    direct children (not all descendants) is what distinguishes it from a
    simple bulged stem.
    """
    if not pairs:
        return False
    partner: dict[int, int] = {}
    for i, j in pairs:
        partner[min(i, j)] = max(i, j)
        partner[max(i, j)] = min(i, j)
    open_stack: list[int] = []
    children: dict[int, int] = defaultdict(int)
    for pos in range(max(partner) + 1):
        q = partner.get(pos)
        if q is None:
            continue
        if q > pos:
            if open_stack:
                children[open_stack[-1]] += 1
            open_stack.append(pos)
        elif open_stack and open_stack[-1] == q:
            open_stack.pop()
    return any(v >= 2 for v in children.values())


def identity(aligner: PairwiseAligner, a: str, b: str) -> float:
    if not a or not b:
        return 0.0
    aln = aligner.align(a, b)[0]
    matches = sum(
        1
        for (x, y) in zip(aln[0], aln[1])
        if x == y and x != "-"
    )
    return matches / max(len(a), len(b))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--split", default="splits/s1prime_filtered/test_novel.jsonl")
    ap.add_argument("--metadata", default="experiments/001_e1prime/seed_family_metadata.json")
    ap.add_argument("--out-jsonl", default="splits/test_novel_eval.jsonl")
    ap.add_argument("--out-catalogue", default="experiments/001_e1prime/test_novel_catalogue.json")
    ap.add_argument("--min-pairs", type=int, default=4)
    ap.add_argument("--max-per-family", type=int, default=10)
    ap.add_argument("--max-identity", type=float, default=0.90)
    ap.add_argument("--seed", type=int, default=17)
    ap.add_argument("--min-clans", type=int, default=150)
    args = ap.parse_args()

    meta = json.loads((ROOT / args.metadata).read_text())
    recs = read_jsonl(ROOT / args.split)
    by_family: dict[str, list[dict]] = defaultdict(list)
    for rec in recs:
        by_family[rec["family"]].append(rec)

    eligible = [
        f
        for f in sorted(by_family)
        if meta.get(f, {}).get("ss_source") == "published"
        and meta.get(f, {}).get("consensus_pairs", 0) >= args.min_pairs
    ]

    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -2.0
    aligner.extend_gap_score = -0.5

    rng = random.Random(args.seed)
    selected: list[dict] = []
    families = {}
    for fam in eligible:
        pool = sorted(by_family[fam], key=lambda r: r["id"])
        rng.shuffle(pool)
        kept: list[dict] = []
        for rec in pool:
            if len(kept) >= args.max_per_family:
                break
            if all(
                identity(aligner, rec["seq"], other["seq"]) < args.max_identity
                for other in kept
            ):
                kept.append(rec)
        selected.extend(kept)
        m = meta[fam]
        families[fam] = {
            "clan": m["clan"],
            "name": m["name"],
            "type_raw": m["type"],
            "type": type_bucket(m["type"]),
            "consensus_pairs": m["consensus_pairs"],
            "ss_source_raw": m["ss_source_raw"],
            "n_in_split": len(pool),
            "n_selected": len(kept),
            "multiloop": has_multiloop(
                [tuple(p) for p in kept[0]["pairs"]] if kept else []
            ),
        }
        print(f"{fam} {m['name']}: {len(kept)}/{len(pool)}", flush=True)

    clans = {v["clan"] for v in families.values()}
    out_jsonl = ROOT / args.out_jsonl
    out_jsonl.parent.mkdir(parents=True, exist_ok=True)
    with open(out_jsonl, "w") as fh:
        for rec in selected:
            fh.write(json.dumps(rec) + "\n")
    seq_hash = hashlib.sha256(out_jsonl.read_bytes()).hexdigest()

    catalogue = {
        "protocol": "PREREGISTRATION.md Amendment 1.6",
        "source_split": str(args.split),
        "params": {
            k: v
            for k, v in vars(args).items()
            if k not in {"out_jsonl", "out_catalogue"}
        },
        "n_families": len(families),
        "n_clans": len(clans),
        "n_sequences": len(selected),
        "min_clans_gate": args.min_clans,
        "gate_passed": len(clans) >= args.min_clans,
        "sequences_sha256": seq_hash,
        "families": families,
    }
    cat_path = ROOT / args.out_catalogue
    cat_path.parent.mkdir(parents=True, exist_ok=True)
    cat_path.write_text(json.dumps(catalogue, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                "families": len(families),
                "clans": len(clans),
                "sequences": len(selected),
                "gate_passed": catalogue["gate_passed"],
                "sha256": seq_hash,
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
