#!/usr/bin/env python3
"""Find a low-identity same-family homolog for each training sequence.

This builds condition (c) of E1': replace a sequence's distal context with the
corresponding region of a *homolog* that is still the same family but below 60%
identity. Same-family means the fold is informative, low identity means the
model cannot copy the answer from a near-duplicate.

Procedure per family: cluster at 60% identity with CD-HIT, pick a candidate from
a different cluster, then verify the candidate with a global alignment
(Biopython) and keep it only if identity falls in [min, max). The verified
identity is stored so the condition can be described exactly rather than
assumed.
"""

from __future__ import annotations

import argparse
import json
import random
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

from Bio import Align
from Bio.Align import substitution_matrices


def read_jsonl(path: Path) -> list[dict]:
    with open(path) as fh:
        return [json.loads(line) for line in fh if line.strip()]


def run_cdhit(
    seqs: dict[str, str], outdir: Path, identity: float, cdhit: str
) -> dict[str, str]:
    """id -> cluster representative id, from one CD-HIT run."""
    fa = outdir / "in.fa"
    with open(fa, "w") as fh:
        for rid, seq in seqs.items():
            fh.write(f">{rid}\n{seq}\n")
    prefix = outdir / "clust"
    cmd = [
        cdhit,
        "-i",
        str(fa),
        "-o",
        str(prefix),
        "-c",
        str(identity),
        "-n",
        "5",
        "-d",
        "0",
        "-T",
        "2",
        "-M",
        "0",
    ]
    subprocess.run(cmd, check=True, capture_output=True)
    clusters: dict[str, str] = {}
    current: str | None = None
    with open(str(prefix) + ".clstr") as fh:
        for line in fh:
            if line.startswith(">Cluster"):
                current = line.strip().split()[-1]
                continue
            if current is None:
                continue
            rid = line.split(">")[1].split("...")[0].strip()
            clusters[rid] = current
    return clusters


def identity(aligner: Align.PairwiseAligner, a: str, b: str) -> float:
    aln = aligner.align(a, b)[0]
    matches = 0
    for (s1, e1), (s2, e2) in zip(aln.aligned[0], aln.aligned[1]):
        matches += sum(x == y for x, y in zip(a[s1:e1], b[s2:e2]))
    return matches / max(1, min(len(a), len(b)))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default="splits/s1prime_filtered/train.jsonl")
    ap.add_argument("--out", default="data/processed/homologs/homolog_map.json")
    ap.add_argument("--identity", type=float, default=0.60)
    ap.add_argument("--min-identity", type=float, default=0.30)
    ap.add_argument("--min-cluster-size", type=int, default=2)
    ap.add_argument("--max-verify", type=int, default=6)
    ap.add_argument("--seed", type=int, default=17)
    ap.add_argument("--tmp", default="data/processed/homologs/tmp")
    ap.add_argument(
        "--cdhit",
        default="tools/cdhit-4.8.1/cd-hit-est",
        help="path to the cd-hit-est binary built from source",
    )
    args = ap.parse_args()

    rng = random.Random(args.seed)
    recs = read_jsonl(Path(args.input))
    by_family: dict[str, list[dict]] = defaultdict(list)
    for rec in recs:
        by_family[rec["family"]].append(rec)

    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("NUC.4.4")
    aligner.open_gap_score = -6.0
    aligner.extend_gap_score = -1.0
    aligner.mode = "global"

    tmp_root = Path(args.tmp)
    mapping: dict[str, dict] = {}
    stats = {"families": 0, "with_candidate": 0, "families_skipped": 0}
    for fam, group in sorted(by_family.items()):
        if len(group) < 2 * args.min_cluster_size:
            stats["families_skipped"] += 1
            continue
        stats["families"] += 1
        outdir = tmp_root / fam.replace("/", "_")
        outdir.mkdir(parents=True, exist_ok=True)
        clusters = run_cdhit(
            {r["id"]: r["seq"] for r in group}, outdir, args.identity, args.cdhit
        )
        members: dict[str, list[dict]] = defaultdict(list)
        for rec in group:
            c = clusters.get(rec["id"])
            if c is not None:
                members[c].append(rec)
        if len(members) < 2:
            continue
        cluster_ids = sorted(members)
        for rec in group:
            own = clusters.get(rec["id"])
            pool = [c for c in cluster_ids if c != own]
            if not pool:
                continue
            rng.shuffle(pool)
            chosen = None
            for cand_cluster in pool[: args.max_verify]:
                cand = rng.choice(members[cand_cluster])
                if cand["id"] == rec["id"]:
                    continue
                ident = identity(aligner, rec["seq"], cand["seq"])
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

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(mapping, indent=1))
    report = {
        "input": args.input,
        "n_records": len(recs),
        "n_mapped": len(mapping),
        "coverage": round(len(mapping) / max(1, len(recs)), 4),
        "identity_range": [args.min_identity, args.identity],
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
