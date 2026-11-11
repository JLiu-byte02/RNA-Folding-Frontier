"""Leakage audit for RNA secondary-structure splits.

The audit is deliberately conservative and cheap enough to run before any model
training:

1. exact sequence duplicates,
2. near-duplicates via shared k-mers,
3. family / clan overlap between splits,
4. basic distribution summaries.
"""

from __future__ import annotations

import hashlib
import json
from collections import Counter, defaultdict
from dataclasses import asdict
from pathlib import Path
from typing import Iterable, Sequence

from .io import Record


def seq_hash(seq: str) -> str:
    return hashlib.sha256(seq.encode()).hexdigest()[:16]


def summarize(records: Sequence[Record]) -> dict:
    lengths = [len(r.seq) for r in records]
    pairs = [len(r.pairs) for r in records]
    families = Counter(r.family or "unknown" for r in records)
    return {
        "n": len(records),
        "length_min": min(lengths) if lengths else 0,
        "length_median": sorted(lengths)[len(lengths) // 2] if lengths else 0,
        "length_max": max(lengths) if lengths else 0,
        "pairs_total": sum(pairs),
        "families": len(families),
        "top_families": families.most_common(10),
    }


def exact_duplicates(splits: dict[str, Sequence[Record]]) -> dict:
    by_hash: dict[str, list[tuple[str, str]]] = defaultdict(list)
    for split, records in splits.items():
        for rec in records:
            by_hash[seq_hash(rec.seq)].append((split, rec.rid))
    cross = {
        h: members
        for h, members in by_hash.items()
        if len({m[0] for m in members}) > 1
    }
    return {
        "n_cross_split_duplicate_sequences": len(cross),
        "examples": [
            {"hash": h, "members": members[:6]} for h, members in list(cross.items())[:20]
        ],
    }


def exact_duplicates_vs_reference(
    reference: Sequence[Record], query: Sequence[Record]
) -> dict:
    ref_hashes = {seq_hash(r.seq) for r in reference}
    hits = [q.rid for q in query if seq_hash(q.seq) in ref_hashes]
    return {
        "n_query": len(query),
        "n_hits": len(hits),
        "hit_rate": round(len(hits) / len(query), 4) if query else 0.0,
        "examples": hits[:20],
    }


def _kmers(seq: str, k: int) -> set[str]:
    return {seq[i : i + k] for i in range(max(0, len(seq) - k + 1))}


def near_duplicates(
    reference: Sequence[Record],
    query: Sequence[Record],
    k: int = 15,
    min_containment: float = 0.5,
    max_pairs: int = 20000,
) -> dict:
    """Find query sequences largely contained in reference sequences.

    Containment is measured on k-mer sets. Candidate pairs are generated with an
    inverted index, then verified with exact Jaccard/containment.
    """
    ref_kmers = [(_kmers(r.seq, k), r) for r in reference]
    index: dict[str, list[int]] = defaultdict(list)
    for idx, (kmers, _) in enumerate(ref_kmers):
        for kmer in kmers:
            index[kmer].append(idx)

    hits: list[dict] = []
    query_with_hit: set[str] = set()
    for qrec in query:
        qset = _kmers(qrec.seq, k)
        if not qset:
            continue
        overlap: Counter[int] = Counter()
        for kmer in qset:
            for idx in index.get(kmer, ()):  # type: ignore[arg-type]
                overlap[idx] += 1
        for idx, shared in overlap.most_common(5):
            rset, rrec = ref_kmers[idx]
            containment = shared / min(len(qset), len(rset))
            if containment < min_containment:
                continue
            union = len(qset | rset)
            hits.append(
                {
                    "query": qrec.rid,
                    "reference": rrec.rid,
                    "query_len": qrec.length,
                    "reference_len": rrec.length,
                    "containment": round(containment, 3),
                    "jaccard": round(shared / union, 3) if union else 0.0,
                }
            )
            query_with_hit.add(qrec.rid)
            if len(hits) >= max_pairs:
                break
        if len(hits) >= max_pairs:
            break
    hits.sort(key=lambda h: h["jaccard"], reverse=True)
    return {
        "k": k,
        "min_containment": min_containment,
        "n_query": len(query),
        "n_query_with_hit": len(query_with_hit),
        "hit_rate": round(len(query_with_hit) / len(query), 4) if query else 0.0,
        "n_hits": len(hits),
        "hits": hits,
    }


def metadata_overlap(
    train: Sequence[Record], test: Sequence[Record], field: str
) -> dict:
    train_values = {getattr(r, field) for r in train if getattr(r, field)}
    test_values = {getattr(r, field) for r in test if getattr(r, field)}
    shared = train_values & test_values
    return {
        "field": field,
        "train_values": len(train_values),
        "test_values": len(test_values),
        "shared_values": len(shared),
        "shared_fraction_of_test": round(
            len(shared) / len(test_values), 4
        )
        if test_values
        else 0.0,
        "examples": sorted(shared)[:20],
    }


def run_audit(
    splits: dict[str, Sequence[Record]],
    reference_split: str,
    pairwise_k: int = 15,
) -> dict:
    train = splits[reference_split]
    report: dict = {
        "reference_split": reference_split,
        "summaries": {name: summarize(recs) for name, recs in splits.items()},
        "exact_duplicates": exact_duplicates(splits),
        "per_test": {},
    }
    for name, recs in splits.items():
        if name == reference_split:
            continue
        report["per_test"][name] = {
            "exact_duplicates_vs_reference": exact_duplicates_vs_reference(
                train, recs
            ),
            "near_duplicates_vs_reference": near_duplicates(
                train, recs, k=pairwise_k
            ),
            "family_overlap": metadata_overlap(train, recs, "family"),
            "clan_overlap": metadata_overlap(train, recs, "clan"),
        }
    return report


def to_markdown(report: dict) -> str:
    lines = ["# Leakage audit", ""]
    lines.append(f"Reference split: `{report['reference_split']}`")
    lines.append("")
    lines.append("| split | n | len median | len max | families |")
    lines.append("|---|---:|---:|---:|---:|")
    for name, summary in report["summaries"].items():
        lines.append(
            f"| {name} | {summary['n']} | {summary['length_median']} | "
            f"{summary['length_max']} | {summary['families']} |"
        )
    lines.append("")
    exact = report["exact_duplicates"]["n_cross_split_duplicate_sequences"]
    lines.append(f"Cross-split exact duplicate sequences: **{exact}**")
    lines.append("")
    for name, payload in report["per_test"].items():
        lines.append(f"## {name}")
        ex = payload["exact_duplicates_vs_reference"]
        lines.append(
            f"- exact duplicates vs reference: {ex['n_hits']} "
            f"({ex['hit_rate']:.2%} of test sequences)"
        )
        nd = payload["near_duplicates_vs_reference"]
        fam = payload["family_overlap"]
        clan = payload["clan_overlap"]
        lines.append(
            f"- near-duplicate query sequences (k={nd['k']}): "
            f"{nd['n_query_with_hit']} ({nd['hit_rate']:.2%}); "
            f"{nd['n_hits']} pairs"
        )
        lines.append(
            f"- shared families: {fam['shared_values']} "
            f"({fam['shared_fraction_of_test']:.1%} of test families)"
        )
        lines.append(
            f"- shared clans: {clan['shared_values']} "
            f"({clan['shared_fraction_of_test']:.1%} of test clans)"
        )
        lines.append("")
    return "\n".join(lines)


def save_report(report: dict, out_dir: str | Path) -> None:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    (out / "audit.json").write_text(json.dumps(report, indent=2))
    (out / "audit.md").write_text(to_markdown(report))

