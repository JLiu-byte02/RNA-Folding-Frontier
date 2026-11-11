#!/usr/bin/env python3
"""Extract per-family structure provenance from the Rfam SEED alignments.

Scoring a predicted structure against a consensus that was itself produced by a
folding tool is partly circular and favours thermodynamic methods, so the
test_novel primary table is restricted to families whose consensus structure is
published. This reads the Stockholm ``#=GF SS`` (structure source) and ``#=GF TP``
(RNA type) lines, counts base pairs in ``#=GF SS_cons``, and joins the result to
the frozen S1' split so the size of that primary subset is known before any model
sees the test data.
"""

from __future__ import annotations

import argparse
import gzip
import json
from collections import Counter, defaultdict
from pathlib import Path

PAIR_CHARS = {"(": ")", "[": "]", "{": "}", "<": ">"}


def count_consensus_pairs(ss_cons: str) -> int:
    stack: list[str] = []
    closes = {v: k for k, v in PAIR_CHARS.items()}
    n = 0
    for ch in ss_cons:
        if ch in PAIR_CHARS:
            stack.append(ch)
        elif ch in closes:
            want = closes[ch]
            if stack and stack[-1] == want:
                stack.pop()
                n += 1
    return n


def read_seed(path: Path) -> dict[str, dict]:
    fams: dict[str, dict] = {}
    cur: dict | None = None
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="latin-1") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("#=GF AC"):
                cur = {"accession": line.split(None, 2)[2].strip(), "rows": 0}
                fams[cur["accession"]] = cur
            elif line.startswith("#=GF ID") and cur is not None:
                cur["name"] = line.split(None, 2)[2].strip()
            elif line.startswith("#=GF TP") and cur is not None:
                cur["type"] = line.split(None, 2)[2].strip()
            elif line.startswith("#=GC SS_cons") and cur is not None:
                cur["ss_cons"] = line.split(None, 2)[2].strip()
            elif line.startswith("#=GF SS") and cur is not None:
                cur["ss_source"] = line.split(None, 2)[2].strip()
            elif line.startswith("#=GF BM") and cur is not None:
                cur.setdefault("build_method", line.split(None, 2)[2].strip())
            elif line and not line.startswith("#") and cur is not None:
                cur["rows"] = cur.get("rows", 0) + 1
            elif line == "//":
                cur = None
    return fams


def classify_source(raw: str) -> str:
    low = (raw or "").lower()
    if low.startswith("published") or "pmid" in low:
        return "published"
    if low.startswith("predicted") or "infernal" in low or low.startswith("consensus"):
        return "predicted"
    if not low:
        return "absent"
    return "other"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", default="data/raw/rfam/Rfam.seed.gz")
    ap.add_argument("--families", default="data/processed/rfam/families.tsv")
    ap.add_argument("--filtered-split", default="splits/s1prime_filtered")
    ap.add_argument("--out", default="experiments/001_e1prime/novel_label_metadata.json")
    ap.add_argument("--min-pairs", type=int, default=4)
    ap.add_argument("--cap-per-family", type=int, default=10)
    args = ap.parse_args()

    fams = read_seed(Path(args.seed))
    family_meta: dict[str, dict] = {}
    for line in Path(args.families).read_text().splitlines()[1:]:
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        acc, name, clan = parts[0], parts[1], parts[2]
        meta = fams.get(acc, {})
        family_meta[acc] = {
            "name": name,
            "clan": clan or acc,
            "type": meta.get("type", ""),
            "ss_source_raw": meta.get("ss_source", ""),
            "ss_source": classify_source(meta.get("ss_source", "")),
            "consensus_pairs": count_consensus_pairs(meta.get("ss_cons", "")),
            "seed_rows": meta.get("rows", 0),
        }

    by_clan_type: Counter = Counter()
    for meta in family_meta.values():
        by_clan_type[(meta["ss_source"], meta["type"])] += 1

    novel = [json.loads(l) for l in open(Path(args.filtered_split) / "test_novel.jsonl")]
    per_family: dict[str, list[dict]] = defaultdict(list)
    for rec in novel:
        per_family[rec["family"]].append(rec)

    primary_fams, eligible_fams = [], []
    primary_seqs = 0
    for fam, recs in per_family.items():
        meta = family_meta.get(fam)
        if meta is None:
            continue
        if meta["consensus_pairs"] >= args.min_pairs:
            eligible_fams.append(fam)
            if meta["ss_source"] == "published":
                primary_fams.append(fam)
                primary_seqs += min(len(recs), args.cap_per_family)

    type_counts = Counter(
        family_meta[f]["type"] for f in primary_fams if f in family_meta
    )
    report = {
        "families_parsed": len(fams),
        "ss_source_counts_all_families": dict(
            Counter(m["ss_source"] for m in family_meta.values())
        ),
        "ss_source_by_type_all_families": {
            f"{k[0]}|{k[1]}": v for k, v in by_clan_type.most_common(40)
        },
        "test_novel": {
            "sequences": len(novel),
            "families": len(per_family),
            "eligible_families_min_pairs": len(eligible_fams),
            "primary_published_families": len(primary_fams),
            "primary_sequence_budget_capped": primary_seqs,
            "primary_type_counts": dict(type_counts.most_common()),
        },
        "rule": {
            "min_pairs_in_ss_cons": args.min_pairs,
            "cap_per_family": args.cap_per_family,
            "requires": "ss_source == published",
        },
        "below_150_clans": len(primary_fams) < 150,
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2))
    (out.parent / "seed_family_metadata.json").write_text(
        json.dumps(family_meta, indent=0)
    )
    print(json.dumps(report["test_novel"], indent=2))
    print(json.dumps({"below_150_clans": report["below_150_clans"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
