#!/usr/bin/env python3
"""The two audit counts Amendment 1.5 asks for, from the -E 10 rescan.

(i)  How many NOHIT pool sequences hit a *test* clan at E <= 10.
(ii) How many training sequences that already have a GA-passing hit also carry
     a secondary hit to a test clan at E <= 10, and what fraction of training
     that is.  The escalation threshold is 0.5% of training.

"Test clan" means a clan deliberately held out of training: the withheld clans or
the S1'-novel clans.  The frozen scan used --cut_ga, so by construction it should
report almost none of these inside training; this measures how much of that is
the cutoff's doing rather than the sequence's.
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def read_tsv(path: Path) -> list[dict]:
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(header, line.rstrip("\n").split("\t"))) for line in fh if line.strip()]


def read_jsonl(path: Path) -> list[dict]:
    with open(path) as fh:
        return [json.loads(line) for line in fh if line.strip()]


def parse_relaxed(tbl_dir: Path, min_hits_only: bool = True) -> dict[str, list[tuple[str, str, float]]]:
    """query -> [(family_accession, clan_label, evalue)]. `!` = first hit only."""
    hits: dict[str, list[tuple[str, str, float]]] = defaultdict(list)
    for tbl in sorted(tbl_dir.glob("*.tbl")):
        with open(tbl) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                f = line.split()
                if len(f) < 16:
                    continue
                # columns: target, acc, query, qacc, mdl, mfrom, mto, sfrom, sto,
                #          strand, trunc, pass, gc, bias, score, evalue, inc, desc
                if min_hits_only and f[11] != "1":
                    continue
                query, acc, ev = f[2], f[1], float(f[15])
                hits[query].append((acc, "", ev))
    return hits


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--relaxed-tbl", default="data/processed/cmscan_relaxed/tbl")
    ap.add_argument("--assignments", default="data/processed/cmscan/assignments.tsv")
    ap.add_argument("--family-clan", default="data/processed/cmscan/family_clan.tsv")
    ap.add_argument("--splits", default="splits/s1prime_filtered")
    ap.add_argument("--pool", default="data/processed/bprna/pool.jsonl")
    ap.add_argument("--evalue", type=float, default=10.0)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    splits = Path(args.splits)
    train = read_jsonl(splits / "train.jsonl")
    withheld = read_jsonl(splits / "test_withheld.jsonl")
    novel = read_jsonl(splits / "test_novel.jsonl")
    test_clans = {r["clan"] for r in withheld} | {r["clan"] for r in novel}
    train_ids = {r["id"] for r in train}

    fam2clan: dict[str, str] = {}
    fc = ROOT / args.family_clan
    if fc.exists():
        for row in read_tsv(fc):
            fam2clan[row["family_accession"]] = row["clan"]
    ga_rows = read_tsv(ROOT / args.assignments)
    for row in ga_rows:
        fam2clan.setdefault(row["family_accession"], row["clan"])
    ga_hit = {row["query"] for row in ga_rows}

    relaxed = parse_relaxed(ROOT / args.relaxed_tbl)

    def test_clan_hits(query: str) -> set[str]:
        out = set()
        for acc, _clan, ev in relaxed.get(query, []):
            if ev > args.evalue:
                continue
            clan = fam2clan.get(acc)
            if clan and clan in test_clans:
                out.add(clan)
        return out

    pool_ids = [r["id"] for r in read_jsonl(ROOT / args.pool)] if (ROOT / args.pool).exists() else []
    nohit = [q for q in pool_ids if q not in ga_hit]
    nohit_hits = {q: test_clan_hits(q) for q in nohit}
    n_nohit_touch = sum(1 for v in nohit_hits.values() if v)

    train_touch = {}
    for q in sorted(train_ids):
        if q in ga_hit:
            t = test_clan_hits(q)
            if t:
                train_touch[q] = sorted(t)

    report = {
        "relaxed_tbl": str(args.relaxed_tbl),
        "evalue_threshold": args.evalue,
        "n_test_clans": len(test_clans),
        "n_pool": len(pool_ids),
        "n_pool_nohit": len(nohit),
        "n_training": len(train_ids),
        "n_training_with_ga_hit": len([q for q in train_ids if q in ga_hit]),
        "audit_i_nohit_touching_test_clan": n_nohit_touch,
        "audit_i_fraction_of_nohit": round(n_nohit_touch / max(1, len(nohit)), 4),
        "audit_ii_training_touching_test_clan": len(train_touch),
        "audit_ii_fraction_of_training": round(len(train_touch) / max(1, len(train_ids)), 5),
        "audit_ii_escalation_threshold": 0.005,
        "audit_ii_escalate": len(train_touch) / max(1, len(train_ids)) > 0.005,
        "examples": dict(list(train_touch.items())[:20]),
    }
    out = ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
