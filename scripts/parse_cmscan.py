#!/usr/bin/env python3
"""Collapse cmscan tblout shards into one family-assignment table.

Output columns: source, query, family_accession, clan, evalue, score, bias.
Sequences with no hit above the reporting threshold get family ``NOHIT``.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def load_clans(path: Path) -> dict[str, str]:
    clan_of: dict[str, str] = {}
    if not path.exists():
        return clan_of
    for line in path.read_text().splitlines():
        parts = line.split("\t")
        if len(parts) >= 2 and parts[0].startswith("CL"):
            clan_of[parts[1]] = parts[0]
    return clan_of


def load_names(path: Path) -> dict[str, str]:
    names: dict[str, str] = {}
    if not path.exists():
        return names
    for line in path.read_text().splitlines():
        parts = line.split("\t")
        if len(parts) >= 2:
            names[parts[0]] = parts[1]
    return names


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tbl-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--clan-membership", default="data/raw/rfam/clan_membership.txt")
    ap.add_argument("--family-table", default="data/processed/rfam/families.tsv")
    args = ap.parse_args()

    clan_of = load_clans(Path(args.clan_membership))
    names = load_names(Path(args.family_table))

    best: dict[str, dict] = {}
    seen: set[str] = set()
    n_lines = 0
    for tbl in sorted(Path(args.tbl_dir).glob("*.tbl")):
        with open(tbl) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                n_lines += 1
                fields = line.split()
                if len(fields) < 18:
                    continue
                target_name = fields[0]
                query_name = fields[2]
                accession = fields[1]
                try:
                    evalue = float(fields[16])
                    score = float(fields[15])
                    bias = float(fields[17])
                except ValueError:
                    continue
                seen.add(query_name)
                cur = best.get(query_name)
                if cur is None or score > cur["score"]:
                    best[query_name] = {
                        "query": query_name,
                        "target_name": target_name,
                        "family_accession": accession,
                        "family_name": names.get(accession, target_name),
                        "clan": clan_of.get(accession, "SINGLETON"),
                        "evalue": evalue,
                        "score": score,
                        "bias": bias,
                    }

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    cols = [
        "source",
        "query",
        "family_accession",
        "family_name",
        "clan",
        "evalue",
        "score",
        "bias",
    ]
    with open(out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        writer.writeheader()
        for query, row in sorted(best.items()):
            source = query.split(":", 1)[0] if ":" in query else "unknown"
            writer.writerow({**row, "source": source})

    summary = {
        "tbl_lines": n_lines,
        "queries_with_hit": len(seen),
        "queries_written": len(best),
        "by_source": {},
    }
    for row in best.values():
        source = row["query"].split(":", 1)[0] if ":" in row["query"] else "unknown"
        summary["by_source"].setdefault(source, {"hit": 0, "families": 0})
        summary["by_source"][source]["hit"] += 1
    (out.parent / "parse_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
