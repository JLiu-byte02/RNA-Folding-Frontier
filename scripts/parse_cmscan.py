#!/usr/bin/env python3
"""Collapse cmscan tblout shards into assignment tables.

Outputs in the directory of ``--out``:

``assignments.tsv``
    Best-scoring family per query (columns: source, query, target_name,
    family_accession, family_name, clan, evalue, score, bias).  Queries with no
    reported hit get ``NOHIT``.
``clan_hits.tsv``
    Best hit per (query, clan).  This is what clan-level leakage analysis must
    use: a sequence can hit several families of one clan and a single-label
    table would hide that.
``all_hits.tsv``
    Every reported hit, with alignment coordinates, for overlap/consistency
    checks downstream.

With ``--only-done`` the parser considers only shards that carry a ``.done``
marker and restricts NOHIT rows to the sequences covered by those shards, so a
partial run can be validated without inventing false negatives.
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
    ap.add_argument(
        "--input-fasta",
        default="data/processed/cmscan/input.fasta",
        help="used to emit explicit NOHIT rows",
    )
    ap.add_argument(
        "--shards-dir",
        default=None,
        help="defaults to <tbl-dir>/../shards; used to resolve covered ids",
    )
    ap.add_argument(
        "--only-done",
        action="store_true",
        help="restrict to shards with a .done marker (partial-run validation)",
    )
    args = ap.parse_args()

    clan_of = load_clans(Path(args.clan_membership))
    names = load_names(Path(args.family_table))

    tbl_dir = Path(args.tbl_dir)
    shards_dir = Path(args.shards_dir) if args.shards_dir else tbl_dir.parent / "shards"
    if args.only_done:
        tbls = [p.with_suffix(".tbl") for p in sorted(tbl_dir.glob("*.done"))]
        tbls = [p for p in tbls if p.exists() and p.stat().st_size > 0]
    else:
        tbls = sorted(tbl_dir.glob("*.tbl"))

    best: dict[str, dict] = {}
    clan_best: dict[tuple[str, str], dict] = {}
    all_hits: list[dict] = []
    seen: set[str] = set()
    n_lines = 0
    for tbl in tbls:
        with open(tbl) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                n_lines += 1
                fields = line.split()
                # Infernal --tblout (per-target) columns:
                # 0 target name (CM), 1 target acc, 2 query name (sequence),
                # 3 query acc, ..., 9 strand, 10 trunc, 11 pass, 12 gc,
                # 13 bias, 14 score, 15 E-value, 16 inc, 17 description
                if len(fields) < 17:
                    continue
                target_name = fields[0]
                query_name = fields[2]
                accession = fields[1]
                try:
                    evalue = float(fields[15])
                    score = float(fields[14])
                    bias = float(fields[13])
                except ValueError:
                    continue
                seen.add(query_name)
                clan = clan_of.get(accession, accession)
                all_hits.append(
                    {
                        "query": query_name,
                        "target_name": target_name,
                        "family_accession": accession,
                        "clan": clan,
                        "evalue": evalue,
                        "score": score,
                        "bias": bias,
                        "mdl_from": fields[5],
                        "mdl_to": fields[6],
                        "seq_from": fields[7],
                        "seq_to": fields[8],
                        "strand": fields[9],
                        "shard": tbl.stem,
                    }
                )
                key = (query_name, clan)
                cur_clan = clan_best.get(key)
                if cur_clan is None or score > cur_clan["score"]:
                    clan_best[key] = {
                        "query": query_name,
                        "clan": clan,
                        "family_accession": accession,
                        "family_name": names.get(accession, target_name),
                        "evalue": evalue,
                        "score": score,
                        "n_hits_in_clan": 0,
                    }
                cur = best.get(query_name)
                if cur is None or score > cur["score"]:
                    best[query_name] = {
                        "query": query_name,
                        "target_name": target_name,
                        "family_accession": accession,
                        "family_name": names.get(accession, target_name),
                        "clan": clan,
                        "evalue": evalue,
                        "score": score,
                        "bias": bias,
                    }

    for hit in all_hits:
        clan_best[(hit["query"], hit["clan"])]["n_hits_in_clan"] += 1

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    cols = [
        "source",
        "query",
        "target_name",
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

        all_ids: list[str] = []
        if args.only_done:
            for tbl in tbls:
                shard_fa = shards_dir / f"{tbl.stem}.fa"
                if not shard_fa.exists():
                    continue
                with open(shard_fa) as fh:
                    for line in fh:
                        if line.startswith(">"):
                            all_ids.append(line[1:].strip().split()[0])
        else:
            inp = Path(args.input_fasta)
            if inp.exists():
                with open(inp) as fh:
                    for line in fh:
                        if line.startswith(">"):
                            all_ids.append(line[1:].strip().split()[0])
        nohit = [q for q in all_ids if q not in best]
        for query in nohit:
            source = query.split(":", 1)[0] if ":" in query else "unknown"
            writer.writerow(
                {
                    "source": source,
                    "query": query,
                    "target_name": "NOHIT",
                    "family_accession": "NOHIT",
                    "family_name": "NOHIT",
                    "clan": "NOHIT",
                    "evalue": "",
                    "score": "",
                    "bias": "",
                }
            )

    with open(out.parent / "clan_hits.tsv", "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "query",
                "clan",
                "family_accession",
                "family_name",
                "evalue",
                "score",
                "n_hits_in_clan",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for key in sorted(clan_best):
            writer.writerow(clan_best[key])

    with open(out.parent / "all_hits.tsv", "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "query",
                "target_name",
                "family_accession",
                "clan",
                "evalue",
                "score",
                "bias",
                "mdl_from",
                "mdl_to",
                "seq_from",
                "seq_to",
                "strand",
                "shard",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for hit in sorted(all_hits, key=lambda h: (h["query"], -h["score"])):
            writer.writerow(hit)

    clans_per_query: dict[str, int] = {}
    for query, _clan in clan_best:
        clans_per_query[query] = clans_per_query.get(query, 0) + 1

    summary = {
        "tables_parsed": len(tbls),
        "tbl_lines": n_lines,
        "queries_with_hit": len(seen),
        "queries_written": len(best),
        "queries_nohit": len(all_ids) - len(best) if all_ids else 0,
        "queries_covered": len(all_ids),
        "clan_pairs": len(clan_best),
        "distinct_clans_hit": len({c for _, c in clan_best}),
        "multi_clan_queries": sum(1 for n in clans_per_query.values() if n > 1),
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
