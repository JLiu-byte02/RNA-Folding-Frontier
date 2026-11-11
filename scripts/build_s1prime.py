#!/usr/bin/env python3
"""Build the clan-held-out splits S1'-novel and S1'-withheld.

Rules (PI Round 2, verbatim intent):

* Pool = bpRNA-1m(90) union Rfam 15.0 SEED sequences.
* Family assignment comes from one uniform assigner: cmscan against Rfam 15.0
  with ``--rfam --cut_ga --nohmmonly``. A family maps to a clan via the release's
  clan membership; a family without a clan is its own singleton clan.
* A clan is **novel** if no bpRNA sequence is assigned to any family in it. Its
  sequences therefore come only from Rfam seeds, and the whole clan is removed
  from training.
* A clan is **withheld** if bpRNA does support it: we remove it from training on
  purpose and use its bpRNA + seed sequences as an OOD test with the *same* label
  provenance as training.
* Labels: bpRNA pairs where they exist, otherwise the seed alignment's SS_cons
  projected onto the sequence. Policy: canonical (AU/UA/GC/CG) plus GU/UG.

The split is written with a manifest of sha256 hashes and the clan lists so it
can be frozen before any model sees it.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import random
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rff.analysis.clans import load_family_table  # noqa: E402
from rff.data.io import pairs_from_dbn  # noqa: E402

CANONICAL = {("A", "U"), ("U", "A"), ("G", "C"), ("C", "G"), ("G", "U"), ("U", "G")}


def load_assignments(path: Path) -> dict[str, dict]:
    rows: dict[str, dict] = {}
    if not path.exists():
        return rows
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != len(header):
                continue
            row = dict(zip(header, parts))
            rows[row["query"]] = row
    return rows


def filter_pairs(seq: str, pairs) -> list[list[int]]:
    keep = []
    for i, j in pairs:
        if (seq[i], seq[j]) in CANONICAL:
            keep.append([min(i, j), max(i, j)])
    return sorted({tuple(p) for p in keep})


def read_fasta(path: Path) -> list[tuple[str, str]]:
    out: list[tuple[str, str]] = []
    name, buf = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    out.append((name, "".join(buf)))
                name, buf = line[1:].strip().split()[0], []
            else:
                buf.append(line.strip())
    if name:
        out.append((name, "".join(buf)))
    return out


def sha256_lines(items: list[str]) -> str:
    h = hashlib.sha256()
    for item in sorted(items):
        h.update(item.encode())
    return h.hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--assignments", default="data/processed/cmscan/assignments.tsv")
    ap.add_argument("--bprna-pool", default="data/processed/bprna/pool.jsonl")
    ap.add_argument("--rfam-seeds", default="data/processed/rfam/seed_sequences.fa")
    ap.add_argument(
        "--rfam-structures", default="data/processed/rfam/seed_structures.json"
    )
    ap.add_argument("--out", default="splits/s1prime")
    ap.add_argument("--family-table", default="data/processed/rfam/families.tsv")
    ap.add_argument("--novel-clan-fraction", type=float, default=0.20)
    ap.add_argument("--seed", type=int, default=17)
    args = ap.parse_args()

    rng = random.Random(args.seed)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    assignments = load_assignments(ROOT / args.assignments)
    if not assignments:
        print("no assignments yet; nothing to do")
        return 1

    # Clan membership must come from the Rfam release, not from whatever the
    # bpRNA scan happened to hit: a clan whose members are never hit by bpRNA is
    # precisely a *novel* clan, and deriving its map from hits would drop it.
    family_table = load_family_table(ROOT / args.family_table)
    if not family_table:
        print(f"missing family table {args.family_table}")
        return 1
    fam_clans = {acc: meta["clan"] for acc, meta in family_table.items()}

    # bpRNA sequences: keep those with a family assignment
    bprna_by_clan: dict[str, list[dict]] = defaultdict(list)
    bprna_unassigned = 0
    bprna_nohit = 0
    bprna_scanned = 0
    with open(ROOT / args.bprna_pool) as fh:
        for line in fh:
            rec = json.loads(line)
            row = assignments.get(f"bprna:{rec['id']}")
            if row is None:
                bprna_unassigned += 1
                continue
            bprna_scanned += 1
            if row["family_accession"] == "NOHIT" or row["clan"] == "NOHIT":
                # NOHIT is "no family above the Rfam gathering threshold", not a
                # family. Folding it into a clan would create a pseudo-cluster.
                bprna_nohit += 1
                continue
            pairs = filter_pairs(rec["seq"], rec["pairs"])
            clan = fam_clans.get(row["family_accession"]) or row["clan"]
            bprna_by_clan[clan].append(
                {
                    "id": f"bprna:{rec['id']}",
                    "seq": rec["seq"],
                    "pairs": pairs,
                    "family": row["family_accession"],
                    "clan": clan,
                    "source": "bprna",
                }
            )

    # Rfam seeds: family known by construction, labels from SS_cons
    structures = json.loads((ROOT / args.rfam_structures).read_text())
    seeds_by_clan: dict[str, list[dict]] = defaultdict(list)
    for header, seq in read_fasta(ROOT / args.rfam_seeds):
        acc, name = header.split("|", 1) if "|" in header else (header, "")
        meta = family_table.get(acc)
        if meta is None:
            continue
        clan = meta["clan"]
        dbn = structures.get(acc, {}).get(name, "")
        pairs = []
        if dbn and len(dbn) == len(seq):
            try:
                pairs = filter_pairs(seq, pairs_from_dbn(dbn))
            except ValueError:
                pairs = []
        seeds_by_clan[clan].append(
            {
                "id": f"rfam:{acc}|{name}",
                "seq": seq,
                "pairs": pairs,
                "family": acc,
                "clan": clan,
                "source": "rfam",
            }
        )

    all_clans = set(bprna_by_clan) | set(seeds_by_clan)
    novel = sorted(c for c in all_clans if not bprna_by_clan.get(c))
    supported = sorted(c for c in all_clans if bprna_by_clan.get(c))

    n_withhold = max(1, int(round(len(supported) * args.novel_clan_fraction)))
    withheld = sorted(rng.sample(supported, n_withhold))
    retained = [c for c in supported if c not in set(withheld)]

    train, novel_test, withheld_test = [], [], []
    for clan in retained:
        train.extend(bprna_by_clan.get(clan, []))
    for clan in novel:
        novel_test.extend(seeds_by_clan.get(clan, []))
    for clan in withheld:
        withheld_test.extend(bprna_by_clan.get(clan, []))
        withheld_test.extend(seeds_by_clan.get(clan, []))

    # ID test: hold out 10% of retained bpRNA families, sequence-deduped later
    fams = sorted({rec["family"] for rec in train})
    id_fams = set(rng.sample(fams, max(1, len(fams) // 10)))
    id_test = [rec for rec in train if rec["family"] in id_fams]
    train = [rec for rec in train if rec["family"] not in id_fams]

    def write(name: str, records: list[dict]) -> str:
        path = out / f"{name}.jsonl"
        with open(path, "w") as fh:
            for rec in records:
                fh.write(json.dumps(rec) + "\n")
        return sha256_lines([rec["id"] for rec in records])

    hashes = {
        "train": write("train", train),
        "id_test": write("id_test", id_test),
        "test_novel": write("test_novel", novel_test),
        "test_withheld": write("test_withheld", withheld_test),
    }

    manifest = {
        "seed": args.seed,
        "n_clans_total": len(all_clans),
        "n_clans_novel": len(novel),
        "n_clans_withheld": len(withheld),
        "n_clans_retained": len(retained),
        "novel_clans": novel,
        "withheld_clans": withheld,
        "counts": {
            "train": len(train),
            "id_test": len(id_test),
            "test_novel": len(novel_test),
            "test_withheld": len(withheld_test),
            "bprna_scanned": bprna_scanned,
            "bprna_pool_file_not_scanned": bprna_unassigned,
            "bprna_nohit": bprna_nohit,
        },
        "pair_families": {
            "train": dict(Counter(rec["source"] for rec in train)),
            "test_novel": dict(Counter(rec["source"] for rec in novel_test)),
            "test_withheld": dict(Counter(rec["source"] for rec in withheld_test)),
        },
        "hashes": hashes,
        "rules": {
            "pairs": (
                "canonical AU/UA/GC/CG plus GU/UG; non-canonical dropped; "
                "crossing (pseudoknotted) pairs RETAINED "
                "(corrected 2026-11-12, docs/corrections.md)"
            ),
            "novelty": "clan has no bpRNA sequence assigned by cmscan --rfam --cut_ga",
            "nohit_policy": "sequences without a cut_ga hit are excluded, never treated as a clan",
            "clan_source": args.family_table,
            "withheld_fraction": args.novel_clan_fraction,
        },
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2))
    print(json.dumps({k: v for k, v in manifest.items() if k != "novel_clans"}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
