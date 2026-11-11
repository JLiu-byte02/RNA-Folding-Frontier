#!/usr/bin/env python3
"""Turn the Rfam SEED alignment dump into the tables the split builder needs.

Outputs (under ``--out``):
  families.tsv        accession, name, clan, description, n_seed
  clan_map.tsv        clan -> family accession (long format)
  seed_sequences.fa   one record per (family, seed sequence), ungapped
  seed_structures.json  accession -> {seq_name: dot-bracket projected from SS_cons}

The SEED dump does **not** carry ``#=GF CLAN``, so clan membership is read from
the release's ``clan_membership.txt`` (Rfam 15.0: 147 clans, 463 of 4,178
families). ``SS_cons`` gives the consensus secondary structure in alignment
columns; projecting it onto each ungapped sequence yields per-sequence
dot-bracket labels of the same type that bpRNA derives from comparative
analysis, which is what S1' needs for the families that bpRNA-1m does not cover.
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
from pathlib import Path


def open_maybe_gzip(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, encoding="utf-8", errors="replace")


def project_ss(ss_cons: str, aligned: str) -> str:
    """Project a consensus structure onto one ungapped sequence."""
    out = []
    for col, ch in enumerate(ss_cons):
        if col >= len(aligned):
            break
        if aligned[col] in "-._~":
            continue
        out.append(ch if ch in "()[]{}<>" else ".")
    return "".join(out)


def parse_stockholm(handle):
    meta: dict[str, str] = {}
    seqs: dict[str, str] = {}
    ss_cons = ""
    for line in handle:
        line = line.rstrip("\n")
        if line.startswith("//"):
            yield meta, seqs, ss_cons
            meta, seqs, ss_cons = {}, {}, ""
            continue
        if not line or line.startswith("#=GF"):
            if line.startswith("#=GF"):
                parts = line[5:].strip().split(None, 1)
                if len(parts) == 2:
                    meta.setdefault(parts[0], parts[1].strip())
            continue
        if line.startswith("#=GC"):
            tag, _, value = line[5:].partition(" ")
            if tag.strip() == "SS_cons":
                ss_cons += value.strip()
            continue
        if line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        name, chunk = parts[0], parts[1]
        if set(chunk) <= set(".-_~"):
            continue
        seqs[name] = seqs.get(name, "") + chunk


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", default="data/raw/rfam/Rfam.seed.gz")
    ap.add_argument("--clan-membership", default="data/raw/rfam/clan_membership.txt")
    ap.add_argument("--clan-table", default="data/raw/rfam/clan.txt")
    ap.add_argument("--out", default="data/processed/rfam")
    ap.add_argument("--max-seqs-per-family", type=int, default=0)
    args = ap.parse_args()

    src = Path(args.seed)
    clan_of: dict[str, str] = {}
    clan_membership = Path(args.clan_membership)
    if clan_membership.exists():
        for line in clan_membership.read_text().splitlines():
            parts = line.split("\t")
            if len(parts) >= 2 and parts[0].startswith("CL"):
                clan_of[parts[1]] = parts[0]
    clan_names: dict[str, str] = {}
    clan_table = Path(args.clan_table)
    if clan_table.exists():
        for line in clan_table.read_text().splitlines():
            parts = line.split("\t")
            if len(parts) >= 2:
                clan_names[parts[0]] = parts[1]

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    fam_rows = []
    clan_rows = []
    structures: dict[str, dict[str, str]] = {}

    with open_maybe_gzip(src) as handle, open(out / "seed_sequences.fa", "w") as fa:
        for meta, seqs, ss_cons in parse_stockholm(handle):
            acc = meta.get("AC", "").split()[0] if meta.get("AC") else ""
            if not acc:
                continue
            name = meta.get("ID", acc).split()[0]
            clan = clan_of.get(acc, "SINGLETON")
            desc = meta.get("DE", "").strip()
            fam_rows.append((acc, name, clan, desc, len(seqs)))
            if clan != "SINGLETON":
                clan_rows.append((clan, acc, name))

            fam_struct: dict[str, str] = {}
            items = list(seqs.items())
            if args.max_seqs_per_family:
                items = items[: args.max_seqs_per_family]
            for seq_name, aligned in items:
                ungapped = re.sub(r"[-._~]", "", aligned).upper().replace("T", "U")
                if not ungapped:
                    continue
                fa.write(f">{acc}|{seq_name}\n{ungapped}\n")
                if ss_cons:
                    fam_struct[seq_name] = project_ss(ss_cons, aligned)
            structures[acc] = fam_struct

    with open(out / "families.tsv", "w") as fh:
        fh.write("accession\tname\tclan\tdescription\tn_seed\n")
        for row in fam_rows:
            fh.write("\t".join(str(x) for x in row) + "\n")

    with open(out / "clan_map.tsv", "w") as fh:
        fh.write("clan\taccession\tname\n")
        for row in clan_rows:
            fh.write("\t".join(row) + "\n")

    (out / "seed_structures.json").write_text(json.dumps(structures))

    n_clans = len({r[0] for r in clan_rows})
    print(f"families: {len(fam_rows)}")
    print(f"families in clans: {len(clan_rows)} across {n_clans} clans")
    print(f"singletons: {sum(1 for r in fam_rows if r[2] == 'SINGLETON')}")
    print(f"seed sequences written: {sum(1 for _ in open(out / 'seed_sequences.fa')) // 2}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
