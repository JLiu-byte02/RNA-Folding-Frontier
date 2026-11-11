#!/usr/bin/env python3
"""Check cmscan assignments against labels we already trust.

ArchiveII's record IDs start with the family name, so the assigner can be
validated without extra annotation: a 5s_* record must be assigned to RF00001 or
another 5S rRNA family, a tRNA_* record to the tRNA clan, and so on. A high
agreement rate is evidence that the assigner and the parser are correct; a low
one means one of them is broken.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rff.analysis.clans import load_family_table  # noqa: E402

EXPECTED = {
    "5s": {"RF00001"},
    "16s": {"RF00177", "RF01959", "RF01960", "RF02542"},
    "23s": {"RF02541", "RF02540", "RF02543"},
    "tRNA": {"RF00005", "RF01852"},
    "tmRNA": {"RF00023", "RF01849", "RF01850", "RF01851", "RF02544"},
    "RNaseP": {"RF00009", "RF00010", "RF00011", "RF00373", "RF02357"},
    "srp": {
        "RF00017",
        "RF00169",
        "RF01502",
        "RF01570",
        "RF01854",
        "RF01855",
        "RF01856",
        "RF01857",
        "RF04183",
    },
    "grp1": {"RF00028"},
    "grp2": {"RF00029"},
    "telomerase": {"RF00024", "RF00025", "RF01050"},
}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--assignments", default="data/processed/cmscan/assignments.tsv")
    ap.add_argument("--family-table", default="data/processed/rfam/families.tsv")
    ap.add_argument("--out", default="experiments/000_e0_leakage/cmscan_validation")
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    family_table = load_family_table(ROOT / args.family_table)
    clan_of = {acc: meta["clan"] for acc, meta in family_table.items()}
    expected_clans = {
        prefix: {clan_of.get(acc, acc) for acc in accs}
        for prefix, accs in EXPECTED.items()
    }

    per_family: dict[str, Counter] = defaultdict(Counter)
    examples: dict[str, list[str]] = defaultdict(list)
    confusions: dict[str, Counter] = defaultdict(Counter)
    with open(args.assignments) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != len(header):
                continue
            row = dict(zip(header, parts))
            if row["source"] != "archiveii":
                continue
            rid = row["query"].split(":", 1)[1]
            prefix = rid.split("_", 1)[0]
            if prefix not in EXPECTED:
                per_family[prefix]["n"] += 1
                continue
            per_family[prefix]["n"] += 1
            if row["family_accession"] == "NOHIT":
                per_family[prefix]["nohit"] += 1
                continue
            if row["family_accession"] in EXPECTED[prefix]:
                per_family[prefix]["agree"] += 1
            else:
                per_family[prefix]["disagree"] += 1
                confusions[prefix][row["family_accession"]] += 1
                if len(examples[prefix]) < 5:
                    examples[prefix].append(
                        f"{rid}->{row['family_accession']}({row['family_name']})"
                    )
            if row["clan"] in expected_clans[prefix]:
                per_family[prefix]["clan_agree"] += 1
            else:
                per_family[prefix]["clan_disagree"] += 1

    report = {
        "per_family": {k: dict(v) for k, v in sorted(per_family.items())},
        "top_confusions": {
            k: dict(v.most_common(5)) for k, v in sorted(confusions.items())
        },
        "examples_of_disagreement": dict(examples),
    }
    known = {k: v for k, v in per_family.items() if k in EXPECTED}
    total = sum(v["n"] for v in known.values())
    scored = sum(v["agree"] + v["disagree"] for v in known.values())
    agree = sum(v["agree"] for v in known.values())
    clan_agree = sum(v["clan_agree"] for v in known.values())
    report["overall"] = {
        "n": total,
        "n_scored": scored,
        "n_nohit": sum(v["nohit"] for v in known.values()),
        "agree": agree,
        "agreement_rate": round(agree / scored, 4) if scored else 0.0,
        "clan_agreement_rate": round(clan_agree / scored, 4) if scored else 0.0,
        "nohit_rate": round(
            sum(v["nohit"] for v in known.values()) / total, 4
        )
        if total
        else 0.0,
    }
    (out / "summary.json").write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2)[:4000])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
