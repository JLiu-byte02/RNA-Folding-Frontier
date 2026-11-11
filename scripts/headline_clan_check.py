"""Does BPfold's five 'unseen family' sets leak at clan level?

BPfold highlights Skip-rope, Twister, c-di-GMP, RAGATH and Cobalamin as unseen
families.  Each maps to one or more Rfam families.  Two failure modes matter:

``same_family``
    a training sequence is assigned to the held-out family itself;
``clan_sibling``
    a training sequence is assigned to a *different* family of the same clan,
    which shares the consensus structure and is the level at which these
    questions are biologically equivalent.

Single-label assignment hides the second mode, so this reads ``clan_hits.tsv``
(best hit per query x clan) instead of ``assignments.tsv``.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rff.analysis.clans import load_family_table  # noqa: E402

HEADLINE = {
    "skipping-rope": ["RF02924"],
    "Twister": ["RF02684", "RF02681"],
    "c-di-GMP": ["RF01051", "RF01786"],
    "RAGATH-18": ["RF03064"],
    "Cobalamin": ["RF00174"],
}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--clan-hits", default="data/processed/cmscan/clan_hits.tsv")
    ap.add_argument("--family-table", default="data/processed/rfam/families.tsv")
    ap.add_argument("--out", default="experiments/000_e0_leakage/headline_clans")
    ap.add_argument("--sources", default="bprna,ts0")
    args = ap.parse_args()

    fams = load_family_table(ROOT / args.family_table)
    clan_of = {acc: meta["clan"] for acc, meta in fams.items()}
    members: dict[str, list[str]] = defaultdict(list)
    for acc, clan in clan_of.items():
        members[clan].append(acc)

    targets = {}
    for name, accs in HEADLINE.items():
        targets[name] = {
            "accs": set(accs),
            "clans": {clan_of.get(a, a) for a in accs},
        }

    sources = {s.strip() for s in args.sources.split(",") if s.strip()}
    stats: dict[str, dict[str, Counter]] = {
        name: {s: Counter() for s in sources} for name in HEADLINE
    }
    sibling_families: dict[str, Counter] = {name: Counter() for name in HEADLINE}
    examples: dict[str, dict[str, list[str]]] = {
        name: defaultdict(list) for name in HEADLINE
    }
    source_totals: Counter = Counter()

    with open(args.clan_hits) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            source = row["query"].split(":", 1)[0]
            if source in sources:
                source_totals[source] += 1
            for name, tgt in targets.items():
                if row["clan"] not in tgt["clans"]:
                    continue
                bucket = stats[name][source]
                if row["family_accession"] in tgt["accs"]:
                    bucket["same_family"] += 1
                else:
                    bucket["clan_sibling"] += 1
                    sibling_families[name][
                        f"{row['family_accession']}({row['family_name']})"
                    ] += 1
                if len(examples[name][source]) < 5:
                    examples[name][source].append(
                        f"{row['query']}->{row['family_accession']}"
                        f"({row['family_name']},{row['clan']},score={row['score']})"
                    )

    report = {
        "sources": sorted(sources),
        "queries_per_source": dict(source_totals),
        "mapping": {
            name: {
                "rfam": sorted(accs),
                "clans": sorted(targets[name]["clans"]),
                "clan_members": {
                    c: sorted(members.get(c, [])) for c in targets[name]["clans"]
                },
            }
            for name, accs in HEADLINE.items()
        },
        "training_hits": {
            name: {s: dict(c) for s, c in per_source.items() if c}
            for name, per_source in stats.items()
        },
        "sibling_families_hit": {
            name: dict(c.most_common(10))
            for name, c in sibling_families.items()
            if c
        },
        "examples": {
            name: {s: v for s, v in per_source.items() if v}
            for name, per_source in examples.items()
            if any(per_source.values())
        },
        "verdict": {},
    }
    for name, per_source in stats.items():
        same = sum(c["same_family"] for c in per_source.values())
        sib = sum(c["clan_sibling"] for c in per_source.values())
        report["verdict"][name] = {
            "same_family_hits": same,
            "clan_sibling_hits": sib,
            "leaks": bool(same or sib),
        }

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "summary.json").write_text(json.dumps(report, indent=2))
    print(json.dumps(report["verdict"], indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
