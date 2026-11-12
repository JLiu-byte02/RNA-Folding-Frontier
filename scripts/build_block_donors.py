#!/usr/bin/env python3
"""Full-length donor strings for the P2-block intervention (Amendment 2 A2.2).

``build_conditions.py`` writes sequences in which only the region *outside one
central window* is replaced. P2-block needs the complementary object: a donor
value for **every** position, so that any block can be embedded with its own
native flank and donor content elsewhere. This script writes exactly that.

Donors, per record, all length-preserving:

``b``  dinucleotide-preserving shuffle of the whole sequence.
``c``  the aligned same-family homolog below 60% identity, projected onto the
       target; positions whose column aligns to a gap keep the native base.
       Only records with a verified homolog are written, matching
       ``build_conditions.py``.
``d``  a contiguous same-length slice from a sequence of a different clan
       (tiled from the start if the donor is shorter than the target).
``e``  ``N`` at every position.

The donor file is ``{"<id>": "<donor string>"}``. Keeping it separate from the
embedding cache means the (expensive) embedding pass can be re-run for a
different block/flank geometry without regenerating or re-drawing donors, which
keeps the comparison between geometries on the same donor content.
"""

from __future__ import annotations

import argparse
import json
import random
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from build_conditions import dinuc_shuffle, positional_map  # noqa: E402


def load_jsonl(path: Path) -> list[dict]:
    return [json.loads(line) for line in open(path) if line.strip()]


def build_donors(
    condition: str,
    recs: list[dict],
    *,
    rng: random.Random,
    donor_pool: list[dict] | None = None,
    homolog_map: dict[str, dict] | None = None,
    seed_validated: bool = True,
) -> dict[str, str]:
    donors: dict[str, str] = {}
    by_clan: dict[str, list[str]] = defaultdict(list)
    for rec in donor_pool if donor_pool is not None else recs:
        by_clan[rec.get("clan", "?")].append(rec["seq"])
    clan_keys = sorted(by_clan)

    if condition == "c":
        if not homolog_map:
            raise SystemExit("condition c needs --homolog-map")
        seq_by_id = {r["id"]: r["seq"] for r in (donor_pool or recs)}
        for rec in recs:
            hit = homolog_map.get(rec["id"])
            if not hit:
                continue
            donor = seq_by_id.get(hit["homolog"])
            if donor is None:
                continue
            target = rec["seq"]
            mapping = positional_map(target, donor)
            out = [
                donor[j] if 0 <= j < len(donor) else target[pos]
                for pos, j in enumerate(mapping)
            ]
            donors[rec["id"]] = "".join(out)
        if not donors:
            raise SystemExit(
                "condition c produced no donors: the homolog map's targets are "
                "not present in --donor-pool. For the per-split homolog maps, "
                "pass the same split as the donor pool."
            )
        return donors

    for rec in recs:
        seq = rec["seq"]
        if condition == "b":
            donors[rec["id"]] = dinuc_shuffle(seq, rng)
        elif condition == "d":
            own = rec.get("clan", "?")
            others = [c for c in clan_keys if c != own]
            if not others:
                raise SystemExit("condition d needs at least two clans in the pool")
            donor = rng.choice(by_clan[rng.choice(others)])
            if len(donor) >= len(seq):
                start = rng.randrange(0, len(donor) - len(seq) + 1)
                donors[rec["id"]] = donor[start : start + len(seq)]
            else:
                donors[rec["id"]] = (donor * (len(seq) // len(donor) + 1))[: len(seq)]
        elif condition == "e":
            donors[rec["id"]] = "N" * len(seq)
        else:
            raise ValueError(f"unknown condition {condition!r}")
    return donors


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--split", required=True, help="jsonl of the target records")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--conditions", default="b,d,e")
    ap.add_argument("--seed", type=int, default=17)
    ap.add_argument("--donor-pool", default=None, help="jsonl to draw (d) donors from")
    ap.add_argument("--homolog-map", default=None, help="for condition c")
    args = ap.parse_args()

    recs = load_jsonl(Path(args.split))
    pool = load_jsonl(Path(args.donor_pool)) if args.donor_pool else None
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    rng = random.Random(args.seed)

    for cond in [c.strip() for c in args.conditions.split(",") if c.strip()]:
        donors = build_donors(
            cond,
            recs,
            rng=rng,
            donor_pool=pool,
            homolog_map=json.loads(Path(args.homolog_map).read_text())
            if args.homolog_map and cond == "c"
            else None,
        )
        for rec in recs:
            if rec["id"] in donors:
                assert len(donors[rec["id"]]) == len(rec["seq"]), rec["id"]
        path = out_dir / f"{cond}.json"
        path.write_text(json.dumps(donors))
        changed = sum(
            1
            for rec in recs
            if rec["id"] in donors
            and sum(a != b for a, b in zip(rec["seq"], donors[rec["id"]])) > 0
        )
        print(
            f"{cond}: {len(donors)}/{len(recs)} records, "
            f"{changed} with >=1 changed base -> {path}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
