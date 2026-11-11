#!/usr/bin/env python3
"""Cross-check the polynomial maximum-nested-subset solver on the real labels.

PI Round 4g §2 asked for three things before the DP replaces the old solver:

1. agreement with the exhaustive solver on every crossing component of 16 pairs
   or fewer, over the actual label sources;
2. a statement of what happened to components larger than 16;
3. a report of how many components exist at all, so the reader can see the DP is
   not being exercised only on toy input.

Point 2 is answered both ways: the old code path is reproduced here (greedy by
shortest span, then lexicographic) and compared with the DP on any component
that exceeds the old limit. The reproduction uses ``crosses`` exactly as the old
code did, so it also shows whether that path could emit a base with two
partners.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import re
import sys
import time
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rff.analysis.labels import (  # noqa: E402
    conflict_components,
    conflicts,
    crosses,
    max_nested_subset_enumeration,
    max_nested_subset_of_arcs,
)
from rff.data.io import pairs_from_dbn  # noqa: E402

_spec = importlib.util.spec_from_file_location(
    "build_rfam_index", ROOT / "scripts" / "build_rfam_index.py"
)
rfam_builder = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(rfam_builder)

CANONICAL = {("A", "U"), ("U", "A"), ("G", "C"), ("C", "G"), ("G", "U"), ("U", "G")}
EXACT_LIMIT = 16


def is_matching(pairs) -> bool:
    used = [p for pair in pairs for p in pair]
    return len(used) == len(set(used))


def old_greedy(component: list[tuple[int, int]]) -> set[tuple[int, int]]:
    """The pre-Round-4g fallback, reproduced verbatim in behaviour."""
    chosen: set[tuple[int, int]] = set()
    for p in sorted(component, key=lambda q: (q[1] - q[0], q)):
        if not any(crosses(p, kept) for kept in chosen):
            chosen.add(p)
    return chosen


def bprna_structures(path: Path):
    data = json.loads(path.read_text(encoding="utf-8"))
    for rid, rec in data.items():
        seq = rec.get("sequence") or ""
        pairs = {(int(i), int(j)) for i, j in rec.get("structure", []) if i != j}
        canonical = {
            (i, j)
            for i, j in pairs
            if 0 <= i < len(seq) and 0 <= j < len(seq) and (seq[i], seq[j]) in CANONICAL
        }
        yield rid, canonical


def rfam_structures(path: Path, max_per_family: int):
    with rfam_builder.open_maybe_gzip(path) as handle:
        for meta, seqs, ss_cons in rfam_builder.parse_stockholm(handle):
            acc = meta.get("AC", "").split()[0] if meta.get("AC") else ""
            if not acc or not ss_cons:
                continue
            items = list(seqs.items())
            if max_per_family:
                items = items[:max_per_family]
            for seq_name, aligned in items:
                ungapped = re.sub(r"[-._~]", "", aligned).upper().replace("T", "U")
                dbn = rfam_builder.project_ss(ss_cons, aligned)
                if not ungapped or len(ungapped) != len(dbn):
                    continue
                try:
                    raw = pairs_from_dbn(dbn)
                except ValueError:
                    continue
                canonical = {
                    (i, j) for i, j in raw if (ungapped[i], ungapped[j]) in CANONICAL
                }
                yield f"{acc}|{seq_name}", canonical


def scan(name: str, structures, limit: int, stats: dict) -> None:
    sizes = Counter()
    n_seq = 0
    n_with_conflict = 0
    n_components = 0
    n_checked = 0
    mismatches = []
    large = []
    n_large_greedy_worse = 0
    n_large_greedy_invalid = 0
    total_large_gap = 0
    started = time.time()
    for rid, pairs in structures:
        n_seq += 1
        if limit and n_seq > limit:
            break
        free, components = conflict_components(pairs)
        if components:
            n_with_conflict += 1
        for component in components:
            n_components += 1
            sizes[len(component)] += 1
            solved = max_nested_subset_of_arcs(component)
            if len(component) <= EXACT_LIMIT:
                n_checked += 1
                reference = max_nested_subset_enumeration(component)
                if reference != solved:
                    mismatches.append(
                        {
                            "id": rid,
                            "n_arcs": len(component),
                            "reference": sorted(reference),
                            "dp": sorted(solved),
                        }
                    )
            else:
                greedy = old_greedy(component)
                if len(greedy) < len(solved):
                    n_large_greedy_worse += 1
                    total_large_gap += len(solved) - len(greedy)
                if not is_matching(greedy):
                    n_large_greedy_invalid += 1
                large.append(
                    {
                        "id": rid,
                        "n_arcs": len(component),
                        "dp_size": len(solved),
                        "old_greedy_size": len(greedy),
                        "old_greedy_is_valid_matching": is_matching(greedy),
                        "dp_is_valid_matching": is_matching(solved),
                    }
                )
            if not is_matching(solved):
                mismatches.append({"id": rid, "error": "DP returned a base with two partners"})

    stats[name] = {
        "n_sequences_scanned": min(n_seq, limit) if limit else n_seq,
        "sequence_limit": limit,
        "n_sequences_with_a_conflict": n_with_conflict,
        "n_components": n_components,
        "component_size_histogram": {str(k): v for k, v in sorted(sizes.items())},
        "max_component_size": max(sizes) if sizes else 0,
        "n_components_compared_with_enumeration": n_checked,
        "n_mismatches": len(mismatches),
        "mismatch_examples": mismatches[:5],
        "n_components_above_16": len(large),
        "n_above_16_where_old_greedy_was_suboptimal": n_large_greedy_worse,
        "total_pairs_lost_by_old_greedy": total_large_gap,
        "n_above_16_where_old_greedy_was_not_a_valid_matching": n_large_greedy_invalid,
        "above_16": large[:10],
        "seconds": round(time.time() - started, 1),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bprna", default="data/raw/bprna/data.json")
    ap.add_argument("--rfam", default="data/raw/rfam/Rfam.seed.gz")
    ap.add_argument("--max-per-family", type=int, default=20)
    ap.add_argument("--limit", type=int, default=0, help="0 = every sequence")
    ap.add_argument("--out", default="experiments/007_nested_dp/nested_subset_check.json")
    args = ap.parse_args()

    stats: dict = {}
    if (ROOT / args.bprna).exists():
        scan("bprna", bprna_structures(ROOT / args.bprna), args.limit, stats)
    else:
        stats["bprna"] = {"status": "input missing"}
    if (ROOT / args.rfam).exists():
        scan(
            "rfam_seed",
            rfam_structures(ROOT / args.rfam, args.max_per_family),
            args.limit,
            stats,
        )
    else:
        stats["rfam_seed"] = {"status": "input missing"}

    report = {
        "purpose": "does the polynomial solver agree with exhaustive search on real components?",
        "solver": "rff.analysis.labels.max_nested_subset_of_arcs",
        "reference": "max_nested_subset_enumeration, limit 16",
        "conflicts": "crossing OR shared endpoint",
        "rfam_cap_per_family": args.max_per_family,
        "sources": stats,
    }
    out = ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"out": str(out)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
