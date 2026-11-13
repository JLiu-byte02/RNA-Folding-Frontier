#!/usr/bin/env python3
"""Coverage of a model dump against its split, plus v3 label sparsity.

Read-only. Answers the three questions the common-cohort policy depends on:

1. Which split ids never reached the model dump (context truncation), how long
   are they, and which families lose every member?
2. How many dump rows carry empty or very sparse v3 truth?
3. Is there a clean length gap between the longest dumped and shortest missing
   sequence, i.e. is truncation the only reason for absence?

Usage:
  python scripts/dump_coverage.py --split SPLIT.jsonl --dump DUMP.jsonl \
      [--out REPORT.json] [--label NAME]

No prediction is read or written; this cannot leak test outcomes.
"""

import argparse
import json
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path


def read_jsonl(path):
    rows = []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if line:
                rows.append(json.loads(line))
    return rows


def _quantiles(values):
    if not values:
        return None
    values = sorted(values)
    return {
        "min": values[0],
        "p25": values[len(values) // 4],
        "median": statistics.median(values),
        "p75": values[(3 * len(values)) // 4],
        "max": values[-1],
    }


def coverage(split_path, dump_path, label=None):
    split = read_jsonl(split_path)
    dump = read_jsonl(dump_path)
    if len({r["id"] for r in split}) != len(split):
        raise SystemExit("duplicate ids in split")
    if len({r["id"] for r in dump}) != len(dump):
        raise SystemExit("duplicate ids in dump")

    seq_len = {r["id"]: len(r["seq"]) for r in split}
    family = {r["id"]: r["family"] for r in split}
    dump_ids = {r["id"] for r in dump}
    missing = [rid for rid in seq_len if rid not in dump_ids]
    extra = sorted(dump_ids - set(seq_len))

    fam_members = Counter(family.values())
    fam_missing = Counter(family[rid] for rid in missing)
    fully_missing = {
        fam: {"n_missing": fam_missing[fam], "n_total": fam_members[fam]}
        for fam in fam_missing
        if fam_missing[fam] == fam_members[fam]
    }
    partial = {
        fam: {"n_missing": fam_missing[fam], "n_total": fam_members[fam]}
        for fam in fam_missing
        if fam_missing[fam] < fam_members[fam]
    }

    n_empty = n_lt4 = n_ge4 = 0
    pair_counts = []
    for row in dump:
        truth = row.get("truth_v3")
        if truth is None:
            truth = row.get("pairs") or []
        n = len(truth)
        pair_counts.append(n)
        if n == 0:
            n_empty += 1
        if n < 4:
            n_lt4 += 1
        else:
            n_ge4 += 1

    report = {
        "label": label or Path(split_path).stem,
        "split_path": str(split_path),
        "dump_path": str(dump_path),
        "n_split": len(split),
        "n_split_families": len(fam_members),
        "n_dump": len(dump),
        "n_missing_from_dump": len(missing),
        "missing_fraction": round(len(missing) / max(1, len(split)), 6),
        "n_extra_in_dump": len(extra),
        "missing_length": _quantiles([seq_len[r] for r in missing]),
        "present_length_max": max((seq_len[r] for r in dump_ids if r in seq_len), default=None),
        "present_length_min": min((seq_len[r] for r in dump_ids if r in seq_len), default=None),
        "length_gap_ok": (
            max((seq_len[r] for r in dump_ids if r in seq_len), default=0)
            < min((seq_len[r] for r in missing), default=10**9)
        ),
        "families_fully_missing": dict(sorted(fully_missing.items())),
        "families_partially_missing": dict(sorted(partial.items())),
        "v3_n_empty_truth": n_empty,
        "v3_n_lt4_pairs": n_lt4,
        "v3_n_ge4_pairs": n_ge4,
        "v3_pair_count_quantiles": _quantiles(pair_counts),
    }
    return report


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--split", required=True)
    ap.add_argument("--dump", required=True)
    ap.add_argument("--out", default=None)
    ap.add_argument("--label", default=None)
    args = ap.parse_args(argv)

    report = coverage(args.split, args.dump, args.label)
    text = json.dumps(report, indent=2, sort_keys=True)
    if args.out:
        out = Path(args.out)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(text + "\n", encoding="utf-8")
    print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
