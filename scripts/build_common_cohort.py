#!/usr/bin/env python3
"""Build one explicit common cohort for the external-baseline table.

The defect this fixes: ``score_external_baseline.py`` silently *skips* a
sequence whose prediction was not produced (EternaFold's 1200 nt cap, the 111
RNAstructure IUPAC-skips, MXfold2's own length cap). Every baseline therefore
ended up scored on a different effective cohort, and the family-macro column
was not comparable across rows.

This script writes the intersection as an explicit id list. Every method -- our
own heads included -- is then rescored with ``--restrict-ids`` against that
file, so ``n_missing_prediction`` must be 0 and the remaining exclusion reasons
must be identical in every row.

Cohort definition (all conditions must hold):
  1. the id appears in our model dump, so our own heads can be scored;
  2. ``truth_v3`` is not null (Amendment 8 parse-failure rule) and has at
     least ``--min-true-pairs`` pairs *before* the nested projection, exactly
     as ``rff.analysis.block_eval.score_dump`` applies it;
  3. a prediction exists in every ``--preds`` file handed to this script.

Writes:
  --out          JSONL of ``{"id", "family", "length"}`` sorted by id
  --report       JSON summary: per-reason exclusion counts and, per tool, the
                 ids it lacks together with their length distribution
  --excluded-out JSONL of ``{"id", "reason", "detail"}`` for audit
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


def _resolve(path: str) -> Path:
    p = Path(path)
    return p if p.is_absolute() else ROOT / p


def refuse_embargoed(path: Path) -> None:
    if "test_novel" in str(path).lower():
        raise SystemExit(f"refusing {path}: test_novel is embargoed until Amendment 14")


def read_jsonl(path: Path) -> list[dict]:
    with open(path, encoding="utf-8") as fh:
        return [json.loads(line) for line in fh if line.strip()]


def write_jsonl(path: str, rows: list[dict]) -> None:
    out = _resolve(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("".join(json.dumps(r, sort_keys=True) + "\n" for r in rows), encoding="utf-8")


def quantiles(values: list[int]) -> dict:
    if not values:
        return {"n": 0}
    arr = np.asarray(values, dtype=float)
    return {
        "n": int(arr.size),
        "min": int(arr.min()),
        "p50": float(np.percentile(arr, 50)),
        "p90": float(np.percentile(arr, 90)),
        "max": int(arr.max()),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--split", required=True)
    ap.add_argument("--dump", required=True, help="our model dump for the same split")
    ap.add_argument("--preds", nargs="+", required=True, help="prediction JSONL per tool")
    ap.add_argument("--tool-names", nargs="+", default=None, help="labels for --preds")
    ap.add_argument("--min-true-pairs", type=int, default=4)
    ap.add_argument("--out", required=True)
    ap.add_argument("--report", required=True)
    ap.add_argument("--excluded-out", default=None)
    args = ap.parse_args()

    split_path = _resolve(args.split)
    dump_path = _resolve(args.dump)
    refuse_embargoed(split_path)
    refuse_embargoed(dump_path)
    pred_paths = [_resolve(p) for p in args.preds]
    for p in pred_paths:
        refuse_embargoed(p)
        if not p.exists():
            raise SystemExit(f"missing prediction file: {p}")
    names = args.tool_names or [p.stem for p in pred_paths]
    if len(names) != len(pred_paths):
        raise SystemExit("--tool-names must match --preds in length")

    split_rows = {row["id"]: row for row in read_jsonl(split_path)}
    dump_rows = {row["id"]: row for row in read_jsonl(dump_path)}
    if len(dump_rows) != len(read_jsonl(dump_path)):
        raise SystemExit("duplicate ids in dump")

    pred_ids: dict[str, set[str]] = {}
    for name, path in zip(names, pred_paths):
        ids = [row["id"] for row in read_jsonl(path)]
        dupes = len(ids) - len(set(ids))
        if dupes:
            raise SystemExit(f"{name}: {dupes} duplicate ids in {path}")
        pred_ids[name] = set(ids)

    status = Counter()
    excluded: dict[str, tuple[str, str]] = {}
    eligible: dict[str, dict] = {}
    for rid, row in dump_rows.items():
        truth = row.get("truth_v3")
        status[str(row.get("truth_v3_status"))] += 1
        if truth is None:
            excluded[rid] = ("parse_failure", str(row.get("truth_v3_status")))
            continue
        if len(truth) < args.min_true_pairs:
            reason = "empty_truth" if not truth else "lt_min_true_pairs"
            excluded[rid] = (reason, f"n_true={len(truth)}")
            status[f"excluded_{reason}"] += 1
            continue
        eligible[rid] = row

    missing_by_tool: dict[str, list[str]] = {
        name: [rid for rid in eligible if rid not in pred_ids[name]] for name in names
    }
    common = set(eligible)
    for name in names:
        common -= set(missing_by_tool[name])
        for rid in missing_by_tool[name]:
            excluded.setdefault(rid, ("missing_prediction", name))

    split_only = [rid for rid in split_rows if rid not in dump_rows]
    out_rows = []
    for rid in sorted(common):
        row = eligible[rid]
        src = split_rows.get(rid)
        out_rows.append(
            {
                "id": rid,
                "family": row.get("family") if src is None else src.get("family"),
                "length": int(row.get("length", len(row.get("seq", "")))),
            }
        )

    report = {
        "split": args.split,
        "dump": args.dump,
        "preds": {name: str(path) for name, path in zip(names, pred_paths)},
        "min_true_pairs": args.min_true_pairs,
        "n_split_rows": len(split_rows),
        "n_dump_rows": len(dump_rows),
        "n_split_not_in_dump": len(split_only),
        "split_not_in_dump": {
            "count": len(split_only),
            "length": quantiles(
                [
                    int(split_rows[r].get("length", len(split_rows[r]["seq"])))
                    for r in split_only
                ]
            ),
        },
        "truth_v3_status": dict(status),
        "n_label_eligible": len(eligible),
        "per_tool": {
            name: {
                "n_ids": len(pred_ids[name]),
                "n_missing_from_label_eligible": len(missing_by_tool[name]),
                "missing_length": quantiles(
                    [int(eligible[r]["length"]) for r in missing_by_tool[name]]
                ),
                "n_extra_ids_not_in_dump": len(pred_ids[name] - set(dump_rows)),
            }
            for name in names
        },
        "n_common": len(out_rows),
        "n_families": len({row["family"] for row in out_rows}),
        "common_length": quantiles([row["length"] for row in out_rows]),
        "n_excluded_total": len(excluded),
        "excluded_by_reason": dict(Counter(reason for reason, _ in excluded.values())),
    }

    write_jsonl(args.out, out_rows)
    _resolve(args.report).parent.mkdir(parents=True, exist_ok=True)
    _resolve(args.report).write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    if args.excluded_out:
        write_jsonl(
            args.excluded_out,
            [
                {"id": rid, "reason": reason, "detail": detail}
                for rid, (reason, detail) in sorted(excluded.items())
            ],
        )
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())