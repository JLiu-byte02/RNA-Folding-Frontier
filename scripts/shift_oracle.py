"""Oracle upper bound for shift correction (diagnostic, not a method).

For every false-positive pair within Manhattan distance <=2 of an unpredicted
true pair, assume an oracle could snap it onto that truth pair. Recompute F1 for
that sequence and aggregate. This is an upper bound: it ignores crossing
constraints, competing FPs mapping to the same truth pair (deduplicated), and
any cost of choosing shifts.

If the oracle gain is small, the near-miss lead is a dead end; if large, a
shift-corrected decoder is worth a cycle. No claim is made until the PI sets a
bar.

    python scripts/shift_oracle.py --split ... --cohort ... \
        --head-preds ... --out ...
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from rff.analysis.labels import max_nested_subset  # noqa: E402
import score_external_baseline as seb  # noqa: E402


def _resolve(path: str) -> Path:
    p = Path(path)
    return p if p.is_absolute() else ROOT / p


def parse_preds(path: Path) -> dict[str, set[tuple[int, int]]]:
    out = {}
    for row in seb.read_jsonl(path):
        if "pairs" in row:
            out[row["id"]] = {(int(a), int(b)) for a, b in row["pairs"]}
        elif "dbn" in row:
            from rff.data.io import pairs_from_dbn

            out[row["id"]] = set(pairs_from_dbn(row["dbn"]))
        else:
            raise SystemExit(f"{path}: row has neither dbn nor pairs")
    return out


def f1(tp: int, fp: int, fn: int) -> float:
    return seb.f1_from_counts(tp, fp, fn)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--split", required=True)
    ap.add_argument("--cohort", required=True)
    ap.add_argument("--head-preds", required=True, help="deployed head predictions (JSONL with dbn or pairs)")
    ap.add_argument("--labels", default=seb.LABELS_V3)
    ap.add_argument("--min-true-pairs", type=int, default=4)
    ap.add_argument("--radius", type=int, default=2)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    split_rows = {r["id"]: r for r in seb.read_jsonl(_resolve(args.split))}
    cohort = seb.read_jsonl(_resolve(args.cohort))
    labels = seb.load_labels_v3(_resolve(args.labels))
    preds = parse_preds(_resolve(args.head_preds))

    fam_rows: dict[str, list[tuple[float, float]]] = defaultdict(list)
    tot = {"tp": 0, "fp": 0, "fn": 0, "otp": 0, "ofp": 0}
    n_shiftable = 0
    n_fp_total = 0
    for row in cohort:
        rid = row["id"]
        src = split_rows.get(rid)
        key = rid[5:] if rid.startswith("rfam:") else rid
        if key in labels:
            raw = labels[key]
            if raw is None:
                continue
        elif src is not None:
            raw = {(int(a), int(b)) for a, b in src["pairs"]}
        else:
            continue
        if len(raw) < args.min_true_pairs or rid not in preds:
            continue
        truth = max_nested_subset(raw)
        pred = preds[rid]
        tp = len(pred & truth)
        fp = len(pred - truth)
        fn = len(truth - pred)
        # recoverable: unpredicted truth pairs with at least one FP within radius
        recoverable = set()
        for (i, j) in pred - truth:
            best = None
            best_d = args.radius + 1
            for (a, b) in truth - pred:
                d = abs(i - a) + abs(j - b)
                if d < best_d:
                    best_d, best = d, (a, b)
            if best is not None and best_d <= args.radius:
                recoverable.add(best)
        shift = len(recoverable)
        n_shiftable += shift
        n_fp_total += fp
        otp, ofp = tp + shift, fp - shift
        tot["tp"] += tp; tot["fp"] += fp; tot["fn"] += fn
        tot["otp"] += otp; tot["ofp"] += ofp
        fam = row.get("family") or (src["family"] if src else "NA")
        fam_rows[fam].append((f1(tp, fp, fn), f1(otp, ofp, fn)))

    base_macro = float(np.mean([np.mean([a for a, _ in v]) for v in fam_rows.values()]))
    orac_macro = float(np.mean([np.mean([b for _, b in v]) for v in fam_rows.values()]))
    out = {
        "split": args.split,
        "n_families": len(fam_rows),
        "radius": args.radius,
        "head_pred_file": args.head_preds,
        "n_fp": n_fp_total,
        "n_shiftable": n_shiftable,
        "shiftable_frac_of_fp": round(n_shiftable / n_fp_total, 4) if n_fp_total else 0,
        "baseline": {
            "macro_f1": round(base_macro, 4),
            "micro_f1": round(f1(tot["tp"], tot["fp"], tot["fn"]), 4),
        },
        "oracle_shift": {
            "macro_f1": round(orac_macro, 4),
            "micro_f1": round(f1(tot["otp"], tot["ofp"], tot["fn"]), 4),
        },
        "delta_macro": round(orac_macro - base_macro, 4),
    }
    dest = _resolve(args.out)
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps(out, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
