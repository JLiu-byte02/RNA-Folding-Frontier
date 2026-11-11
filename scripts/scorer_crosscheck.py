#!/usr/bin/env python3
"""Cross-check the project scorer against an independent implementation.

PI Round 4e: "score the same RNAfold predictions with our scorer and with one
independent implementation, and require exact agreement per sequence."

This runs on the saved RNAfold predictions for the local ArchiveII mirror
(`experiments/004_turner_control/archiveii_rnafold_predictions.jsonl`). Any
disagreement is printed with the sequence id and both counts; the run exits
non-zero if a single sequence disagrees, so it can gate later work.

The ±1 slip variant is computed and reported separately. It is a sensitivity
diagnostic (how much F1 a one-base offset destroys), never a primary number.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rff.analysis.scoring import prf  # noqa: E402
from rff.decoders.vienna_soft import parse_dotbracket as project_parse  # noqa: E402

_spec = importlib.util.spec_from_file_location(
    "independent_scorer", ROOT / "scripts" / "independent_scorer.py"
)
independent = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(independent)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--archive", default="data/raw/archiveII/data.json")
    ap.add_argument(
        "--predictions",
        default="experiments/004_turner_control/archiveii_rnafold_predictions.jsonl",
    )
    ap.add_argument("--out", default="experiments/004_turner_control/scorer_crosscheck.json")
    args = ap.parse_args()

    data = json.loads((ROOT / args.archive).read_text(encoding="utf-8"))
    preds = [
        json.loads(line)
        for line in (ROOT / args.predictions).read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]

    per_family: dict[str, list[int]] = defaultdict(lambda: [0, 0, 0])
    slip: dict[int, list[int]] = {1: [0, 0, 0], -1: [0, 0, 0]}
    mismatches = []
    n = 0
    for row in preds:
        rec = data.get(row["id"])
        if rec is None:
            mismatches.append({"id": row["id"], "error": "id not in archive"})
            continue
        truth = {tuple(p) for p in rec["structure"]}
        db = row["db"]
        a = project_parse(db)
        b = independent.pairs_from_db(db)
        if a != b:
            mismatches.append(
                {"id": row["id"], "error": "pair set", "project": sorted(a), "independent": sorted(b)}
            )
            continue
        ours = prf(*independent.counts(a, truth))
        theirs = independent.counts(a, truth)
        if (ours["tp"], ours["fp"], ours["fn"]) != theirs:
            mismatches.append(
                {
                    "id": row["id"],
                    "error": "counts",
                    "project": [ours["tp"], ours["fp"], ours["fn"]],
                    "independent": list(theirs),
                }
            )
            continue
        f1_theirs = independent.f1_from_counts(*theirs)
        if abs(float(f1_theirs) - ours["f1"]) > 1e-12:
            mismatches.append(
                {"id": row["id"], "error": "f1", "project": ours["f1"], "independent": str(f1_theirs)}
            )
            continue
        fam = row["family"]
        tp, fp, fn = theirs
        bucket = per_family[fam]
        bucket[0] += tp
        bucket[1] += fp
        bucket[2] += fn
        for delta, acc in slip.items():
            s_tp, s_fp, s_fn = independent.counts(independent.shift_pairs(a, delta), truth)
            acc[0] += s_tp
            acc[1] += s_fp
            acc[2] += s_fn
        n += 1

    macro_project = sum(
        (prf(*triple)["f1"] for triple in per_family.values()), start=0.0
    ) / len(per_family)
    macro_independent = independent.family_macro_f1(per_family)
    micro = prf(
        sum(v[0] for v in per_family.values()),
        sum(v[1] for v in per_family.values()),
        sum(v[2] for v in per_family.values()),
    )

    def slip_entry(delta: int) -> dict:
        tp, fp, fn = slip[delta]
        s = prf(tp, fp, fn)
        return {
            "shift": delta,
            "precision": round(s["precision"], 4),
            "recall": round(s["recall"], 4),
            "f1": round(s["f1"], 4),
        }

    report = {
        "status": "scorer cross-check, pipeline exercise (not a protocol number)",
        "n_sequences_compared": n,
        "n_mismatches": len(mismatches),
        "mismatches": mismatches[:20],
        "agreement": {
            "pair_sets_exact": len(mismatches) == 0,
            "family_macro_f1_project": round(float(macro_project), 12),
            "family_macro_f1_independent": str(macro_independent),
            "family_macro_f1_equal": abs(float(macro_independent) - macro_project) < 1e-12,
            "micro_f1": round(micro["f1"], 4),
        },
        "sensitivity_pm1": [slip_entry(1), slip_entry(-1)],
        "pseudoknot_policy": (
            "no removal: true pairs are scored as given, including crossing ones; "
            "predictions are always non-crossing, so a pseudoknotted true pair that "
            "is predicted in the wrong register counts as FP+FN. No separate "
            "pseudoknot document exists; this records the implemented behaviour."
        ),
        "inputs": {"archive": args.archive, "predictions": args.predictions},
    }
    out = ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: report[k] for k in ("n_sequences_compared", "n_mismatches", "agreement")}, indent=1))
    return 0 if not mismatches else 1


if __name__ == "__main__":
    raise SystemExit(main())
