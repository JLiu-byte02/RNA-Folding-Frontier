#!/usr/bin/env python3
"""Turner-only control on the local ArchiveII mirror.

This is deliberately **not** a protocol number. The frozen split files live only
on the cluster (Amendment 4 §A4.5), so nothing run here can be compared with the
model table. What this run does establish, on 3,370 real sequences:

* the ArchiveII reader, the dot-bracket parser and the pair scorer agree with
  each other end to end;
* RNAfold's family-macro F1 on this mirror is available as a reference for the
  size of the "classical method" baseline;
* the pair-count-vs-length relation Flamm et al. describe can be checked on a
  corpus we control.

Output is written with ``"status": "pipeline exercise, not a protocol number"``
so it cannot be mistaken for an evaluation result later.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rff.analysis.scoring import prf  # noqa: E402
from rff.decoders.vienna_soft import parse_dotbracket  # noqa: E402


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def length_band(n: int) -> str:
    if n < 80:
        return "lt80"
    if n < 150:
        return "80-149"
    if n < 300:
        return "150-299"
    return "ge300"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--archive", default="data/raw/archiveII/data.json")
    ap.add_argument("--out", default="experiments/004_turner_control/archiveii_rnafold.json")
    ap.add_argument("--limit", type=int, default=0, help="0 = all records")
    ap.add_argument("--save-db", default="", help="optional jsonl of id -> predicted dot-bracket")
    args = ap.parse_args()

    import RNA

    path = ROOT / args.archive
    data = json.loads(path.read_text(encoding="utf-8"))
    items = list(data.items())
    if args.limit:
        items = items[: args.limit]

    per_family: dict[str, list[int]] = defaultdict(lambda: [0, 0, 0])
    per_band: dict[str, list[int]] = defaultdict(lambda: [0, 0, 0])
    rows = []
    db_rows = []
    n_skipped = 0
    for rid, rec in items:
        seq = rec.get("sequence")
        truth = {tuple(p) for p in rec.get("structure", [])}
        if not seq:
            n_skipped += 1
            continue
        family = rid.split("_")[0]
        structure, energy = RNA.fold(seq)
        pred = parse_dotbracket(structure)
        if args.save_db:
            db_rows.append({"id": rid, "family": family, "length": len(seq), "db": structure})
        tp = len(pred & truth)
        fp = len(pred - truth)
        fn = len(truth - pred)
        for bucket, key in ((per_family, family), (per_band, length_band(len(seq)))):
            bucket[key][0] += tp
            bucket[key][1] += fp
            bucket[key][2] += fn
        rows.append(
            {
                "id": rid,
                "family": family,
                "length": len(seq),
                "energy": round(float(energy), 3),
                "n_true": len(truth),
                "n_pred": len(pred),
                "tp": tp,
                "fp": fp,
                "fn": fn,
            }
        )

    def table(acc: dict[str, list[int]]) -> dict:
        out = {}
        for key, (tp, fp, fn) in sorted(acc.items()):
            s = prf(tp, fp, fn)
            out[key] = {
                "n_pairs_true": tp + fn,
                "n_pairs_pred": tp + fp,
                "tp": tp,
                "fp": fp,
                "fn": fn,
                "precision": round(s["precision"], 4),
                "recall": round(s["recall"], 4),
                "f1": round(s["f1"], 4),
            }
        return out

    fam = table(per_family)
    bands = table(per_band)
    macro = sum(v["f1"] for v in fam.values()) / len(fam) if fam else float("nan")
    micro = prf(sum(r["tp"] for r in rows), sum(r["fp"] for r in rows), sum(r["fn"] for r in rows))
    pred_counts = [r["n_pred"] for r in rows]
    true_counts = [r["n_true"] for r in rows]

    report = {
        "status": "pipeline exercise, not a protocol number (Amendment 4 A4.5)",
        "why": "the frozen split files exist only on the cluster; this runs on the local ArchiveII mirror",
        "tool": "RNAfold",
        "vienna_version": RNA.__version__,
        "inputs": {
            "archive": args.archive,
            "archive_sha256": sha256_file(path),
            "n_records_in_file": len(data),
            "n_scored": len(rows),
            "n_skipped": n_skipped,
        },
        "summary": {
            "family_macro_f1": round(float(macro), 4),
            "micro": {
                "precision": round(micro["precision"], 4),
                "recall": round(micro["recall"], 4),
                "f1": round(micro["f1"], 4),
                "tp": micro["tp"],
                "fp": micro["fp"],
                "fn": micro["fn"],
            },
            "mean_pairs_pred": round(sum(pred_counts) / len(pred_counts), 3) if pred_counts else None,
            "mean_pairs_true": round(sum(true_counts) / len(true_counts), 3) if true_counts else None,
        },
        "per_family": fam,
        "per_length_band": bands,
        "caveats": [
            "ArchiveII is contaminated by design (61.2% exact duplicates of the BPfold training union), so this number is not a generalisation estimate.",
            "The mirror has 3,370 records while the cmscan input used in E0 had 3,966; the difference is unresolved.",
            "Canonical/GU filtering is not applied here: labels are scored as distributed.",
        ],
    }
    out = ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2) + "\n")
    if args.save_db:
        db_path = ROOT / args.save_db
        db_path.parent.mkdir(parents=True, exist_ok=True)
        db_path.write_text("".join(json.dumps(r) + "\n" for r in db_rows))
    print(json.dumps(report["summary"], indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
