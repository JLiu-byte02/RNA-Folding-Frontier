"""Export a split JSONL to BPSEQ + .lst for MXfold2 retraining (A13.1 control).

MXfold2's folding layer is a nested dynamic program, so crossing pairs cannot
be represented.  We therefore project each reference structure onto its maximum
nested subset with the same polynomial solver used everywhere else in this
project (``rff.analysis.labels.max_nested_subset``) and record exactly how many
pairs that projection removed.  No other label transformation is applied.

The export also writes a manifest that pins the input split hash, the cap, and
the output list hash so a later reader can tell which sequences MXfold2 saw.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rff.analysis.labels import max_nested_subset  # noqa: E402


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def safe_name(idx: int) -> str:
    return f"seq{idx:07d}"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--split", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--max-len", type=int, default=600)
    ap.add_argument("--min-pairs", type=int, default=1)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--stride", type=int, default=1,
                    help="keep every Nth row (before the length cap); spread a subset over the split")
    ap.add_argument("--manifest", default=None)
    args = ap.parse_args()

    split = Path(args.split)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    with open(split, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if line:
                rows.append(json.loads(line))

    exported = []
    mapping = []
    n_too_long = 0
    n_few_pairs = 0
    n_crossing_rows = 0
    pairs_before = 0
    pairs_after = 0
    for idx, row in enumerate(rows):
        if args.limit and len(exported) >= args.limit:
            break
        if args.stride > 1 and idx % args.stride:
            continue
        seq = row["seq"]
        if len(seq) > args.max_len:
            n_too_long += 1
            continue
        pairs = {tuple(sorted(p)) for p in row.get("pairs", [])}
        pairs_before += len(pairs)
        nested = max_nested_subset(pairs)
        pairs_after += len(nested)
        if len(nested) != len(pairs):
            n_crossing_rows += 1
        if len(nested) < args.min_pairs:
            n_few_pairs += 1
            continue
        name = safe_name(idx)
        path = out_dir / f"{name}.bpseq"
        partner = {i: j for i, j in nested}
        partner.update({j: i for i, j in nested})
        with open(path, "w", encoding="utf-8", newline="\n") as out:
            for k, base in enumerate(seq, start=1):
                j = partner.get(k - 1, 0)
                out.write(f"{k} {base} {j + 1 if j else 0}\n")
        exported.append(path)
        mapping.append({"id": row["id"], "file": path.name, "family": row.get("family"),
                        "clan": row.get("clan"), "length": len(seq),
                        "n_pairs_full": len(pairs), "n_pairs_nested": len(nested)})

    lst = out_dir / "train.lst"
    with open(lst, "w", encoding="utf-8", newline="\n") as out:
        for path in exported:
            out.write(f"{path}\n")
    with open(out_dir / "mapping.jsonl", "w", encoding="utf-8", newline="\n") as out:
        for rec in mapping:
            out.write(json.dumps(rec, sort_keys=True) + "\n")

    manifest = {
        "split": str(split),
        "split_sha256": sha256_file(split),
        "n_rows_in_split": len(rows),
        "n_exported": len(exported),
        "n_skipped_too_long": n_too_long,
        "n_skipped_few_pairs": n_few_pairs,
        "max_len": args.max_len,
        "min_pairs": args.min_pairs,
        "limit": args.limit,
        "stride": args.stride,
        "n_rows_with_crossing_pairs": n_crossing_rows,
        "pairs_before_nested_projection": pairs_before,
        "pairs_after_nested_projection": pairs_after,
        "lst_sha256": sha256_file(lst),
        "exported_len_mean": round(sum(m["length"] for m in mapping) / len(mapping), 2) if mapping else 0.0,
        "exported_len_p50": sorted(m["length"] for m in mapping)[len(mapping) // 2] if mapping else 0,
        "exported_len_max": max((m["length"] for m in mapping), default=0),
        "label_policy": "frozen split pairs projected to the maximum nested "
                        "subset; no sequence-level filtering beyond length and "
                        "min_pairs",
    }
    dest = Path(args.manifest) if args.manifest else out_dir / "manifest.json"
    dest.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
