#!/usr/bin/env python3
"""Compare bpRNA and ArchiveII annotations of identical sequences.

For every ArchiveII sequence that appears verbatim in bpRNA TR0, compare the two
dot-bracket / bpseq annotations at the base-pair level. The resulting agreement
F1 is an empirical ceiling on any model's F1 for those sequences, and the
disagreement breakdown (shifted register vs. genuinely different pairs) says
whether the noise is annotation drift or structure disagreement.

Usage:
    python scripts/label_discordance.py --out experiments/000_e0_leakage/label_discordance
"""

from __future__ import annotations

import argparse
import json
import statistics
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rff.analysis.labels import classify_disagreement, pair_prf  # noqa: E402
from rff.data.io import load_archiveii_json, load_path  # noqa: E402


def load_archiveii(path: Path):
    """Load ArchiveII from the HF JSON dump or from a directory of .bpseq files."""
    if path.is_dir():
        records = load_path(path)
        for rec in records:
            rec.family = rec.rid.split("_", 1)[0]
        return records
    return load_archiveii_json(path)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--archiveii", default="data/raw/archiveII/data.json")
    ap.add_argument(
        "--bpseq-dir", default="data/raw/bpfold/BPfold_data/bpRNA/TR0"
    )
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    bp = load_path(ROOT / args.bpseq_dir)
    print(f"bpRNA records: {len(bp)}", flush=True)
    index: dict[str, list] = defaultdict(list)
    for rec in bp:
        index[rec.seq].append(rec)

    arch = load_archiveii(ROOT / args.archiveii)
    print(f"ArchiveII records: {len(arch)}", flush=True)

    rows: list[dict] = []
    per_family: dict[str, list[float]] = defaultdict(list)
    for rec in arch:
        matches = index.get(rec.seq)
        if not matches:
            continue
        scored = []
        for cand in matches:
            m = pair_prf(rec.pairs, cand.pairs)
            scored.append((m["f1"], cand, m))
        scored.sort(key=lambda t: -t[0])
        best_f1, best, best_m = scored[0]
        cats = classify_disagreement(rec.pairs, best.pairs)
        rows.append(
            {
                "archiveii_id": rec.rid,
                "family": rec.family,
                "length": rec.length,
                "n_bprna_matches": len(matches),
                "best_bprna_id": best.rid,
                "archiveii_pairs": len(rec.pairs),
                "bprna_pairs": len(best.pairs),
                "pair_f1": round(best_f1, 4),
                "precision": round(best_m["precision"], 4),
                "recall": round(best_m["recall"], 4),
                "shared_pairs": cats["shared"],
                "shifted_pairs": cats["shifted"],
                "archiveii_only_pairs": cats["a_only"],
                "bprna_only_pairs": cats["b_only"],
                "identical": best_f1 == 1.0,
            }
        )
        per_family[rec.family or "unknown"].append(best_f1)

    f1s = [r["pair_f1"] for r in rows]
    summary = {
        "n_archiveii": len(arch),
        "n_exact_duplicate_sequences": len(rows),
        "n_identical_annotations": sum(1 for r in rows if r["identical"]),
        "identity_rate": round(
            sum(1 for r in rows if r["identical"]) / len(rows), 4
        )
        if rows
        else 0.0,
        "mean_pair_f1": round(statistics.mean(f1s), 4) if f1s else 0.0,
        "median_pair_f1": round(statistics.median(f1s), 4) if f1s else 0.0,
        "p10_pair_f1": round(sorted(f1s)[int(0.1 * len(f1s))], 4) if f1s else 0.0,
        "n_f1_below_0.8": sum(1 for f in f1s if f < 0.8),
        "n_f1_below_0.5": sum(1 for f in f1s if f < 0.5),
        "per_family": {
            fam: {
                "n": len(v),
                "mean_pair_f1": round(statistics.mean(v), 4),
                "median_pair_f1": round(statistics.median(v), 4),
                "n_identical": sum(1 for f in v if f == 1.0),
            }
            for fam, v in sorted(per_family.items())
        },
    }

    (out_dir / "rows.json").write_text(json.dumps(rows, indent=1))
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2))

    lines = [
        "# Label agreement on identical sequences (ArchiveII vs bpRNA TR0)",
        "",
        f"- ArchiveII records: {summary['n_archiveii']}",
        f"- Exact duplicate sequences: {summary['n_exact_duplicate_sequences']}",
        f"- Identical annotations: {summary['n_identical_annotations']} "
        f"({summary['identity_rate']:.1%})",
        f"- Pair-level F1 (ArchiveII vs bpRNA): mean {summary['mean_pair_f1']}, "
        f"median {summary['median_pair_f1']}, p10 {summary['p10_pair_f1']}",
        f"- Sequences with F1 < 0.8: {summary['n_f1_below_0.8']}; "
        f"F1 < 0.5: {summary['n_f1_below_0.5']}",
        "",
        "| family | n | mean F1 | median F1 | identical |",
        "| --- | --- | --- | --- | --- |",
    ]
    for fam, v in summary["per_family"].items():
        lines.append(
            f"| {fam} | {v['n']} | {v['mean_pair_f1']} | {v['median_pair_f1']} | "
            f"{v['n_identical']} |"
        )
    (out_dir / "report.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
