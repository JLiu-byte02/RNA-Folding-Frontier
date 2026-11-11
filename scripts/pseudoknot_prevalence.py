#!/usr/bin/env python3
"""How often do the label sources actually contain crossing pairs?

Amendment 6 needs a measured number before it freezes the primary-F1 transform:
if crossing pairs are rare, the nested reduction changes little; if they are
common in one source, the primary and secondary F1 will differ and the paper has
to say so.

Sources, both reachable locally without the cluster:

* the pinned bpRNA-1m mirror (``data/raw/bprna/data.json``);
* Rfam 15.0 SEED alignments (``data/raw/rfam/Rfam.seed.gz``), projected onto each
  seed sequence exactly as ``scripts/build_rfam_index.py`` does.

The Rfam projection is imported from the builder rather than re-written, so the
number describes the labels the split builder would actually produce.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rff.analysis.labels import pseudoknot_pairs  # noqa: E402
from rff.data.io import pairs_from_dbn  # noqa: E402

_spec = importlib.util.spec_from_file_location(
    "build_rfam_index", ROOT / "scripts" / "build_rfam_index.py"
)
rfam_builder = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(rfam_builder)

CANONICAL = {("A", "U"), ("U", "A"), ("G", "C"), ("C", "G"), ("G", "U"), ("U", "G")}


def summarise(counts: dict[str, list[int]]) -> dict:
    total_pairs = sum(v[0] for v in counts.values())
    total_crossing = sum(v[1] for v in counts.values())
    with_cross = sum(1 for v in counts.values() if v[1] > 0)
    return {
        "n_sequences": len(counts),
        "n_sequences_with_a_crossing_pair": with_cross,
        "fraction_sequences_with_crossing": round(with_cross / len(counts), 6) if counts else None,
        "n_pairs": total_pairs,
        "n_crossing_pairs": total_crossing,
        "fraction_pairs_crossing": round(total_crossing / total_pairs, 6) if total_pairs else None,
    }


def bprna(path: Path) -> dict[str, list[int]]:
    data = json.loads(path.read_text(encoding="utf-8"))
    out: dict[str, list[int]] = {}
    for rid, rec in data.items():
        seq = rec.get("sequence") or ""
        pairs = {(int(i), int(j)) for i, j in rec.get("structure", []) if i != j}
        canonical = {(i, j) for i, j in pairs if 0 <= i < len(seq) and 0 <= j < len(seq)
                     and (seq[i], seq[j]) in CANONICAL}
        out[rid] = [len(canonical), len(pseudoknot_pairs(canonical))]
    return out


def rfam(path: Path, max_per_family: int) -> dict[str, list[int]]:
    out: dict[str, list[int]] = {}
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
                if not ungapped or len(ungapped) != len(
                    rfam_builder.project_ss(ss_cons, aligned)
                ):
                    continue
                dbn = rfam_builder.project_ss(ss_cons, aligned)
                try:
                    raw = pairs_from_dbn(dbn)
                except ValueError:
                    continue
                canonical = {(i, j) for i, j in raw if (ungapped[i], ungapped[j]) in CANONICAL}
                out[f"{acc}|{seq_name}"] = [len(canonical), len(pseudoknot_pairs(canonical))]
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bprna", default="data/raw/bprna/data.json")
    ap.add_argument("--rfam", default="data/raw/rfam/Rfam.seed.gz")
    ap.add_argument("--max-per-family", type=int, default=50)
    ap.add_argument("--out", default="experiments/005_label_policy/pseudoknot_prevalence.json")
    args = ap.parse_args()

    report = {"status": "label-source measurement, not a model result", "sources": {}}
    if (ROOT / args.bprna).exists():
        bp = bprna(ROOT / args.bprna)
        report["sources"]["bprna"] = summarise(bp)
    if (ROOT / args.rfam).exists():
        rf = rfam(ROOT / args.rfam, args.max_per_family)
        report["sources"]["rfam_seed"] = summarise(rf)
    report["note"] = (
        "crossing pairs are counted after the canonical/GU filter, which is the "
        "state in which labels reach the scorer; Rfam seeds are capped per family "
        "for cost, and the cap is recorded"
    )
    report["rfam_cap_per_family"] = args.max_per_family

    out = ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report["sources"], indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
