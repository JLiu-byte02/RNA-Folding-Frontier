#!/usr/bin/env python3
"""Build pair-level projected Rfam seed labels (``labels_v3``).

PI Round 4h §3. The v2 artifact projects the consensus *characters* and then
parses the projected string. When a gap removes one endpoint of a helix, that
can silently pair the surviving endpoint with a different closing column. v3
parses ``SS_cons`` in alignment-column coordinates, maps both endpoints to
ungapped sequence coordinates, drops the pair if either column is a gap, and
then applies the frozen canonical-pair filter.

The output is a full replacement map for every seed row, not a delta against a
delta. A missing key would otherwise be ambiguous between "unchanged" and
"parse failed"; values are either a canonical pair list or ``null`` for a
consensus row that could not be parsed unambiguously.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.util
import json
import re
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rff.data.io import pairs_from_dbn  # noqa: E402
from rff.data.wuss import (  # noqa: E402
    canonical_pairs,
    pairs_from_wuss,
    project_pairs,
    project_wuss,
)

_spec = importlib.util.spec_from_file_location(
    "build_rfam_index", ROOT / "scripts" / "build_rfam_index.py"
)
rfam_builder = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(rfam_builder)


def bucket(n: int) -> str:
    if n == 0:
        return "0"
    if n < 4:
        return "1-3"
    if n < 10:
        return "4-9"
    if n < 20:
        return "10-19"
    return "20+"


def read_eval(path: Path) -> dict[str, dict]:
    out: dict[str, dict] = {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if not line.strip():
                continue
            rec = json.loads(line)
            key = str(rec["id"])
            if key.startswith("rfam:"):
                key = key[len("rfam:"):]
            out[key] = rec
    return out


def read_split_dir(path: Path) -> dict[str, dict[str, dict]]:
    """Read membership files, normalising ``rfam:`` prefixes on ids."""
    out: dict[str, dict[str, dict]] = {}
    for split_path in sorted(path.iterdir()):
        if split_path.suffix not in {".jsonl", ".json", ".txt", ".tsv"}:
            continue
        if split_path.name.startswith("manifest"):
            continue
        records: dict[str, dict] = {}
        with open(split_path, encoding="utf-8") as fh:
            if split_path.suffix == ".json":
                payload = json.load(fh)
                if isinstance(payload, list):
                    items = [(None, rec) for rec in payload if isinstance(rec, dict)]
                elif isinstance(payload, dict):
                    items = list(payload.items())
                else:
                    raise ValueError(f"{split_path}: unsupported JSON topology")
            else:
                items = []
                for line in fh:
                    line = line.strip()
                    if not line:
                        continue
                    if line.startswith("{"):
                        items.append((None, json.loads(line)))
                    else:
                        items.append((line.split("\t")[0], {}))
            for fallback_key, rec in items:
                if not isinstance(rec, dict):
                    rec = {}
                value = rec.get("id") or rec.get("name") or fallback_key
                if not value:
                    continue
                key = str(value)
                if key.startswith("rfam:"):
                    key = key[len("rfam:"):]
                records[key] = rec
        if records:
            out[split_path.name] = records
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", default="data/raw/rfam/Rfam.seed.gz")
    ap.add_argument("--max-per-family", type=int, default=0)
    ap.add_argument("--eval-jsonl", default="experiments/001_e1prime/test_novel_eval.jsonl")
    ap.add_argument("--split-dir", default="", help="optional: per-split membership files")
    ap.add_argument("--out-dir", default="experiments/011_labels_v3")
    args = ap.parse_args()

    out_dir = ROOT / args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    eval_path = ROOT / args.eval_jsonl
    eval_recs = read_eval(eval_path) if eval_path.exists() else {}
    eval_seen: set[str] = set()
    eval_v1_mismatch: list[str] = []
    eval_v1_empty_rescued = 0
    eval_v3_parse_errors = 0

    pairs_out: dict[str, list[list[int]] | None] = {}
    n_seqs = 0
    n_length_mismatch = 0
    n_v1_parse_ok = 0
    n_v1_parse_fail = 0
    n_v2_parse_ok = 0
    n_v2_parse_fail = 0
    n_v3_parse_ok = 0
    n_v3_parse_fail = 0
    n_v2_v3_compared = 0
    n_v2_v3_disagreements = 0
    n_v2_only = 0
    n_v3_only = 0
    n_changed_vs_v2_effective = 0
    pairs_v1_effective = 0
    pairs_v2_effective = 0
    pairs_v3 = 0
    v3_dist = Counter()
    v2_dist = Counter()
    disagreement_examples: list[dict] = []

    with rfam_builder.open_maybe_gzip(ROOT / args.seed) as handle:
        for meta, seqs, ss_cons in rfam_builder.parse_stockholm(handle):
            acc = meta.get("AC", "").split()[0] if meta.get("AC") else ""
            if not acc or not ss_cons:
                continue
            items = list(seqs.items())
            if args.max_per_family:
                items = items[: args.max_per_family]
            for name, aligned in items:
                ungapped = re.sub(r"[-._~]", "", aligned).upper().replace("T", "U")
                if not ungapped:
                    continue
                n_seqs += 1
                key = f"{acc}|{name}"

                if len(aligned) != len(ss_cons):
                    n_length_mismatch += 1
                    p1 = p2 = p3 = None
                else:
                    dbn_v1 = rfam_builder.project_ss(ss_cons, aligned)
                    if len(dbn_v1) != len(ungapped):
                        n_length_mismatch += 1
                        p1 = p2 = p3 = None
                    else:
                        try:
                            p1 = canonical_pairs(ungapped, pairs_from_dbn(dbn_v1))
                            n_v1_parse_ok += 1
                        except ValueError:
                            n_v1_parse_fail += 1
                            p1 = None

                        try:
                            p2 = canonical_pairs(
                                ungapped,
                                pairs_from_wuss(project_wuss(ss_cons, aligned)),
                            )
                            n_v2_parse_ok += 1
                        except ValueError:
                            n_v2_parse_fail += 1
                            p2 = None

                        try:
                            p3 = canonical_pairs(
                                ungapped, project_pairs(ss_cons, aligned)
                            )
                            n_v3_parse_ok += 1
                        except ValueError:
                            n_v3_parse_fail += 1
                            p3 = None

                p1_eff = p1 if p1 is not None else set()
                p2_eff = p2 if p2 is not None else p1_eff
                pairs_v1_effective += len(p1_eff)
                pairs_v2_effective += len(p2_eff)
                v2_dist[bucket(len(p2_eff))] += 1

                if p2 is not None and p3 is not None:
                    n_v2_v3_compared += 1
                    if p2 != p3:
                        n_v2_v3_disagreements += 1
                        n_v2_only += len(p2 - p3)
                        n_v3_only += len(p3 - p2)
                        if len(disagreement_examples) < 20:
                            disagreement_examples.append(
                                {
                                    "key": key,
                                    "v2_only": sorted(p2 - p3),
                                    "v3_only": sorted(p3 - p2),
                                }
                            )

                if p3 is None:
                    pairs_out[key] = None
                else:
                    n_v3_dist = len(p3)
                    pairs_v3 += n_v3_dist
                    v3_dist[bucket(n_v3_dist)] += 1
                    if p3 != p2_eff:
                        n_changed_vs_v2_effective += 1
                    pairs_out[key] = [list(p) for p in sorted(p3)]

                if key in eval_recs:
                    eval_seen.add(key)
                    rec = eval_recs[key]
                    rec_p1 = {(int(i), int(j)) for i, j in rec.get("pairs", [])}
                    if rec_p1 != p1_eff:
                        eval_v1_mismatch.append(key)
                    if not rec_p1 and p3:
                        eval_v1_empty_rescued += 1
                    if p3 is None:
                        eval_v3_parse_errors += 1

    payload = json.dumps(pairs_out, sort_keys=True, separators=(",", ":")).encode("utf-8")
    blob = gzip.compress(payload, mtime=0)
    pairs_path = out_dir / "labels_v3_pairs.json.gz"
    pairs_path.write_bytes(blob)

    eval_audit: dict = {"status": "not requested"} if not eval_path.exists() else {}
    if eval_recs:
        n_eval = len(eval_recs)
        v1_empty = v1_lt4 = v3_ok = v3_empty = v3_lt4 = 0
        v1_dist = Counter()
        v3_dist_eval = Counter()
        for key, rec in eval_recs.items():
            p1 = {(int(i), int(j)) for i, j in rec.get("pairs", [])}
            v1_dist[bucket(len(p1))] += 1
            v1_empty += not p1
            v1_lt4 += len(p1) < 4
            p3 = pairs_out.get(key, "__missing__")
            if p3 == "__missing__":
                continue
            if p3 is None:
                continue
            v3_ok += 1
            v3_dist_eval[bucket(len(p3))] += 1
            v3_empty += not p3
            v3_lt4 += len(p3) < 4
        dropout = v3_lt4 / n_eval if n_eval else 0.0
        eval_audit = {
            "eval_jsonl": str(eval_path.relative_to(ROOT)),
            "n_records": n_eval,
            "n_keys_missing_from_seed": n_eval - len(eval_seen),
            "n_v1_pairs_mismatch_against_reprojection": len(eval_v1_mismatch),
            "v1_pair_count_distribution": dict(sorted(v1_dist.items())),
            "n_v1_empty_truth": v1_empty,
            "n_v1_lt4_pairs": v1_lt4,
            "n_v3_parse_ok": v3_ok,
            "n_v3_parse_errors": eval_v3_parse_errors,
            "v3_pair_count_distribution": dict(sorted(v3_dist_eval.items())),
            "n_v3_empty_truth": v3_empty,
            "n_v3_lt4_pairs": v3_lt4,
            "n_v1_empty_rescued_by_v3": eval_v1_empty_rescued,
            "n_excluded_by_lt4_rule": v3_lt4,
            "excluded_fraction": v3_lt4 / n_eval if n_eval else 0.0,
            "dropout_fraction_vs_current_subset": dropout,
            "lt4_gate_exceeded_10pct": (v3_lt4 / n_eval > 0.10) if n_eval else False,
        }

    summary = {
        "purpose": "pair-level projection of Rfam consensus structure (labels_v3)",
        "seed_source": "data/raw/rfam/Rfam.seed.gz",
        "max_per_family": args.max_per_family,
        "n_sequences": n_seqs,
        "n_length_mismatch": n_length_mismatch,
        "n_v1_parse_ok": n_v1_parse_ok,
        "n_v1_parse_fail": n_v1_parse_fail,
        "n_v2_parse_ok": n_v2_parse_ok,
        "n_v2_parse_fail": n_v2_parse_fail,
        "n_v3_parse_ok": n_v3_parse_ok,
        "n_v3_parse_fail": n_v3_parse_fail,
        "parse_rate_v3": n_v3_parse_ok / max(1, n_seqs),
        "pairs_effective_v1": pairs_v1_effective,
        "pairs_effective_v2": pairs_v2_effective,
        "pairs_v3": pairs_v3,
        "pairs_added_v3_vs_v2_effective": pairs_v3 - pairs_v2_effective,
        "n_changed_vs_v2_effective": n_changed_vs_v2_effective,
        "n_v2_v3_compared": n_v2_v3_compared,
        "n_v2_v3_disagreements": n_v2_v3_disagreements,
        "pairs_v2_only_on_compared_rows": n_v2_only,
        "pairs_v3_only_on_compared_rows": n_v3_only,
        "disagreement_examples": disagreement_examples,
        "v2_effective_pair_count_distribution": dict(sorted(v2_dist.items())),
        "v3_pair_count_distribution": dict(sorted(v3_dist.items())),
        "n_v3_zero_pairs": v3_dist.get("0", 0),
        "n_v3_lt4_pairs": v3_dist.get("0", 0) + v3_dist.get("1-3", 0),
        "eval_subset_audit": eval_audit,
        "pairs_file": str(pairs_path.relative_to(ROOT)),
        "pairs_bytes": len(blob),
        "pairs_bytes_uncompressed": len(payload),
        "pairs_sha256": hashlib.sha256(blob).hexdigest(),
    }
    if args.split_dir:
        splits = read_split_dir(ROOT / args.split_dir)
        per_split = {}
        for name, records in splits.items():
            n = len(records)
            v1_empty = v1_lt4 = v3_ok = v3_empty = v3_lt4 = missing = 0
            v1_known = 0
            for key, rec in records.items():
                p3 = pairs_out.get(key, "__missing__")
                if p3 == "__missing__":
                    missing += 1
                    continue
                if p3 is None:
                    continue
                v3_ok += 1
                v3_empty += not p3
                v3_lt4 += len(p3) < 4
                if "pairs" in rec:
                    v1_known += 1
                    p1 = {(int(i), int(j)) for i, j in rec["pairs"]}
                    v1_empty += not p1
                    v1_lt4 += len(p1) < 4
            per_split[name] = {
                "n_records": n,
                "n_ids_missing_from_seed": missing,
                "n_v1_pairs_known": v1_known,
                "n_v1_empty_truth": v1_empty,
                "n_v1_lt4_pairs": v1_lt4,
                "n_v3_parse_ok": v3_ok,
                "n_v3_parse_errors": n - missing - v3_ok,
                "n_v3_empty_truth": v3_empty,
                "n_v3_lt4_pairs": v3_lt4,
                "n_excluded_by_lt4_rule": v3_lt4,
                "excluded_fraction": v3_lt4 / n if n else 0.0,
            }
        summary["per_split"] = per_split
    else:
        summary["per_split"] = {
            "status": "not computed: pass --split-dir with the server membership lists"
        }
    (out_dir / "labels_v3_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "n_sequences": n_seqs,
                "v3_parse_rate": summary["parse_rate_v3"],
                "pairs_v2_effective": pairs_v2_effective,
                "pairs_v3": pairs_v3,
                "v2_v3_disagreements": n_v2_v3_disagreements,
                "v3_lt4_fraction": summary["n_v3_lt4_pairs"] / max(1, n_seqs),
                "eval": eval_audit,
                "pairs_bytes": len(blob),
                "pairs_sha256": summary["pairs_sha256"],
            }
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
