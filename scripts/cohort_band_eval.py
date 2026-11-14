#!/usr/bin/env python3
"""Length-band stratified version of ``cohort_paired_eval.py``.

The common cohort caps at 597 nt, and families differ wildly in length, so a
single family-macro number can hide a length-dependent reversal.  This script
recomputes per-sequence F1 under the same truth policy as the cohort evaluation
and aggregates family-macro F1 inside fixed length bands, with a family
bootstrap inside each band for the paired deltas.

Descriptive only: bands are chosen before looking at the numbers and no band is
used to select a model or a threshold.

Usage mirrors cohort_paired_eval:

    python scripts/cohort_band_eval.py --split ... --cohort ... \
        --method RNAfold=... --method head_mean=... --compare head_mean:RNAfold \
        --out band.json
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


def load_preds(path: Path) -> dict[str, set[tuple[int, int]]]:
    out: dict[str, set[tuple[int, int]]] = {}
    for row in seb.read_jsonl(path):
        if "pairs" in row:
            out[row["id"]] = {(int(a), int(b)) for a, b in row["pairs"]}
        elif "dbn" in row:
            from rff.data.io import pairs_from_dbn

            out[row["id"]] = set(pairs_from_dbn(row["dbn"]))
        else:
            raise SystemExit(f"{path}: row {row.get('id')!r} has neither pairs nor dbn")
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--split", required=True)
    ap.add_argument("--cohort", required=True)
    ap.add_argument("--labels", default=seb.LABELS_V3)
    ap.add_argument("--method", action="append", required=True, help="NAME=PATH")
    ap.add_argument("--mean-of", action="append", default=[])
    ap.add_argument("--compare", action="append", default=[], help="A:B")
    ap.add_argument("--bands", default="0,80,150,300,600")
    ap.add_argument("--min-true-pairs", type=int, default=4)
    ap.add_argument("--n-boot", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=20261113)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    split_path = _resolve(args.split)
    cohort_path = _resolve(args.cohort)
    seb.refuse_embargoed(split_path)
    seb.refuse_embargoed(cohort_path)

    split_rows = {r["id"]: r for r in seb.read_jsonl(split_path)}
    cohort = seb.read_jsonl(cohort_path)
    labels = seb.load_labels_v3(_resolve(args.labels))

    preds = {}
    for spec in args.method:
        name, _, path = spec.partition("=")
        preds[name] = load_preds(_resolve(path))

    kept, family, length = [], {}, {}
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
        if len(raw) < args.min_true_pairs:
            continue
        if any(rid not in preds[name] for name in preds):
            continue
        kept.append(rid)
        family[rid] = row["family"] if row.get("family") is not None else src["family"]
        length[rid] = int(row["length"]) if row.get("length") is not None else len(src["seq"])

    per_seq: dict[str, dict[str, float]] = {}
    for name, table in preds.items():
        f1s = {}
        for rid in kept:
            src = split_rows.get(rid)
            key = rid[5:] if rid.startswith("rfam:") else rid
            raw = labels[key] if key in labels else {(int(a), int(b)) for a, b in src["pairs"]}
            truth = max_nested_subset(raw)
            pred = table[rid]
            t = len(truth & pred)
            f1s[rid] = seb.f1_from_counts(t, len(pred - truth), len(truth - pred))
        per_seq[name] = f1s

    for spec in args.mean_of:
        name, _, parts = spec.partition("=")
        plist = [p for p in parts.split(",") if p]
        per_seq[name] = {rid: float(np.mean([per_seq[p][rid] for p in plist])) for rid in kept}

    edges = [int(x) for x in args.bands.split(",")]
    bands = list(zip(edges[:-1], edges[1:]))
    rng = np.random.default_rng(args.seed)
    result = {"split": args.split, "cohort": args.cohort, "bands": [], "n_sequences": len(kept)}

    for lo, hi in bands:
        ids = [rid for rid in kept if lo <= length[rid] < hi]
        if not ids:
            result["bands"].append({"band": [lo, hi], "n_sequences": 0})
            continue
        fams = sorted({family[rid] for rid in ids})
        fam_idx = {f: k for k, f in enumerate(fams)}
        acc = {name: np.zeros(len(fams)) for name in per_seq}
        cnt = np.zeros(len(fams))
        for rid in ids:
            k = fam_idx[family[rid]]
            cnt[k] += 1
            for name in per_seq:
                acc[name][k] += per_seq[name][rid]
        fam_means = {name: acc[name] / cnt for name in per_seq}
        idx = rng.integers(0, len(fams), size=(args.n_boot, len(fams)))
        draws = {name: fam_means[name][idx].mean(axis=1) for name in fam_means}
        methods = {}
        for name, means in fam_means.items():
            lo_ci, hi_ci = np.percentile(draws[name], [2.5, 97.5])
            methods[name] = {
                "family_macro_f1": round(float(means.mean()), 4),
                "ci95": [round(float(lo_ci), 4), round(float(hi_ci), 4)],
            }
        paired = []
        for spec in args.compare:
            a, _, b = spec.partition(":")
            d = draws[a] - draws[b]
            lo_ci, hi_ci = np.percentile(d, [2.5, 97.5])
            paired.append({
                "a": a, "b": b,
                "delta": round(float(fam_means[a].mean() - fam_means[b].mean()), 4),
                "ci95": [round(float(lo_ci), 4), round(float(hi_ci), 4)],
                "frac_families_a_better": round(float((fam_means[a] > fam_means[b]).mean()), 3),
            })
        result["bands"].append({
            "band": [lo, hi], "n_sequences": len(ids), "n_families": len(fams),
            "methods": methods, "paired": paired,
        })
        print(f"band {lo}-{hi}: {len(ids)} seqs / {len(fams)} families")
        for name, m in methods.items():
            print(f"  {name:12s} {m['family_macro_f1']:.4f} [{m['ci95'][0]:.4f},{m['ci95'][1]:.4f}]")
        for p in paired:
            print(f"  {p['a']} - {p['b']}: {p['delta']:+.4f} [{p['ci95'][0]:+.4f},{p['ci95'][1]:+.4f}]")

    out = _resolve(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
