#!/usr/bin/env python3
"""Paired family-block evaluation of several predictors on one explicit cohort.

Two defects this closes:

1. Marginal intervals do not compare.  ``score_external_baseline.py`` reports
   each method's family-macro F1 with its own bootstrap interval; overlapping
   intervals are routinely (and wrongly) read as "no difference".  Here the
   difference A - B is bootstrapped by resampling *families*, using the same
   resample for both members, so the interval is on the quantity of interest.

2. Cohort drift.  Every method is scored on the id list written by
   ``build_common_cohort.py`` and nothing else; a missing prediction is an
   error here, not a silent skip.

Label policy is imported from ``score_external_baseline`` (one source of truth):
rfam rows use ``labels_v3``, bpRNA rows use their frozen ``pairs``;
``--min-true-pairs`` applies to the raw set before the nested projection.

Outputs family-level mean F1 per method (so the bootstrap can be re-derived
without the predictions) plus the paired deltas requested with ``--compare``.
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
    ap.add_argument("--mean-of", action="append", default=[], help="NAME=part1,part2")
    ap.add_argument("--compare", action="append", default=[], help="A:B")
    ap.add_argument("--min-true-pairs", type=int, default=4)
    ap.add_argument("--n-boot", type=int, default=10000)
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

    paths: dict[str, Path] = {}
    for spec in args.method:
        name, _, path = spec.partition("=")
        paths[name] = _resolve(path)

    preds = {name: load_preds(path) for name, path in paths.items()}

    kept: list[str] = []
    family: dict[str, str] = {}
    n_lt_min = n_no_pred = n_parsefail = 0
    truth_cache: dict[str, set[tuple[int, int]]] = {}
    for row in cohort:
        rid = row["id"]
        src = split_rows.get(rid)
        key = rid[5:] if rid.startswith("rfam:") else rid
        if key in labels:
            raw = labels[key]
            if raw is None:
                n_parsefail += 1
                continue
        elif src is not None:
            raw = {(int(a), int(b)) for a, b in src["pairs"]}
        else:
            raise SystemExit(f"{rid}: no label and not in split")
        if len(raw) < args.min_true_pairs:
            n_lt_min += 1
            continue
        if any(rid not in preds[name] for name in preds):
            n_no_pred += 1
            continue
        truth_cache[rid] = max_nested_subset(raw)
        family[rid] = row["family"] if row.get("family") is not None else src["family"]
        kept.append(rid)

    if not kept:
        raise SystemExit("empty cohort after filtering")
    if n_no_pred:
        raise SystemExit(f"{n_no_pred} cohort rows lack at least one prediction; cohort is not common")

    fams = sorted(set(family.values()))
    fam_index = {f: k for k, f in enumerate(fams)}

    per_seq: dict[str, dict[str, float]] = {}
    totals: dict[str, list[int]] = {}
    for name, table in preds.items():
        f1s: dict[str, float] = {}
        tp = fp = fn = 0
        for rid in kept:
            truth = truth_cache[rid]
            pred = table[rid]
            t = len(truth & pred)
            a = len(pred - truth)
            b = len(truth - pred)
            tp += t
            fp += a
            fn += b
            f1s[rid] = seb.f1_from_counts(t, a, b)
        per_seq[name] = f1s
        totals[name] = [tp, fp, fn]

    for spec in args.mean_of:
        name, _, parts = spec.partition("=")
        plist = [p for p in parts.split(",") if p]
        missing = [p for p in plist if p not in per_seq]
        if missing:
            raise SystemExit(f"--mean-of {name}: unknown parts {missing}")
        per_seq[name] = {
            rid: float(np.mean([per_seq[p][rid] for p in plist])) for rid in kept
        }
        totals[name] = [sum(totals[p][k] for p in plist) for k in range(3)]

    fam_means: dict[str, np.ndarray] = {}
    for name, f1s in per_seq.items():
        acc = np.zeros(len(fams), dtype=float)
        cnt = np.zeros(len(fams), dtype=float)
        for rid in kept:
            k = fam_index[family[rid]]
            acc[k] += f1s[rid]
            cnt[k] += 1.0
        fam_means[name] = acc / cnt

    rng = np.random.default_rng(args.seed)
    idx = rng.integers(0, len(fams), size=(args.n_boot, len(fams)))

    methods: dict[str, dict] = {}
    draws: dict[str, np.ndarray] = {}
    for name, means in fam_means.items():
        draws[name] = means[idx].mean(axis=1)
        lo, hi = np.percentile(draws[name], [2.5, 97.5])
        tp, fp, fn = totals[name]
        methods[name] = {
            "family_macro_f1": round(float(means.mean()), 6),
            "family_macro_f1_ci95": [round(float(lo), 6), round(float(hi), 6)],
            "micro_f1": round(seb.f1_from_counts(tp, fp, fn), 6),
            "tp": tp,
            "fp": fp,
            "fn": fn,
            "per_family_mean_f1": {f: round(float(means[k]), 6) for k, f in enumerate(fams)},
        }

    paired = []
    for spec in args.compare:
        a, _, b = spec.partition(":")
        if a not in fam_means or b not in fam_means:
            raise SystemExit(f"--compare {spec}: unknown method")
        delta_draws = draws[a] - draws[b]
        lo, hi = np.percentile(delta_draws, [2.5, 97.5])
        p_lo = float((delta_draws <= 0).mean())
        p_hi = float((delta_draws >= 0).mean())
        diff = fam_means[a] - fam_means[b]
        paired.append(
            {
                "a": a,
                "b": b,
                "delta": round(float(fam_means[a].mean() - fam_means[b].mean()), 6),
                "ci95": [round(float(lo), 6), round(float(hi), 6)],
                "p_two_sided": round(min(1.0, 2.0 * min(p_lo, p_hi)), 6),
                "frac_families_a_better": round(float((diff > 0).mean()), 4),
                "frac_families_tied": round(float((diff == 0).mean()), 4),
            }
        )

    result = {
        "split": args.split,
        "cohort": args.cohort,
        "label_set": args.labels,
        "min_true_pairs": args.min_true_pairs,
        "n_sequences": len(kept),
        "n_families": len(fams),
        "n_boot": args.n_boot,
        "bootstrap_seed": args.seed,
        "excluded": {
            "lt_min_true_pairs": n_lt_min,
            "parse_failure": n_parsefail,
            "missing_prediction": n_no_pred,
        },
        "methods": methods,
        "paired": paired,
    }
    out = _resolve(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print(f"cohort {len(kept)} seqs / {len(fams)} families")
    for name in fam_means:
        m = methods[name]
        print(
            f"  {name:12s} macro {m['family_macro_f1']:.4f} "
            f"[{m['family_macro_f1_ci95'][0]:.4f},{m['family_macro_f1_ci95'][1]:.4f}] "
            f"micro {m['micro_f1']:.4f}"
        )
    for p in paired:
        print(
            f"  {p['a']} - {p['b']}: {p['delta']:+.4f} "
            f"[{p['ci95'][0]:+.4f},{p['ci95'][1]:+.4f}] p={p['p_two_sided']:.4f} "
            f"fam-win {p['frac_families_a_better']:.3f}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
