#!/usr/bin/env python3
"""Within-family test of the short-vs-long crossover.

The band table shows the head at parity below 150 nt and behind above it, but
that is two marginal comparisons. A reviewer will ask whether the *interaction*
is real. This script restricts to families that have at least one sequence in a
short band and one in a long band, computes each family's mean F1 per method in
each band, and bootstraps the difference of paired deltas:

    D = (a_short - b_short) - (a_long - b_long)

over families. A CI excluding zero means the head's relative performance真的
declines with length inside the same families, not because different families
populate the two bands.

Descriptive/robustness analysis; bands are fixed before the run.
"""

from __future__ import annotations

import argparse
import json
import sys
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
    out = {}
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
    ap.add_argument("--short", default="80,150")
    ap.add_argument("--long", default="150,600")
    ap.add_argument("--min-true-pairs", type=int, default=4)
    ap.add_argument("--n-boot", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=20261113)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    short_lo, short_hi = (int(x) for x in args.short.split(","))
    long_lo, long_hi = (int(x) for x in args.long.split(","))

    split_rows = {r["id"]: r for r in seb.read_jsonl(_resolve(args.split))}
    cohort = seb.read_jsonl(_resolve(args.cohort))
    labels = seb.load_labels_v3(_resolve(args.labels))
    preds = {}
    for spec in args.method:
        name, _, path = spec.partition("=")
        preds[name] = load_preds(_resolve(path))

    kept, family, length = [], {}, {}
    truth = {}
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
        if len(raw) < args.min_true_pairs or any(rid not in preds[n] for n in preds):
            continue
        kept.append(rid)
        family[rid] = row["family"] if row.get("family") is not None else src["family"]
        length[rid] = int(row["length"]) if row.get("length") is not None else len(src["seq"])
        truth[rid] = max_nested_subset(raw)

    f1 = {name: {} for name in preds}
    for name, table in preds.items():
        for rid in kept:
            pred = table[rid]
            t = len(truth[rid] & pred)
            f1[name][rid] = seb.f1_from_counts(t, len(pred - truth[rid]), len(truth[rid] - pred))
    for spec in args.mean_of:
        name, _, parts = spec.partition("=")
        plist = [p for p in parts.split(",") if p]
        f1[name] = {rid: float(np.mean([f1[p][rid] for p in plist])) for rid in kept}

    short_ids = [r for r in kept if short_lo <= length[r] < short_hi]
    long_ids = [r for r in kept if long_lo <= length[r] < long_hi]
    fam_short = {}
    fam_long = {}
    for rid in short_ids:
        fam_short.setdefault(family[rid], []).append(rid)
    for rid in long_ids:
        fam_long.setdefault(family[rid], []).append(rid)
    common = sorted(set(fam_short) & set(fam_long))

    result = {
        "split": args.split,
        "short_band": [short_lo, short_hi],
        "long_band": [long_lo, long_hi],
        "n_sequences_short": len(short_ids),
        "n_sequences_long": len(long_ids),
        "n_families_short": len(fam_short),
        "n_families_long": len(fam_long),
        "n_families_common": len(common),
        "n_boot": args.n_boot,
        "seed": args.seed,
        "interactions": [],
    }
    if not common:
        result["note"] = "no family has sequences in both bands"
        _resolve(args.out).write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
        print(json.dumps(result, indent=2))
        return 0

    # family-level band means
    band_mean = {"short": {}, "long": {}}
    for band, fam_map, ids in (("short", fam_short, short_ids), ("long", fam_long, long_ids)):
        for name in f1:
            band_mean[band][name] = {
                fam: float(np.mean([f1[name][rid] for rid in fam_map[fam]])) for fam in common
            }

    rng = np.random.default_rng(args.seed)
    idx = rng.integers(0, len(common), size=(args.n_boot, len(common)))
    for spec in args.compare:
        a, _, b = spec.partition(":")
        d_short = np.array([band_mean["short"][a][f] - band_mean["short"][b][f] for f in common])
        d_long = np.array([band_mean["long"][a][f] - band_mean["long"][b][f] for f in common])
        D = d_short - d_long
        draws = D[idx].mean(axis=1)
        lo, hi = np.percentile(draws, [2.5, 97.5])
        p_lo = float((draws <= 0).mean())
        p_hi = float((draws >= 0).mean())
        entry = {
            "a": a, "b": b,
            "delta_short": round(float(d_short.mean()), 4),
            "delta_long": round(float(d_long.mean()), 4),
            "interaction": round(float(D.mean()), 4),
            "ci95": [round(float(lo), 4), round(float(hi), 4)],
            "p_two_sided": round(min(1.0, 2.0 * min(p_lo, p_hi)), 4),
            "frac_families_declining": round(float((D < 0).mean()), 3),
        }
        result["interactions"].append(entry)
        print(f"{a} - {b}: short {entry['delta_short']:+.3f} long {entry['delta_long']:+.3f} "
              f"interaction {entry['interaction']:+.3f} [{entry['ci95'][0]:+.3f},{entry['ci95'][1]:+.3f}] "
              f"p={entry['p_two_sided']:.4f}")

    out = _resolve(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
