#!/usr/bin/env python3
"""How much of each sequence is forced at +-w around the sequence centre?

The swap conditions (b/d/e) only rewrite the *distal* region, defined as
everything outside ``[L//2 - w, L//2 + w]``.  For a sequence shorter than
``2w+1`` nothing is distal, so the intervention is vacuous and condition (x)
must be bit-identical to (a).  This script quantifies that, per split, so the
P2 table can be stratified by how much of the input was actually perturbed.
"""

from __future__ import annotations

import argparse
import json
import statistics
from pathlib import Path


def n_distal(length: int, window: int) -> int:
    kept = min(length, length // 2 + window + 1) - max(0, length // 2 - window)
    return length - max(0, kept)


def band(length: int) -> str:
    for upper, name in ((100, "L<100"), (200, "100-199"), (400, "200-399"), (10**9, "L>=400")):
        if length < upper:
            return name
    return "L>=400"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--splits", nargs="+", required=True)
    ap.add_argument("--window", type=int, default=32)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    report = {}
    for path in args.splits:
        name = Path(path).stem
        lengths, fracs, bands = [], [], {}
        for line in open(path):
            if not line.strip():
                continue
            rec = json.loads(line)
            L = len(rec["seq"])
            d = n_distal(L, args.window)
            lengths.append(L)
            fracs.append(d / L)
            b = bands.setdefault(band(L), {"n": 0, "vacuous": 0})
            b["n"] += 1
            b["vacuous"] += int(d == 0)
        fracs.sort()
        report[name] = {
            "n": len(lengths),
            "window": args.window,
            "len_mean": round(statistics.fmean(lengths), 1),
            "len_median": statistics.median(lengths),
            "len_max": max(lengths),
            "vacuous": sum(1 for f in fracs if f == 0),
            "vacuous_frac": round(sum(1 for f in fracs if f == 0) / len(fracs), 4),
            "distal_frac_median": round(statistics.median(fracs), 4),
            "distal_frac_p90": round(fracs[int(0.9 * (len(fracs) - 1))], 4),
            "distal_frac_mean": round(statistics.fmean(fracs), 4),
            "bands": {k: {**v, "vacuous_frac": round(v["vacuous"] / v["n"], 4)} for k, v in sorted(bands.items())},
        }
        r = report[name]
        print(
            f"{name:16s} n={r['n']:6d} L_med={r['len_median']:6.0f} "
            f"vacuous={r['vacuous']:5d} ({r['vacuous_frac']:.1%}) "
            f"distal_frac med={r['distal_frac_median']:.2f} p90={r['distal_frac_p90']:.2f}"
        )

    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(json.dumps(report, indent=2) + "\n")
        print(f"wrote {args.out}")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
