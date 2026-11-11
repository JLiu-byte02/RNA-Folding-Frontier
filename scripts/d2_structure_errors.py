#!/usr/bin/env python3
"""D2 — are the false positives "spurious partners"? (H2, predictions 1-3)

H2 says the distal region supplies sequence that is locally complementary to the
kept region, the 2D head believes it, and the result is a confident long-range
pair that no global constraint ever rejects. Its three checkable consequences:

1. under native inputs, long-range false positives grow with sequence length;
2. false positives are enriched for a stackable complementary run of >= 4 pairs
   that is absent from the true structure;
3. under condition (c), the new false positives link the kept window to a
   distal position that has exactly that complementarity.

Helix length is computed from the sequence itself, so this is a property of the
prediction, not of the model's internals.
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]

PAIRS = set("AU UA GC CG GU UG".split())


def read_jsonl(path: Path) -> list[dict]:
    with open(path) as fh:
        return [json.loads(line) for line in fh if line.strip()]


def helix_length(seq: str, i: int, j: int, cap: int = 40) -> int:
    """Length of the contiguous stack containing pair (i, j)."""
    if seq[i] + seq[j] not in PAIRS:
        return 0
    back = 0
    while (
        back < cap
        and i - back - 1 >= 0
        and j + back + 1 < len(seq)
        and seq[i - back - 1] + seq[j + back + 1] in PAIRS
    ):
        back += 1
    fwd = 0
    while (
        fwd < cap
        and i + fwd + 1 < j - fwd - 1
        and seq[i + fwd + 1] + seq[j - fwd - 1] in PAIRS
    ):
        fwd += 1
    return 1 + back + fwd


def length_band(L: int) -> str:
    for upper, name in ((100, "L<100"), (200, "100-199"), (400, "200-399"), (10**9, "L>=400")):
        if L < upper:
            return name
    return "L>=400"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump-dir", required=True)
    ap.add_argument("--names", nargs="*", default=[])
    ap.add_argument("--split", required=True, help="jsonl with id -> seq")
    ap.add_argument("--threshold", type=float, default=0.5)
    ap.add_argument("--min-helix", type=int, default=4)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    seqs = {r["id"]: r["seq"] for r in read_jsonl(ROOT / args.split)}
    dump_dir = Path(args.dump_dir)
    files = (
        [dump_dir / f"{n}.jsonl" for n in args.names]
        if args.names
        else sorted(dump_dir.glob("*.jsonl"))
    )

    report: dict[str, dict] = {}
    for f in files:
        rows = read_jsonl(f)
        if not rows or "truth" not in rows[0] or "scored" not in rows[0]:
            continue
        stats = {
            "tp_helix": [],
            "fp_helix": [],
            "fp_long_range": 0,
            "fp_total": 0,
            "fp_long_helix_kept_link": 0,
            "by_band": defaultdict(lambda: {"fp": 0, "fn": 0, "tp": 0, "fp_long": 0, "n_seq": 0}),
        }
        for row in rows:
            seq = seqs.get(row["id"])
            if seq is None:
                continue
            truth = {tuple(p) for p in row["truth"]}
            scored = {tuple(int(x) for x in k.split(",")): float(v) for k, v in row["scored"].items()}
            pred = {k for k, v in scored.items() if v >= args.threshold}
            band = length_band(len(seq))
            b = stats["by_band"][band]
            b["n_seq"] += 1
            kept_lo, kept_hi = row["kept"]
            for i, j in pred - truth:
                h = helix_length(seq, i, j)
                stats["fp_helix"].append(h)
                stats["fp_total"] += 1
                b["fp"] += 1
                if abs(i - j) > 128:
                    stats["fp_long_range"] += 1
                    b["fp_long"] += 1
                if h >= args.min_helix and ((kept_lo <= i < kept_hi) != (kept_lo <= j < kept_hi)):
                    stats["fp_long_helix_kept_link"] += 1
            for i, j in truth & pred:
                stats["tp_helix"].append(helix_length(seq, i, j))
                b["tp"] += 1
            b["fn"] += len(truth - pred)

        fp = np.array(stats["fp_helix"]) if stats["fp_helix"] else np.array([0])
        tp = np.array(stats["tp_helix"]) if stats["tp_helix"] else np.array([0])
        report[f.stem] = {
            "n_sequences": len(rows),
            "n_fp": int(stats["fp_total"]),
            "fp_frac_helix_ge_min": round(float((fp >= args.min_helix).mean()), 4),
            "tp_frac_helix_ge_min": round(float((tp >= args.min_helix).mean()), 4),
            "fp_median_helix": float(np.median(fp)),
            "tp_median_helix": float(np.median(tp)),
            "fp_long_range_frac": round(stats["fp_long_range"] / max(1, stats["fp_total"]), 4),
            "fp_helix_spanning_kept_boundary": stats["fp_long_helix_kept_link"],
            "by_length_band": {
                k: {
                    "n_seq": v["n_seq"],
                    "fp": v["fp"],
                    "tp": v["tp"],
                    "fn": v["fn"],
                    "fp_long_range": v["fp_long"],
                }
                for k, v in sorted(stats["by_band"].items())
            },
        }
        print(json.dumps({f.stem: report[f.stem]}), flush=True)

    out = ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
