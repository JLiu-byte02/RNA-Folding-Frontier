#!/usr/bin/env python3
"""Exploratory H2.4 check: are FP helices energetically competitive?"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rff.data.io import pairs_from_dbn  # noqa: E402

try:
    import RNA
except ImportError as exc:  # pragma: no cover
    raise SystemExit("ViennaRNA's Python module `RNA` is required") from exc

CANONICAL = set("AU UA GC CG GU UG".split())
LONG_RANGE = 128


def canonical(seq: str, i: int, j: int) -> bool:
    return 0 <= i < j < len(seq) and seq[i] + seq[j] in CANONICAL


def stack_pairs(seq: str, i: int, j: int, cap: int = 40) -> list[tuple[int, int]]:
    if not canonical(seq, i, j):
        return []
    pairs = {(i, j)}
    a, b = i, j
    while (
        len(pairs) < cap
        and a - 1 >= 0
        and b + 1 < len(seq)
        and canonical(seq, a - 1, b + 1)
    ):
        a -= 1
        b += 1
        pairs.add((a, b))
    a, b = i, j
    while (
        len(pairs) < cap
        and a + 1 < b - 1
        and canonical(seq, a + 1, b - 1)
    ):
        a += 1
        b -= 1
        pairs.add((a, b))
    return sorted(pairs)


def isolated_energy(seq: str, pairs: list[tuple[int, int]]) -> float:
    lo = min(i for i, _ in pairs)
    hi = max(j for _, j in pairs)
    sub = seq[lo : hi + 1]
    db = ["."] * len(sub)
    for i, j in pairs:
        db[i - lo] = "("
        db[j - lo] = ")"
    return float(RNA.energy_of_struct(sub, "".join(db)))


def quantiles(values: list[float]) -> dict:
    if not values:
        return {"n": 0}
    arr = np.asarray(values, dtype=float)
    return {
        "n": int(arr.size),
        "mean": round(float(arr.mean()), 4),
        "median": round(float(np.median(arr)), 4),
        "p10": round(float(np.percentile(arr, 10)), 4),
        "p90": round(float(np.percentile(arr, 90)), 4),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump-dir", required=True)
    ap.add_argument("--split", required=True)
    ap.add_argument("--name", required=True)
    ap.add_argument("--threshold", type=float, default=0.5)
    ap.add_argument("--min-true-pairs", type=int, default=4)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    seqs = {}
    with open(ROOT / args.split, encoding="utf-8") as fh:
        for line in fh:
            if line.strip():
                rec = json.loads(line)
                seqs[rec["id"]] = rec["seq"]

    fp_energies: list[float] = []
    true_energies: list[float] = []
    fp_in_mfe = 0
    fp_any_in_mfe = 0
    n_fp_helices = 0
    n_rows = 0
    n_skipped = 0
    with open(Path(args.dump_dir) / f"{args.name}.jsonl", encoding="utf-8") as fh:
        for line in fh:
            if not line.strip():
                continue
            row = json.loads(line)
            seq = seqs.get(row["id"])
            if seq is None:
                n_skipped += 1
                continue
            truth_raw = row.get("truth_v3")
            if truth_raw is None or len(truth_raw) < args.min_true_pairs:
                continue
            n_rows += 1
            truth = {tuple(p) for p in truth_raw}
            pred = {
                tuple(int(x) for x in key.split(","))
                for key, value in row["scored"].items()
                if float(value) >= args.threshold
            }
            mfe_struct, _ = RNA.fold(seq)
            mfe_pairs = set(pairs_from_dbn(mfe_struct))

            fp_long = {p for p in pred - truth if abs(p[0] - p[1]) > LONG_RANGE}
            fp_helices = {
                tuple(stack)
                for p in fp_long
                for stack in [stack_pairs(seq, *p)]
                if len(stack) >= 4
            }
            true_helices = {
                tuple(stack)
                for p in truth
                if abs(p[0] - p[1]) > LONG_RANGE
                for stack in [stack_pairs(seq, *p)]
                if len(stack) >= 4
            }
            for helix in fp_helices:
                n_fp_helices += 1
                fp_energies.append(isolated_energy(seq, list(helix)))
                hp = set(helix)
                fp_in_mfe += hp.issubset(mfe_pairs)
                fp_any_in_mfe += bool(hp & mfe_pairs)
            for helix in true_helices:
                true_energies.append(isolated_energy(seq, list(helix)))

    result = {
        "name": args.name,
        "split": args.split,
        "n_rows_scored": n_rows,
        "n_sequences_missing": n_skipped,
        "long_range_definition": f"|i-j| > {LONG_RANGE}",
        "helix_definition": "contiguous canonical stack of >=4 pairs",
        "fp_helices": quantiles(fp_energies),
        "true_helices": quantiles(true_energies),
        "n_fp_helices": n_fp_helices,
        "fp_helix_fully_in_rnafold_mfe": fp_in_mfe,
        "fp_helix_fully_in_mfe_fraction": round(fp_in_mfe / n_fp_helices, 4)
        if n_fp_helices
        else None,
        "fp_helix_any_pair_in_mfe": fp_any_in_mfe,
        "fp_helix_any_pair_in_mfe_fraction": round(fp_any_in_mfe / n_fp_helices, 4)
        if n_fp_helices
        else None,
    }
    out = ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"out": args.out, **result}, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
