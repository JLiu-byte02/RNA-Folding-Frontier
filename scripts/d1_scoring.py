#!/usr/bin/env python3
"""D1: did the condition change quality, or only the operating point?

K1 in Amendment 2 needs a specific measurement, so this script produces exactly
that and nothing that could be confused with it:

* three decoders — bare threshold, threshold + one partner, threshold + one
  partner + non-crossing (the production decoder, identical to
  ``pair_map_to_pairs``);
* the pre-registered metric — family-macro F1 (per-sequence F1, averaged inside
  each family, then averaged over families), not pooled pairs;
* an oracle threshold per (condition, decoder, length band), swept over a fixed
  grid so the answer cannot depend on where a quantile happened to land.

The oracle is the point of the exercise: the head was trained with
``pos_weight=20``, so a fixed 0.5 cut is an arbitrary operating point and a
condition can "win" by shifting the score scale. If the e-a gap in the 80-299
band survives each condition choosing its own best threshold, it is not
calibration.

Also reported, at the fixed cut with the production decoder: per-region and
per-distance error splits (D2), predicted/true pair ratio, score quantiles.

Outputs
-------
``<out>``                 summary JSON, one entry per dump file
``<out>.per_sequence/``   one JSONL per dump file, for the paired bootstrap
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402

from rff.analysis.labels import pair_prf  # noqa: E402
from rff.analysis.scoring import (  # noqa: E402
    Region,
    distance_bin,
    increment_sweep,
    length_band,
    prf,
)

ARMS = {
    "bare": dict(partner=False, nesting=False),
    "partner": dict(partner=True, nesting=False),
    "nested": dict(partner=True, nesting=True),
}
BANDS = ("all", "80-149", "150-299", "80-299", "L>=300")
GRID = tuple(round(0.20 + 0.05 * k, 2) for k in range(15))  # 0.20 .. 0.90


def read_jsonl(path: Path) -> list[dict]:
    with open(path) as fh:
        return [json.loads(line) for line in fh if line.strip()]


def family_macro(f1_by_row: list[float], families: list[str]) -> float:
    per_family: dict[str, list[float]] = defaultdict(list)
    for f1, fam in zip(f1_by_row, families):
        per_family[fam].append(f1)
    if not per_family:
        return float("nan")
    return float(np.mean([np.mean(v) for v in per_family.values()]))


def restrict(rows: list[dict], band: str) -> list[int]:
    if band == "all":
        return list(range(len(rows)))
    if band == "80-299":
        return [k for k, r in enumerate(rows) if 80 <= r["length"] < 300]
    return [k for k, r in enumerate(rows) if length_band(r["length"]) == band]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump-dir", required=True)
    ap.add_argument("--names", nargs="*", default=[])
    ap.add_argument("--threshold", type=float, default=0.5)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    dump_dir = Path(args.dump_dir)
    files = (
        [dump_dir / f"{n}.jsonl" for n in args.names]
        if args.names
        else sorted(dump_dir.glob("*.jsonl"))
    )
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    per_seq_dir = out.with_suffix(out.suffix + ".per_sequence")
    per_seq_dir.mkdir(parents=True, exist_ok=True)
    fixed = min(GRID, key=lambda c: abs(c - args.threshold))

    report: dict[str, dict] = {}
    for f in files:
        rows = read_jsonl(f)
        if not rows or "truth" not in rows[0] or "scored" not in rows[0]:
            continue
        for row in rows:
            row["_scored"] = row["scored"]
            row["_truth"] = {tuple(p) for p in row["truth"]}
        families = [r["family"] for r in rows]

        summary: dict = {
            "n_sequences": len(rows),
            "threshold_grid": list(GRID),
            "fixed_threshold": fixed,
            "arms": {},
        }
        f1_cache: dict[str, list[list[float]]] = {}
        oracle: dict[str, dict[str, float]] = {}

        for arm, flags in ARMS.items():
            decodes = increment_sweep(rows, list(GRID), **flags)
            f1 = [
                [pair_prf(r["_truth"], decodes[k][c])["f1"] for k, r in enumerate(rows)]
                for c in range(len(GRID))
            ]
            f1_cache[arm] = f1
            oracle[arm] = {}
            armsum: dict = {"bands": {}}
            for band in BANDS:
                idx = restrict(rows, band)
                if not idx:
                    armsum["bands"][band] = {"n": 0}
                    oracle[arm][band] = 0.5
                    continue
                fams = [families[k] for k in idx]
                curve = [
                    family_macro([f1[c][k] for k in idx], fams) for c in range(len(GRID))
                ]
                best = int(np.argmax(curve))
                oracle[arm][band] = GRID[best]
                armsum["bands"][band] = {
                    "n": len(idx),
                    "at_fixed": {
                        "cut": fixed,
                        "family_macro_f1": round(float(curve[GRID.index(fixed)]), 4),
                    },
                    "at_oracle": {
                        "cut": GRID[best],
                        "family_macro_f1": round(float(curve[best]), 4),
                        "curve": {str(c): round(float(v), 4) for c, v in zip(GRID, curve)},
                        "boundary": best in (0, len(GRID) - 1),
                    },
                }
            summary["arms"][arm] = armsum

        with open(per_seq_dir / f"{f.stem}.jsonl", "w") as fh:
            for k, row in enumerate(rows):
                rec = {
                    "id": row["id"],
                    "family": row["family"],
                    "source": row.get("source"),
                    "length": row["length"],
                    "band": length_band(row["length"]),
                }
                for arm in ARMS:
                    f1 = f1_cache[arm]
                    rec[f"f1_{arm}_fixed"] = round(f1[GRID.index(fixed)][k], 6)
                    rec[f"f1_{arm}_oracle_all"] = round(
                        f1[GRID.index(oracle[arm]["all"])][k], 6
                    )
                    rec[f"f1_{arm}_oracle_80_299"] = (
                        round(f1[GRID.index(oracle[arm]["80-299"])][k], 6)
                        if 80 <= row["length"] < 300
                        else None
                    )
                fh.write(json.dumps(rec) + "\n")

        nested0 = increment_sweep(rows, [fixed], **ARMS["nested"])
        tot = {"tp": 0, "fp": 0, "fn": 0}
        by_region: dict[str, dict[str, int]] = {}
        by_dist: dict[str, dict[str, int]] = {}
        preds, trues = [], []
        for k, row in enumerate(rows):
            pred = set(nested0[k][0])
            preds.append(len(pred))
            trues.append(len(row["_truth"]))
            reg = Region(row["kept"][0], row["kept"][1])
            for name, s in (
                ("tp", row["_truth"] & pred),
                ("fp", pred - row["_truth"]),
                ("fn", row["_truth"] - pred),
            ):
                tot[name] += len(s)
                for i, j in s:
                    by_region.setdefault(reg.where(i, j), {"tp": 0, "fp": 0, "fn": 0})[
                        name
                    ] += 1
                    by_dist.setdefault(distance_bin(i, j), {"tp": 0, "fp": 0, "fn": 0})[
                        name
                    ] += 1
        summary["error_split_fixed_nested"] = {
            **{k: int(v) for k, v in tot.items()},
            **{k: round(v, 4) if isinstance(v, float) else v for k, v in prf(**tot).items()},
            "mean_predicted_pairs": round(float(np.mean(preds)), 2),
            "mean_true_pairs": round(float(np.mean(trues)), 2),
            "pred_to_true_ratio": round(
                float(np.mean(preds) / max(1e-9, np.mean(trues))), 3
            ),
            "by_region": {
                k: {
                    kk: (round(vv, 4) if isinstance(vv, float) else vv)
                    for kk, vv in prf(**v).items()
                }
                for k, v in sorted(by_region.items())
            },
            "by_distance": {
                k: {
                    kk: (round(vv, 4) if isinstance(vv, float) else vv)
                    for kk, vv in prf(**v).items()
                }
                for k, v in sorted(by_dist.items())
            },
        }

        scores = np.concatenate(
            [
                np.fromiter((float(v) for v in r["scored"].values()), dtype=float)
                for r in rows
            ]
        )
        summary["score_distribution"] = {
            "n_cells": int(scores.size),
            "quantiles": {
                str(q): round(float(np.quantile(scores, q)), 4)
                for q in (0.5, 0.9, 0.99, 0.999)
            },
            "frac_ge_fixed": round(float((scores >= fixed).mean()), 4),
        }

        report[f.stem] = summary
        print(
            json.dumps(
                {
                    f.stem: summary["arms"]["partner"]["bands"]["80-299"]["at_oracle"][
                        "family_macro_f1"
                    ]
                }
            )
        )

    out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"out": str(out), "per_sequence": str(per_seq_dir)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
