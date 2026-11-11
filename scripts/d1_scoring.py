#!/usr/bin/env python3
"""D1/D2 diagnostics on the dumped P2 predictions.

D1 answers the question the P2 table cannot: did the condition change the amount
of predicted structure (a calibration/recall shift) or the quality of it at a
fixed operating point? We report predicted-pair counts against true counts,
precision/recall at 0.5, AUPRC, and an oracle threshold per condition.

D2 splits the same errors by where the pair sits relative to the kept window and
by pair distance, and measures how often a false positive has a compensating
complementary run at its ends (the "spurious partner" hypothesis).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402

from rff.analysis.scoring import (  # noqa: E402
    Region,
    average_precision,
    best_threshold,
    distance_bin,
    pair_scores,
    prf,
)


def read_jsonl(path: Path) -> list[dict]:
    with open(path) as fh:
        return [json.loads(line) for line in fh if line.strip()]


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
    report: dict[str, dict] = {}
    for f in files:
        rows = read_jsonl(f)
        if not rows or "truth" not in rows[0] or "scored" not in rows[0]:
            # Not a prediction dump (e.g. a split file sitting in the same dir).
            continue
        agg = {
            "at_threshold": args.threshold,
            "n_sequences": len(rows),
        }
        tot = {"tp": 0, "fp": 0, "fn": 0}
        by_region: dict[str, dict[str, int]] = {}
        by_dist: dict[str, dict[str, int]] = {}
        ys, ss, preds, trues = [], [], [], []
        for row in rows:
            truth = {tuple(p) for p in row["truth"]}
            scored = {tuple(int(x) for x in k.split(",")): float(v) for k, v in row["scored"].items()}
            pred = {k for k, v in scored.items() if v >= args.threshold}
            preds.append(len(pred))
            trues.append(len(truth))
            reg = Region(row["kept"][0], row["kept"][1])
            for name, s in (("tp", truth & pred), ("fp", pred - truth), ("fn", truth - pred)):
                tot[name] += len(s)
                for i, j in s:
                    by_region.setdefault(reg.where(i, j), {"tp": 0, "fp": 0, "fn": 0})[name] += 1
                    by_dist.setdefault(distance_bin(i, j), {"tp": 0, "fp": 0, "fn": 0})[name] += 1
            y, s = pair_scores(truth, scored)
            ys.append(y)
            ss.append(s)

        y_all, s_all = np.concatenate(ys), np.concatenate(ss)
        cut, oracle = best_threshold(y_all, s_all)
        agg.update({k: int(v) for k, v in tot.items()})
        agg.update({k: round(v, 4) if isinstance(v, float) else v for k, v in prf(**tot).items()})
        agg["mean_predicted_pairs"] = round(float(np.mean(preds)), 2)
        agg["mean_true_pairs"] = round(float(np.mean(trues)), 2)
        agg["pred_to_true_ratio"] = round(float(np.mean(preds) / max(1e-9, np.mean(trues))), 3)
        agg["auprc"] = round(average_precision(y_all, s_all), 4)
        agg["oracle_threshold"] = round(cut, 4)
        agg["oracle_f1"] = round(oracle, 4)
        agg["by_region"] = {
            k: {kk: (round(vv, 4) if isinstance(vv, float) else vv) for kk, vv in prf(**v).items()}
            for k, v in sorted(by_region.items())
        }
        agg["by_distance"] = {
            k: {kk: (round(vv, 4) if isinstance(vv, float) else vv) for kk, vv in prf(**v).items()}
            for k, v in sorted(by_dist.items())
        }
        report[f.stem] = agg
        print(json.dumps({f.stem: agg}), flush=True)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
