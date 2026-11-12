#!/usr/bin/env python3
"""K1: is "condition e helps" just a calibration shift?

Amendment 2 A2.4 fixes the test: oracle threshold plus one-partner decoding,
condition e minus condition a, family-macro F1 restricted to the 80-299 nt band.
If that difference is below 0.01 with a bootstrap interval that includes 0, the
"distal context is a net-negative prior" reading is dropped and H2 is not
pursued.

Inputs are the per-sequence files from ``d1_scoring.py``, one per seed per
condition. Seeds are averaged within each sequence before the family means are
taken, so a sequence contributes one number, not one per seed. The bootstrap is
paired over families, matching the P2 table.

All three decoders are printed, because the contrast between them is itself
informative: if only the bare-threshold arm shows a gain, the effect is density,
not quality.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from analyze_p2 import bootstrap_delta  # noqa: E402

ARMS = ("bare", "partner", "nested")


def read_jsonl(path: Path) -> list[dict]:
    with open(path) as fh:
        return [json.loads(line) for line in fh if line.strip()]


def load(files: list[str], arm: str, band: tuple[int, int]) -> dict[str, float]:
    per_id: dict[str, list[float]] = defaultdict(list)
    family: dict[str, str] = {}
    for f in files:
        for row in read_jsonl(Path(f)):
            if not (band[0] <= row["length"] < band[1]):
                continue
            value = row.get(f"f1_{arm}_oracle_80_299")
            if value is None:
                continue
            per_id[row["id"]].append(float(value))
            family[row["id"]] = row["family"]
    by_family: dict[str, list[float]] = defaultdict(list)
    for rid, values in per_id.items():
        by_family[family[rid]].append(float(np.mean(values)))
    return {f: float(np.mean(v)) for f, v in by_family.items()}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--a", nargs="+", required=True, help="per-sequence jsonl, one per seed")
    ap.add_argument("--e", nargs="+", required=True)
    ap.add_argument("--metric-label", default="oracle threshold, 80-299 nt band")
    ap.add_argument("--n-boot", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--kill-delta", type=float, default=0.01)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    result = {"metric": args.metric_label, "band": [80, 300], "arms": {}}
    for arm in ARMS:
        base = load(args.a, arm, (80, 300))
        other = load(args.e, arm, (80, 300))
        stats = bootstrap_delta(base, other, args.n_boot, args.seed)
        result["arms"][arm] = stats
        print(json.dumps({arm: stats}))

    partner = result["arms"]["partner"]
    fires = (
        partner["delta"] is not None
        and partner["delta"] < args.kill_delta
        and partner["lo"] is not None
        and partner["lo"] <= 0.0 <= partner["hi"]
    )
    result["k1_fires"] = bool(fires)
    result["rule"] = (
        f"drop the 'distal context is a net-negative prior' reading iff "
        f"partner-arm e-a < {args.kill_delta} and the 95% CI includes 0"
    )
    if partner["delta"] is None:
        result["verdict"] = "insufficient families"
    elif partner["delta"] < 0.0 and partner["hi"] is not None and partner["hi"] < 0.0:
        result["verdict"] = "e is worse than a; the gain premise is not supported"
    elif fires:
        result["verdict"] = "K1 FIRES - calibration is a sufficient explanation"
    elif partner["delta"] <= 0.0:
        result["verdict"] = "no evidence that e helps"
    else:
        result["verdict"] = "K1 does not fire - the gain survives an oracle cut"
    print(json.dumps({k: result[k] for k in ("k1_fires", "verdict")}))
    if args.out:
        out = Path(args.out)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(result, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
