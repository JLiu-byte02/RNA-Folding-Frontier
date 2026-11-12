#!/usr/bin/env python3
"""Aggregate corrected P2-central per-seed evaluations."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--files", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    merged: dict[str, list[dict]] = {}
    for name in args.files:
        data = json.loads((ROOT / name).read_text(encoding="utf-8"))
        for key, block in data.items():
            merged.setdefault(key, []).extend(block["per_seed"])

    result = {}
    for key, seeds in sorted(merged.items()):
        block = {"n_seeds": len(seeds), "metrics": {}}
        for mode in ("fixed", "oracle"):
            for subset in ("kept_kept", "all"):
                for metric in ("family_macro_f1", "micro_f1"):
                    vals = [
                        s[mode][subset][metric]
                        for s in seeds
                        if s[mode][subset][metric] is not None
                    ]
                    block["metrics"][f"{mode}_{subset}_{metric}"] = (
                        round(float(np.mean(vals)), 6) if vals else None
                    )
        block["oracle_cuts"] = [s["oracle"]["cut"] for s in seeds]
        block["per_seed"] = seeds
        result[key] = block

    out = ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"out": args.out, "blocks": len(result)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
