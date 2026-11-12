#!/usr/bin/env python3
"""Validate and pin an energy-calibration fit before H2.4 scoring."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fit", default="experiments/002_e1prime_runs/energy_fit_coarse.json")
    ap.add_argument("--lockfile", default="state/amendments.lock.json")
    ap.add_argument("--write-lock", action="store_true")
    args = ap.parse_args()

    fit_path = ROOT / args.fit
    raw = fit_path.read_bytes()
    fit = json.loads(raw)
    selected = fit.get("selected") or {}
    s = float(selected.get("s", float("nan")))
    b = float(selected.get("b", float("nan")))
    if not (math.isfinite(s) and s > 0 and math.isfinite(b)):
        raise SystemExit(f"invalid selected calibration: {selected}")
    if fit.get("n_sequences", 0) < 1000:
        raise SystemExit(f"fit used too few sequences: {fit.get('n_sequences')}")
    if "inputs" not in fit:
        raise SystemExit("fit has no input hashes")

    sha = hashlib.sha256(raw).hexdigest()
    result = {
        "path": args.fit,
        "bytes": len(raw),
        "sha256": sha,
        "selected": {"s": s, "b": b},
        "n_sequences": fit.get("n_sequences"),
        "grid": fit.get("grid"),
        "inputs": fit.get("inputs"),
    }
    if args.write_lock:
        lock_path = ROOT / args.lockfile
        lock = json.loads(lock_path.read_text(encoding="utf-8"))
        lock.setdefault("files", {})[args.fit] = {
            "blob_sha256": sha,
            "first_quoted_in": "pending-commit",
            "note": (
                "Frozen (s,b) energy calibration from the training split; "
                "H2.4 must not run until this file is committed."
            ),
        }
        lock_path.write_text(json.dumps(lock, indent=2) + "\n", encoding="utf-8")
        result["lockfile_updated"] = args.lockfile
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
