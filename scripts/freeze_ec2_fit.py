#!/usr/bin/env python3
"""Validate and pin the E-C2 coverage-gate fit before any withheld scoring.

Amendment 13 A13.2 requires the gate parameters to be frozen and hashed before
the withheld split is scored.  This script refuses obviously invalid or
degenerate fits and (optionally) appends the digest to the lock file.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

GRID = {
    "s_max": [0.1875, 0.5, 1.0],
    "alpha": [-1.0, 0.0, 1.0],
    "beta": [0.0, 1.0, 2.0],
    "gamma": [0.0, 0.25, 0.5],
}


def _edge(value: float, axis: str) -> bool:
    vals = GRID[axis]
    return abs(value - min(vals)) < 1e-12 or abs(value - max(vals)) < 1e-12


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fit", default="experiments/002_e1prime_runs/ec2_fit.json")
    ap.add_argument("--lockfile", default="state/amendments.lock.json")
    ap.add_argument("--write-lock", action="store_true")
    args = ap.parse_args()

    fit_path = ROOT / args.fit
    raw = fit_path.read_bytes()
    fit = json.loads(raw)

    if fit.get("amendment") != "AMENDMENT_13.md":
        raise SystemExit(f"unexpected amendment tag: {fit.get('amendment')!r}")
    n = int(fit.get("n_sequences", 0))
    if n < 1000:
        raise SystemExit(f"fit used too few sequences: {n}")
    if not fit.get("head_source_sha256") or not fit.get("pppl_source_sha256"):
        raise SystemExit("fit is missing head/pppl source hashes")
    if fit.get("n_dropped", {}).get("no_pppl"):
        raise SystemExit("some sequences had no pseudo-perplexity; fit is incomplete")

    best_g1 = fit.get("best_g1") or {}
    best_g2 = fit.get("best_g2") or {}
    for name, cfg in (("best_g1", best_g1), ("best_g2", best_g2)):
        for axis in ("s_max", "alpha", "beta"):
            if not math.isfinite(float(cfg.get(axis, float("nan")))):
                raise SystemExit(f"{name} has non-finite {axis}: {cfg.get(axis)}")
        if not math.isfinite(float(cfg.get("mean_f1", float("nan")))):
            raise SystemExit(f"{name} has non-finite mean_f1")
        if name == "best_g2" and not math.isfinite(float(cfg.get("gamma", float("nan")))):
            raise SystemExit("best_g2 has non-finite gamma")

    edges = [
        axis
        for axis in ("s_max", "alpha", "beta")
        if _edge(float(best_g1[axis]), axis)
    ]
    sha = hashlib.sha256(raw).hexdigest()
    result = {
        "path": args.fit,
        "bytes": len(raw),
        "sha256": sha,
        "n_sequences": n,
        "folds": fit.get("folds"),
        "best_g1": best_g1,
        "best_g2": best_g2,
        "head_source_sha256": fit.get("head_source_sha256"),
        "pppl_source_sha256": fit.get("pppl_source_sha256"),
        "g1_on_grid_edge": edges,
        "warning": (
            "best_g1 sits on a grid edge for " + ", ".join(edges) +
            "; the pinned value is a boundary optimum and any refinement must be "
            "pre-registered in a new amendment before withheld scoring"
        ) if edges else None,
    }
    if args.write_lock:
        lock_path = ROOT / args.lockfile
        lock = json.loads(lock_path.read_text(encoding="utf-8"))
        lock.setdefault("files", {})[args.fit] = {
            "blob_sha256": sha,
            "first_quoted_in": "pending-commit",
            "note": (
                "Frozen E-C2 coverage-gate parameters fitted on the training "
                "split only; withheld scoring must not run until this file is "
                "committed (Amendment 13 A13.2)."
            ),
        }
        lock_path.write_text(json.dumps(lock, indent=2) + "\n", encoding="utf-8")
        result["lockfile_updated"] = args.lockfile
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
