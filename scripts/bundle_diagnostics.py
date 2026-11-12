#!/usr/bin/env python3
"""Bundle D1/K1/D2 summaries with per-file hashes for auditability."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PATTERNS = (
    "d1_seed*.json",
    "k1_*.json",
    "d2_seed17_*.json",
    "d2c_seed17_*.json",
    "d1_summary.json",
)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default="experiments/002_e1prime_runs")
    ap.add_argument(
        "--out", default="experiments/002_e1prime_runs/diagnostics_bundle_2026-11-12.json"
    )
    args = ap.parse_args()

    source = ROOT / args.dir
    files: dict[str, dict] = {}
    for pattern in PATTERNS:
        for path in sorted(source.glob(pattern)):
            rel = path.relative_to(ROOT).as_posix()
            raw = path.read_bytes()
            files[rel] = {
                "sha256": hashlib.sha256(raw).hexdigest(),
                "bytes": len(raw),
                "content": json.loads(raw),
            }
    if not files:
        raise SystemExit(f"no result files matched in {source}")

    bundle = {
        "description": (
            "Audit bundle for the D1/K1/D2 diagnostics run on 2026-11-12. "
            "Each entry carries the exact content and sha256 of its source file."
        ),
        "n_files": len(files),
        "files": files,
    }
    out = ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(bundle, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "out": out.relative_to(ROOT).as_posix(),
                "n_files": len(files),
                "bytes": out.stat().st_size,
                "sha256": hashlib.sha256(out.read_bytes()).hexdigest(),
            }
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
