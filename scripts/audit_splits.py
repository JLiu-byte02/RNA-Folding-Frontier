#!/usr/bin/env python3
"""Run the leakage audit from a JSON manifest.

Manifest format:

{
  "reference_split": "train",
  "splits": {
    "train": {"path": "data/...", "family_from": "filename"},
    "test":  {"path": "data/...", "family_from": "filename"}
  }
}
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from rff.data.audit import run_audit, save_report  # noqa: E402
from rff.data.io import load_path  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest")
    parser.add_argument("--out", default="experiments/000_e0_leakage/report")
    args = parser.parse_args()

    manifest = json.loads(Path(args.manifest).read_text())
    splits = {}
    for name, spec in manifest["splits"].items():
        records = load_path(spec["path"])
        family_from = spec.get("family_from")
        if family_from == "bprna":
            for rec in records:
                match = re.search(r"RFAM_(\d+)", rec.rid)
                rec.family = f"RFAM_{match.group(1)}" if match else "unknown"
        elif family_from == "prefix":
            for rec in records:
                rec.family = Path(rec.rid).name.split("_")[0]
        elif family_from == "filename":
            for rec in records:
                if rec.family is None:
                    rec.family = Path(rec.rid).name.split("_")[0]
        for rec in records:
            rec.split = name
        splits[name] = records
        print(f"{name}: {len(records)} records")

    report = run_audit(splits, manifest.get("reference_split", "train"))
    save_report(report, args.out)
    print(f"wrote {args.out}/audit.json and audit.md")


if __name__ == "__main__":
    main()

