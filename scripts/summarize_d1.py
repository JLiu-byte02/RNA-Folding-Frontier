#!/usr/bin/env python3
"""Summarise D1 JSON files into seed-level and mean tables."""

from __future__ import annotations

import argparse
import json
import re
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
NAME = re.compile(
    r"d1_seed(?P<seed>\d+)_(?P<split>id_test|test_withheld)_"
    r"(?P<policy>v1all|v1nonempty|v3all|v3primary|.+)\.json$"
)
CONDITIONS = {
    "a": ("a_id_test", "a_test_withheld"),
    "e": ("e_id_test", "e_test_withheld"),
}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default="experiments/002_e1prime_runs")
    ap.add_argument("--out", default="experiments/002_e1prime_runs/d1_summary.json")
    args = ap.parse_args()

    rows: dict[tuple, list[dict]] = defaultdict(list)
    for path in sorted((ROOT / args.dir).glob("d1_seed*.json")):
        match = NAME.search(path.name)
        if not match:
            continue
        meta = match.groupdict()
        report = json.loads(path.read_text(encoding="utf-8"))
        for cond, names in CONDITIONS.items():
            dump_name = names[0] if meta["split"] == "id_test" else names[1]
            if dump_name not in report:
                continue
            summary = report[dump_name]
            for arm, armsum in summary["arms"].items():
                band = armsum["bands"].get("80-299", {})
                oracle = band.get("at_oracle", {})
                fixed = band.get("at_fixed", {})
                rows[
                    (meta["policy"], meta["split"], cond, arm)
                ].append(
                    {
                        "seed": int(meta["seed"]),
                        "n": band.get("n"),
                        "fixed_cut": fixed.get("cut"),
                        "fixed_f1": fixed.get("family_macro_f1"),
                        "oracle_cut": oracle.get("cut"),
                        "oracle_f1": oracle.get("family_macro_f1"),
                        "oracle_boundary": oracle.get("boundary"),
                    }
                )

    table: dict[str, dict] = {}
    for key, values in sorted(rows.items()):
        policy, split, cond, arm = key
        by_seed = {str(v["seed"]): v for v in sorted(values, key=lambda x: x["seed"])}
        f1s = [v["oracle_f1"] for v in values if v["oracle_f1"] is not None]
        table[f"{policy}|{split}|{cond}|{arm}"] = {
            "seeds": by_seed,
            "mean_oracle_f1": round(sum(f1s) / len(f1s), 6) if f1s else None,
            "n_seeds": len(values),
        }

    deltas: dict[str, dict] = {}
    for policy, split in sorted({(k[0], k[1]) for k in rows}):
        for arm in ("bare", "partner", "nested"):
            a = table.get(f"{policy}|{split}|a|{arm}", {}).get("mean_oracle_f1")
            e = table.get(f"{policy}|{split}|e|{arm}", {}).get("mean_oracle_f1")
            deltas[f"{policy}|{split}|{arm}"] = (
                round(e - a, 6) if a is not None and e is not None else None
            )

    out = {
        "description": "D1 family-macro F1 (80-299 nt); oracle cut per condition/decoder/band",
        "table": table,
        "mean_delta_e_minus_a": deltas,
    }
    path = ROOT / args.out
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"out": str(path), "rows": len(table), "deltas": len(deltas)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
