#!/usr/bin/env python3
"""Cross-check v2/v3 against bpRNA-1m labels on exact-match sequences."""

from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from build_labels_v2 import iter_seeds  # noqa: E402
from rff.analysis.labels import pair_prf  # noqa: E402
from rff.data.io import normalize_seq, pairs_from_dbn  # noqa: E402
from rff.data.wuss import (  # noqa: E402
    canonical_pairs,
    load_labels_v2_delta,
    load_labels_v3,
)


def canonical(seq: str, pairs) -> set[tuple[int, int]]:
    try:
        return canonical_pairs(seq, list(pairs))
    except ValueError:
        return set()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", default="data/raw/rfam/Rfam.seed.gz")
    ap.add_argument("--bprna", default="data/raw/bprna/data.json")
    ap.add_argument("--v3", default="experiments/011_labels_v3/labels_v3_pairs.json.gz")
    ap.add_argument("--v2", default="experiments/008_labels_v2/labels_v2_delta.json.gz")
    ap.add_argument("--out", default="experiments/011_labels_v3/v2_v3_bprna_crosscheck.json")
    ap.add_argument("--max-inspect", type=int, default=10)
    args = ap.parse_args()

    v3_map = load_labels_v3(ROOT / args.v3)
    v2_delta = load_labels_v2_delta(ROOT / args.v2)
    seed_by_seq: dict[str, list[dict]] = defaultdict(list)
    for acc, name, seq, dbn_v1, _dbn_v2 in iter_seeds(ROOT / args.seed, 0):
        key = f"{acc}|{name}"
        v3 = v3_map.get(key)
        if v3 is None:
            continue
        try:
            v1 = canonical(seq, pairs_from_dbn(dbn_v1))
        except ValueError:
            v1 = set()
        v2 = v2_delta.get(key, v1)
        seed_by_seq[seq].append({"key": key, "v1": v1, "v2": v2, "v3": v3})

    bprna = json.loads((ROOT / args.bprna).read_text(encoding="utf-8"))
    comparisons = []
    inspected = []
    for rid, rec in bprna.items():
        seq = normalize_seq(rec.get("sequence", ""))
        if seq not in seed_by_seq:
            continue
        truth = canonical(
            seq, {(int(i), int(j)) for i, j in rec.get("structure", []) if i != j}
        )
        for seed in seed_by_seq[seq]:
            v2_f1 = pair_prf(truth, seed["v2"])["f1"]
            v3_f1 = pair_prf(truth, seed["v3"])["f1"]
            comparisons.append(
                {
                    "bprna_id": rid,
                    "seed_key": seed["key"],
                    "n_truth": len(truth),
                    "n_v2": len(seed["v2"]),
                    "n_v3": len(seed["v3"]),
                    "v2_f1": v2_f1,
                    "v3_f1": v3_f1,
                }
            )
            if (
                seed["v2"] != seed["v3"]
                and truth
                and len(inspected) < args.max_inspect
            ):
                inspected.append(
                    {
                        "bprna_id": rid,
                        "seed_key": seed["key"],
                        "n_truth": len(truth),
                        "n_v2": len(seed["v2"]),
                        "n_v3": len(seed["v3"]),
                        "v2_f1": round(v2_f1, 4),
                        "v3_f1": round(v3_f1, 4),
                        "v2_only": sorted(seed["v2"] - seed["v3"])[:20],
                        "v3_only": sorted(seed["v3"] - seed["v2"])[:20],
                    }
                )

    if not comparisons:
        raise SystemExit("no exact sequence overlap found")
    result = {
        "description": "v2/v3 seed labels against bpRNA-1m labels on exact sequence matches",
        "n_comparisons": len(comparisons),
        "n_unique_bprna_ids": len({c["bprna_id"] for c in comparisons}),
        "n_unique_seed_keys": len({c["seed_key"] for c in comparisons}),
        "mean_v2_f1": round(sum(c["v2_f1"] for c in comparisons) / len(comparisons), 6),
        "mean_v3_f1": round(sum(c["v3_f1"] for c in comparisons) / len(comparisons), 6),
        "v3_better_count": sum(c["v3_f1"] > c["v2_f1"] for c in comparisons),
        "v2_better_count": sum(c["v2_f1"] > c["v3_f1"] for c in comparisons),
        "inspect_v2_v3_differences": inspected,
    }
    out = ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({k: result[k] for k in result if k != "inspect_v2_v3_differences"}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
